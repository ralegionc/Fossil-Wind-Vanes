"""Daily surface winds from AWI-ESM-1-1-LR and MPI-ESM1-2-LR (LGM and pre-industrial)
at the model grid cells nearest each dune tile, reduced to sand flux roses.

Whole files that must be downloaded go to a temporary folder on C: one at a
time and are deleted after extraction. Extracted daily series stay on C:;
only the flux roses (about 1 MB) are written to the project on D:.
"""
from concurrent.futures import ThreadPoolExecutor
import os
from pathlib import Path
import re
import sys
import tempfile

import numpy as np
import pandas as pd
import requests
import xarray as xr

from flux import KAPPA, Z0, flux_rose

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
WORK = Path(os.environ.get("PMIP_TMP", tempfile.gettempdir())) / "fossil_wind_vanes_daily"
SEARCH = "https://esgf.ceda.ac.uk/esg-search/search"
MODELS = ["AWI-ESM-1-1-LR", "MPI-ESM1-2-LR"]
EXPS = ["lgm", "piControl"]
YEARS = 30
# DKRZ, the only server holding MPI-ESM1-2-LR daily LGM winds, slowed to ~40 kB/s on
# 2026-10-12, so that run uses its last 20 years (one file). See docs/PLAN.md.
YEARS_OVERRIDE = {("MPI-ESM1-2-LR", "lgm"): 20}
THRESHOLDS = [3.0, 4.0, 5.0]   # daily-mean 10 m wind speed thresholds (m/s)
TIMES = xr.coders.CFDatetimeCoder(use_cftime=True)


def files_for(model, exp, var, years=YEARS):
    q = dict(project="CMIP6", source_id=model, experiment_id=exp, variable_id=var, table_id="day",
             member_id="r1i1p1f1", type="File", distrib="true", limit=5000,
             format="application/solr+json", fields="title,url,dataset_id")
    docs = requests.get(SEARCH, params=q, timeout=120).json()["response"]["docs"]
    ver = max(d["dataset_id"].split("|")[0].rsplit(".", 1)[-1] for d in docs)
    files = {}
    for d in docs:
        if d["dataset_id"].split("|")[0].rsplit(".", 1)[-1] != ver:
            continue
        f = files.setdefault(d["title"], {"title": d["title"], "http": [], "opendap": []})
        for u in d["url"]:
            url, _, kind = u.split("|")
            if kind == "HTTPServer":
                f["http"].append(url)
            elif kind == "OPENDAP":
                f["opendap"].append(url.removesuffix(".html"))
    span = {t: tuple(int(x[:4]) for x in re.search(r"_(\d{8})-(\d{8})\.nc$", t).groups()) for t in files}
    end = max(s[1] for s in span.values())
    return [files[t] for t in sorted(files) if span[t][1] > end - years]


def dap_ok(url):
    try:
        return requests.head(url + ".dds", timeout=20).status_code == 200
    except requests.RequestException:
        return False


def download(urls, dest, tries=20):
    """Download with size check, resuming partial files (connections drop silently)."""
    for attempt in range(tries):
        for u in urls:
            try:
                total = int(requests.head(u, timeout=60, allow_redirects=True).headers["Content-Length"])
                have = dest.stat().st_size if dest.exists() else 0
                if have == total:
                    return
                if have > total:
                    dest.unlink()
                    have = 0
                hdr = {"Range": f"bytes={have}-"} if have else {}
                with requests.get(u, stream=True, timeout=300, headers=hdr) as r:
                    r.raise_for_status()
                    mode = "ab" if (have and r.status_code == 206) else "wb"
                    with open(dest, mode) as fh:
                        for chunk in r.iter_content(1 << 22):
                            fh.write(chunk)
                if dest.stat().st_size == total:
                    return
                print(f"   short download {dest.stat().st_size}/{total}, resuming", flush=True)
            except (requests.RequestException, KeyError) as err:
                print("   http failed", u.split("/")[2], err, flush=True)
    raise RuntimeError(f"could not download {dest.name}")


def extract(f, var, cells_lat, cells_lon):
    """Daily values (n_cells, n_days) at model cells given by index arrays."""
    src = next((u for u in f["opendap"] if dap_ok(u)), None)
    tmp = None
    if src is None:
        WORK.mkdir(parents=True, exist_ok=True)
        tmp = WORK / f["title"]
        download(f["http"], tmp)
        src = tmp
    ds = xr.open_dataset(src, decode_times=TIMES, cache=False)
    lo, hi = int(cells_lat.min()), int(cells_lat.max())
    band = ds[var].isel(lat=slice(lo, hi + 1)).load()   # one contiguous request
    times = ds["time"].values
    ds.close()
    if tmp is not None:
        tmp.unlink(missing_ok=True)
    vals = band.values[:, cells_lat - lo, cells_lon]     # (time, n_cells)
    return vals.T.astype(np.float32), times


def grid_of(model):
    f = files_for(model, "piControl", "uas")[-1]
    src = next(u for u in f["opendap"] if dap_ok(u))
    ds = xr.open_dataset(src, decode_times=TIMES)
    lat, lon = ds["lat"].values, ds["lon"].values
    ds.close()
    return lat, lon


def main(models=MODELS):
    tiles = pd.read_csv(PROC / "tiles.csv")
    WORK.mkdir(parents=True, exist_ok=True)
    for model in models:
        out = {}
        lat, lon = grid_of(model)
        ti = np.abs(tiles.lat.values[:, None] - lat[None, :]).argmin(1)
        dlon = np.abs(((tiles.lon.values[:, None] % 360) - (lon[None, :] % 360) + 180) % 360 - 180)
        tj = dlon.argmin(1)
        cells, tile_cell = np.unique(np.c_[ti, tj], axis=0, return_inverse=True)
        out[f"{model}|tile_cell"] = tile_cell.ravel()
        out[f"{model}|cell_lat"] = lat[cells[:, 0]]
        out[f"{model}|cell_lon"] = lon[cells[:, 1]]
        print(model, "grid", lat.size, "x", lon.size, "->", len(cells), "cells", flush=True)
        for exp in EXPS:
            years = YEARS_OVERRIDE.get((model, exp), YEARS)
            series = {}
            for var in ("uas", "vas"):
                cache = WORK / f"{model}_{exp}_{var}.npy"
                if cache.exists():
                    series[var] = np.load(cache)[:, -years * 365:]
                    continue
                fl = files_for(model, exp, var, years)
                print(f"  {exp} {var}: {len(fl)} files", flush=True)
                with ThreadPoolExecutor(max_workers=2) as ex:
                    parts = list(ex.map(lambda f: extract(f, var, cells[:, 0], cells[:, 1])[0], fl))
                arr = np.concatenate(parts, axis=1)[:, -years * 365:]
                np.save(cache, arr)
                series[var] = arr
                print(f"  {exp} {var}: {arr.shape[1]} days", flush=True)
            n = min(series["uas"].shape[1], series["vas"].shape[1])
            u, v = series["uas"][:, -n:].astype(float), series["vas"][:, -n:].astype(float)
            out[f"{model}|{exp}|n_days"] = n
            out[f"{model}|{exp}|mean_speed"] = np.hypot(u, v).mean(1)
            for th in THRESHOLDS:
                u_star_t = th * KAPPA / np.log(10.0 / Z0)
                out[f"{model}|{exp}|rose{th:g}"] = flux_rose(u, v, u_star_t=u_star_t)
        dest = PROC / f"daily_roses_{model}.npz"
        np.savez_compressed(dest, **out)
        print("saved", dest, flush=True)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:] or MODELS))
