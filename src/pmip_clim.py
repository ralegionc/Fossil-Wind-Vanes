"""Monthly-mean surface wind climatologies (last 50 years) for each PMIP run.

Files are subset remotely over OPeNDAP where a node offers it. Otherwise they
are downloaded to a temporary folder on C: (the D: drive is nearly full),
reduced, and deleted. Only the 12-month climatology over 40S to 55N is kept.
"""
import json
import os
from pathlib import Path
import sys
import tempfile

import numpy as np
import requests
import xarray as xr

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "data" / "processed" / "pmip_manifest.json"
OUTDIR = ROOT / "data" / "processed" / "pmip"
TMP = Path(os.environ.get("PMIP_TMP", tempfile.gettempdir())) / "fossil_wind_vanes_pmip"
MONTHS = 600
LAT = (-40.0, 55.0)


def dap_ok(url):
    try:
        return requests.head(url + ".dds", timeout=20).status_code == 200
    except requests.RequestException:
        return False


TIMES = xr.coders.CFDatetimeCoder(use_cftime=True)


def open_remote(f, var):
    for u in f["opendap"]:
        if dap_ok(u):
            return xr.open_dataset(u, decode_times=TIMES), None
    TMP.mkdir(parents=True, exist_ok=True)
    for u in f["http"]:
        dest = TMP / f["title"]
        try:
            with requests.get(u, stream=True, timeout=120) as r:
                r.raise_for_status()
                with open(dest, "wb") as fh:
                    for chunk in r.iter_content(1 << 20):
                        fh.write(chunk)
            return xr.open_dataset(dest, decode_times=TIMES, cache=False), dest
        except requests.RequestException as err:
            print("   http failed", u.split("/")[2], err, flush=True)
    raise RuntimeError(f"no working URL for {f['title']}")


def lat_slice(ds):
    lat = ds["lat"].values
    return slice(*LAT) if lat[0] < lat[-1] else slice(LAT[1], LAT[0])


def climatology(sel, var):
    parts = []
    # newest files last; walk backwards until we have MONTHS months
    have = 0
    for f in reversed(sel["files"]):
        ds, tmp = open_remote(f, var)
        need = MONTHS - have
        x = ds[var].sel(lat=lat_slice(ds)).isel(time=slice(-need, None)).load()
        ds.close()  # release the handle before deleting (Windows locks open files)
        if tmp:
            tmp.unlink(missing_ok=True)
        parts.append(x)
        have += x.sizes["time"]
        if have >= MONTHS:
            break
    x = xr.concat(parts[::-1], dim="time").sortby("time")
    clim = x.groupby("time.month").mean("time").astype("float32")
    clim.attrs.update(n_months=int(x.sizes["time"]),
                      first=str(x.time.values[0]), last=str(x.time.values[-1]))
    return clim


def main():
    manifest = json.loads(MANIFEST.read_text())
    OUTDIR.mkdir(parents=True, exist_ok=True)
    items = list(manifest.items())
    if "--reverse" in sys.argv:  # lets a second worker start from the other end
        items = items[::-1]
    for key, sel in items:
        if not sel:
            continue
        model, exp, var = key.split("|")
        out = OUTDIR / f"{model}_{exp}_{var}.nc"
        claim = out.with_suffix(".working")
        if out.exists() or claim.exists():
            continue
        claim.touch()
        try:
            clim = climatology(sel, var)
        except Exception as err:
            print(f"{key}: FAILED {err}", flush=True)
            claim.unlink(missing_ok=True)
            continue
        clim.to_dataset(name=var).to_netcdf(out)
        claim.unlink(missing_ok=True)
        print(f"{key}: {clim.attrs['n_months']} months {clim.attrs['first'][:7]}..{clim.attrs['last'][:7]} "
              f"grid {clim.sizes['lat']}x{clim.sizes['lon']}", flush=True)


if __name__ == "__main__":
    sys.exit(main())
