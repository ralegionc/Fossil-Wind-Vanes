"""Directional spectra of young dunes from Sentinel-2 (test D in docs/PLAN.md).

For each base tile, the central 10 km of the Sentinel-2 L2A red band (10 m) is
read straight from the public AWS bucket for two clear scenes from opposite
seasons, so the sun lights the dunes from different directions. Each scene's
spectrum is saved on a (wavelength x true-north direction) grid; images are not
stored.
"""
from concurrent.futures import ThreadPoolExecutor, as_completed
import math
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import requests

from orient import power_spectrum

os.environ.setdefault("AWS_NO_SIGN_REQUEST", "YES")
os.environ.setdefault("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
os.environ.setdefault("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", ".tif")
os.environ.setdefault("GDAL_HTTP_MAX_RETRY", "4")

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
OUTDIR = PROC / "s2polar"
STAC = "https://earth-search.aws.element84.com/v1/search"
HALF_M = 5000.0
PIX = 10.0
WL_EDGES = np.logspace(np.log10(30.0), np.log10(3000.0), 41)
SEASONS = {"jj": ("06", "07"), "dj": ("12", "01")}


def search(lon, lat, half_deg_lon, half_deg_lat):
    poly = {"type": "Polygon", "coordinates": [[
        [lon - half_deg_lon, lat - half_deg_lat], [lon + half_deg_lon, lat - half_deg_lat],
        [lon + half_deg_lon, lat + half_deg_lat], [lon - half_deg_lon, lat + half_deg_lat],
        [lon - half_deg_lon, lat - half_deg_lat]]]}
    body = {"collections": ["sentinel-2-l2a"], "intersects": poly,
            "datetime": "2021-01-01T00:00:00Z/2024-12-31T23:59:59Z",
            "query": {"eo:cloud_cover": {"lt": 5}},
            "sortby": [{"field": "properties.eo:cloud_cover", "direction": "asc"}], "limit": 200}
    for attempt in range(4):
        try:
            r = requests.post(STAC, json=body, timeout=90)
            r.raise_for_status()
            return r.json()["features"], poly
        except requests.RequestException:
            if attempt == 3:
                raise
    return [], poly


def pick_scenes(items, poly):
    """Clearest scene per season whose footprint holds the whole window."""
    from shapely.geometry import shape
    win = shape(poly)
    out = {}
    for it in items:
        month = it["properties"]["datetime"][5:7]
        for season, months in SEASONS.items():
            href = it["assets"].get("red", {}).get("href", "")
            # some catalogue entries point at JPEG2000 files in a requester-pays bucket
            if not (href.startswith("https://sentinel-cogs") and href.endswith(".tif")):
                continue
            if month in months and season not in out and shape(it["geometry"]).contains(win):
                out[season] = it
    return out


def read_window(item, lon, lat):
    import rasterio
    from rasterio.windows import from_bounds
    from pyproj import Transformer

    p = item["properties"]
    epsg = p.get("proj:epsg") or int(str(p.get("proj:code")).split(":")[-1])
    tr = Transformer.from_crs(4326, epsg, always_xy=True)
    x, y = tr.transform(lon, lat)
    # grid azimuth of true north at the tile centre, for rotating to true north
    xn, yn = tr.transform(lon, lat + 0.01)
    conv = math.degrees(math.atan2(xn - x, yn - y))
    with rasterio.open(item["assets"]["red"]["href"]) as src:
        win = from_bounds(x - HALF_M, y - HALF_M, x + HALF_M, y + HALF_M, src.transform)
        z = src.read(1, window=win, boundless=False).astype("float32")
    return z, conv


def polar(P, kx, ky, conv):
    k = np.hypot(kx, ky)
    with np.errstate(divide="ignore"):
        wl = np.where(k > 0, 1.0 / k, np.inf)
    m = (wl >= WL_EDGES[0]) & (wl < WL_EDGES[-1])
    wi = np.digitize(wl[m], WL_EDGES) - 1
    # grid azimuth -> true azimuth: subtract the grid azimuth of true north
    az = (np.degrees(np.arctan2(kx[m], ky[m])) - conv) % 180.0
    ai = az.astype(int) % 180
    flat = wi * 180 + ai
    s = np.bincount(flat, weights=P[m], minlength=40 * 180).reshape(40, 180)
    n = np.bincount(flat, minlength=40 * 180).reshape(40, 180)
    mean = np.where(n > 0, s / np.maximum(n, 1), np.nan)
    return np.log10(mean + 1e-12).astype(np.float16), n.astype(np.uint16)


def process(row):
    half_lat = HALF_M / 111320.0
    half_lon = half_lat / math.cos(math.radians(row.lat))
    items, poly = search(row.lon, row.lat, half_lon, half_lat)
    scenes = pick_scenes(items, poly)
    out = {"tile": row.tile, "n_items": len(items)}
    for season, it in scenes.items():
        z, conv = read_window(it, row.lon, row.lat)
        if (z == 0).mean() > 0.01 or z.shape[0] < 900 or z.shape[1] < 900:
            out[f"{season}_status"] = "bad_window"
            continue
        P, kx, ky = power_spectrum(z, PIX, PIX)
        logp, n = polar(P, kx, ky, conv)
        np.savez_compressed(OUTDIR / f"{row.tile}_{season}.npz", logp=logp, n=n)
        p = it["properties"]
        out.update({f"{season}_status": "ok", f"{season}_date": p["datetime"][:10],
                    f"{season}_sun_az": p.get("view:sun_azimuth"), f"{season}_sun_el": p.get("view:sun_elevation"),
                    f"{season}_conv": conv, f"{season}_cloud": p.get("eo:cloud_cover")})
    for season in SEASONS:
        out.setdefault(f"{season}_status", "no_scene")
    return out


def main(workers=6):
    from addendum import base_tiles
    from analysis import load
    OUTDIR.mkdir(parents=True, exist_ok=True)
    base = base_tiles(load())
    log = PROC / "s2_tiles.csv"
    done = pd.read_csv(log) if log.exists() else pd.DataFrame(columns=["tile"])
    todo = base[~base.tile.isin(done.tile)]
    print(f"{len(done)} done, {len(todo)} to do", flush=True)
    rows, n = [], 0
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(process, r): r.tile for r in todo.itertuples()}
        for fu in as_completed(futs):
            try:
                rows.append(fu.result())
            except Exception as err:
                print("tile", futs[fu], "failed:", err, flush=True)
            n += 1
            if len(rows) >= 50 or n == len(futs):
                new = pd.DataFrame(rows)
                done = pd.concat([done, new], ignore_index=True) if len(done) else new
                done.to_csv(log, index=False)
                rows = []
                print(f"{n}/{len(futs)}", flush=True)


if __name__ == "__main__":
    sys.exit(main())
