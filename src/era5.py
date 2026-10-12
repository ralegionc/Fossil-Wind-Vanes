"""Sample ERA5 hourly 10 m winds at every dune tile from the public ARCO-ERA5 store.

The store is chunked one hour per chunk over the whole globe, so pulling a full
hourly series at a point means reading every chunk. Instead we draw a seeded
random sample of hours from 1991 to 2020. Sand flux roses are averages over
time, so a uniform random sample gives an unbiased estimate of them, and the
diurnal and seasonal cycles are covered.

Output: data/processed/era5_winds.npz with u, v arrays of shape (n_tiles, n_hours)
and the sampled hour stamps.
"""
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import sys
import time

import gcsfs
import numcodecs
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
STORE = "gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3"
OUT = ROOT / "data" / "processed" / "era5_winds.npz"
PARTIAL = ROOT / "data" / "processed" / "era5_partial"

N_HOURS = 6000
SEED = 20261012
EPOCH = np.datetime64("1900-01-01T00")
START, END = np.datetime64("1991-01-01T00"), np.datetime64("2021-01-01T00")


def sample_hours(n=N_HOURS, seed=SEED):
    lo = int((START - EPOCH) / np.timedelta64(1, "h"))
    hi = int((END - EPOCH) / np.timedelta64(1, "h"))
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(np.arange(lo, hi), size=n, replace=False))


def grid_index(lat, lon):
    """Nearest cell on the 0.25 deg grid (lat 90 to -90, lon 0 to 359.75)."""
    i = np.rint((90.0 - np.asarray(lat)) / 0.25).astype(int)
    j = np.rint((np.mod(lon, 360.0)) / 0.25).astype(int) % 1440
    return i, j


def main():
    tiles = pd.read_csv(ROOT / "data" / "processed" / "tiles.csv")
    gi, gj = grid_index(tiles.lat.values, tiles.lon.values)
    hours = sample_hours()
    PARTIAL.mkdir(parents=True, exist_ok=True)

    fs = gcsfs.GCSFileSystem(token="anon")
    codec = numcodecs.Blosc()

    # sanity check: the time coordinate must be hours since 1900-01-01 with index == hour
    tchunk = np.frombuffer(codec.decode(fs.cat(f"{STORE}/time/0")), dtype="<i8")
    assert tchunk[0] == 0 and tchunk[1] == 1, tchunk[:3]

    def fetch(h):
        f = PARTIAL / f"{h}.npy"
        if f.exists():
            return h
        out = []
        for var in ("10m_u_component_of_wind", "10m_v_component_of_wind"):
            raw = fs.cat(f"{STORE}/{var}/{h}.0.0")
            arr = np.frombuffer(codec.decode(raw), dtype="<f4").reshape(721, 1440)
            out.append(arr[gi, gj])
        np.save(f, np.stack(out).astype(np.float32))
        return h

    t0 = time.time()
    done = 0
    with ThreadPoolExecutor(max_workers=12) as ex:
        futs = [ex.submit(fetch, int(h)) for h in hours]
        for fu in as_completed(futs):
            fu.result()
            done += 1
            if done % 250 == 0:
                rate = done / (time.time() - t0)
                print(f"{done}/{len(hours)} hours, {rate:.1f}/s", flush=True)

    uv = np.stack([np.load(PARTIAL / f"{int(h)}.npy") for h in hours], axis=-1)
    stamps = EPOCH + hours.astype("timedelta64[h]")
    np.savez_compressed(OUT, u=uv[0], v=uv[1], hours=stamps.astype("datetime64[h]").astype(str),
                        tile=tiles.tile.values, grid_i=gi, grid_j=gj)
    print("saved", OUT, uv.shape, f"{time.time() - t0:.0f}s")


if __name__ == "__main__":
    sys.exit(main())
