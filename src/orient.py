"""Directional elevation spectra of each tile from the Copernicus 30 m DEM.

For each Gunn tile we read the DEM window straight from the public AWS bucket
(the DEM itself is never written to disk), take the 2D power spectrum of the
detrended, Hann-windowed elevation, and save it compactly on a polar grid of
wavelength by direction. Crest orientations are measured from these spectra in
spectra.py.

Angles are compass azimuths (clockwise from north) of the crest line, axial,
in [0, 180).
"""
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import math
import os
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "processed" / "dem_tiles.csv"

os.environ.setdefault("AWS_NO_SIGN_REQUEST", "YES")
os.environ.setdefault("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
os.environ.setdefault("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", ".tif")
os.environ.setdefault("GDAL_HTTP_MAX_RETRY", "4")
os.environ.setdefault("GDAL_HTTP_RETRY_DELAY", "2")

HALF = 0.140625          # half of Gunn's 0.28125 deg tile
ARCSEC_M = 30.87         # metres per arcsecond of latitude
MIN_WL = 120.0           # 4 pixels: anything shorter is not resolved


def cop_url(lat_floor, lon_floor):
    ns = "N" if lat_floor >= 0 else "S"
    ew = "E" if lon_floor >= 0 else "W"
    name = f"Copernicus_DSM_COG_10_{ns}{abs(lat_floor):02d}_00_{ew}{abs(lon_floor):03d}_00_DEM"
    return f"/vsicurl/https://copernicus-dem-30m.s3.amazonaws.com/{name}/{name}.tif"


def read_dem(lat, lon, half=HALF):
    """Elevation window centred on (lat, lon); rows run north to south."""
    import rasterio
    from rasterio.merge import merge

    w, s, e, n = lon - half, lat - half, lon + half, lat + half
    srcs = []
    for la in range(math.floor(s), math.floor(n - 1e-9) + 1):
        for lo in range(math.floor(w), math.floor(e - 1e-9) + 1):
            try:
                srcs.append(rasterio.open(cop_url(la, lo)))
            except rasterio.errors.RasterioIOError:
                pass  # ocean tile, not in the bucket
    if not srcs:
        return None
    try:
        res = (1 / 3600, 1 / 3600)
        arr, _ = merge(srcs, bounds=(w, s, e, n), res=res, nodata=np.nan, dtype="float32")
    finally:
        for src in srcs:
            src.close()
    z = arr[0]
    if np.isnan(z).mean() > 0.02:
        return None
    return np.where(np.isnan(z), np.nanmean(z), z)


def power_spectrum(z, dx, dy):
    """Return power, east wavenumber and north wavenumber (cycles per metre)."""
    ny, nx = z.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    A = np.c_[xx.ravel(), yy.ravel(), np.ones(z.size)]
    coef, *_ = np.linalg.lstsq(A, z.ravel(), rcond=None)
    z = z - (A @ coef).reshape(z.shape)
    win = np.outer(np.hanning(ny), np.hanning(nx))
    P = np.abs(np.fft.fft2(z * win)) ** 2
    kx = np.fft.fftfreq(nx, dx)[None, :] * np.ones((ny, 1))
    ky = -np.fft.fftfreq(ny, dy)[:, None] * np.ones((1, nx))  # rows run south
    return P, kx, ky


def angular_profile(P, kx, ky, wl_lo, wl_hi, smooth_deg=5.0):
    """Smoothed axial histogram of spectral power over wave-vector azimuth (1 deg bins)."""
    k = np.hypot(kx, ky)
    with np.errstate(divide="ignore"):
        wl = np.where(k > 0, 1.0 / k, np.inf)
    m = (wl >= wl_lo) & (wl <= wl_hi)
    if m.sum() < 20:
        return None, None, None
    p = P[m]
    wave_az = np.degrees(np.arctan2(kx[m], ky[m])) % 180.0
    hist, _ = np.histogram(wave_az, bins=180, range=(0, 180), weights=p)
    x = np.arange(-30, 31)
    ker = np.exp(-0.5 * (x / smooth_deg) ** 2)
    sm = np.convolve(np.r_[hist[-30:], hist, hist[:30]], ker, mode="same")[30:-30]
    c = np.sum(p * np.exp(2j * np.radians(wave_az)))
    R = np.abs(c) / p.sum()
    share = p.sum() / P[(wl >= MIN_WL) & np.isfinite(wl)].sum()
    return sm, R, share


def band_orientation(P, kx, ky, wl_lo, wl_hi, smooth_deg=5.0):
    """Crest azimuth (deg, axial), concentration R and band power share."""
    sm, R, share = angular_profile(P, kx, ky, wl_lo, wl_hi, smooth_deg)
    if sm is None:
        return np.nan, np.nan, np.nan
    wave_peak = (np.argmax(sm) + 0.5) % 180.0
    return (wave_peak + 90.0) % 180.0, R, share


def secondary_orientation(P, kx, ky, wl_lo, wl_hi, exclude_crest, wedge=15.0, smooth_deg=5.0):
    """Strongest crest set in the band once the wedge around exclude_crest is masked.

    Harmonics of the large dunes sit at the large dunes' own azimuth, so masking
    that wedge leaves only crest sets that are genuinely oblique to them.
    Contrast is peak height over the median of the profile (1 = isotropic).
    """
    sm, _, _ = angular_profile(P, kx, ky, wl_lo, wl_hi, smooth_deg)
    if sm is None or np.isnan(exclude_crest):
        return np.nan, np.nan
    centres = np.arange(180) + 0.5
    crest_of_bin = (centres + 90.0) % 180.0
    d = np.abs(crest_of_bin - exclude_crest) % 180.0
    d = np.minimum(d, 180.0 - d)
    masked = np.where(d > wedge, sm, -np.inf)
    j = int(np.argmax(masked))
    # require a real local maximum, not the shoulder of the masked wedge
    left, right = sm[(j - 1) % 180], sm[(j + 1) % 180]
    if not (sm[j] >= left and sm[j] >= right):
        return np.nan, float(sm[j] / np.median(sm))
    return float(crest_of_bin[j]), float(sm[j] / np.median(sm))


def peak_wavelength(P, kx, ky, wl_lo=MIN_WL, wl_hi=5000.0, nbins=60):
    """Wavelength of the peak of the radially averaged k-weighted spectrum."""
    k = np.hypot(kx, ky)
    m = (k > 1 / wl_hi) & (k < 1 / wl_lo)
    edges = np.logspace(np.log10(1 / wl_hi), np.log10(1 / wl_lo), nbins + 1)
    idx = np.digitize(k[m], edges) - 1
    s = np.bincount(idx, weights=P[m] * k[m] ** 2, minlength=nbins)[:nbins]
    n = np.bincount(idx, minlength=nbins)[:nbins]
    prof = np.where(n > 0, s / np.maximum(n, 1), 0)
    j = np.argmax(prof)
    return 1.0 / np.sqrt(edges[j] * edges[j + 1])


WL_EDGES = np.logspace(np.log10(MIN_WL), np.log10(6000.0), 41)   # 40 wavelength bins
POLAR = ROOT / "data" / "processed" / "polar"


def polar_spectrum(P, kx, ky):
    """Mean power on a (40 wavelength x 180 wave-azimuth) grid, as log10 float16.

    This compact form keeps everything the orientation analysis needs, so band
    choices can be revisited without downloading the DEM again.
    """
    k = np.hypot(kx, ky)
    with np.errstate(divide="ignore"):
        wl = np.where(k > 0, 1.0 / k, np.inf)
    m = (wl >= WL_EDGES[0]) & (wl < WL_EDGES[-1])
    wi = np.digitize(wl[m], WL_EDGES) - 1
    ai = (np.degrees(np.arctan2(kx[m], ky[m])) % 180.0).astype(int) % 180
    flat = wi * 180 + ai
    s = np.bincount(flat, weights=P[m], minlength=40 * 180).reshape(40, 180)
    n = np.bincount(flat, minlength=40 * 180).reshape(40, 180)
    mean = np.where(n > 0, s / np.maximum(n, 1), np.nan)
    return np.log10(mean + 1e-12).astype(np.float16), n.astype(np.uint16)


def process_tile(row):
    out = {"tile": row.tile}
    f = POLAR / f"{row.tile}.npz"
    z = read_dem(row.lat, row.lon)
    if z is None:
        out["status"] = "no_dem"
        return out
    dy = ARCSEC_M
    dx = ARCSEC_M * math.cos(math.radians(row.lat))
    P, kx, ky = power_spectrum(z, dx, dy)
    logp, n = polar_spectrum(P, kx, ky)
    np.savez_compressed(f, logp=logp, n=n)
    out["relief"] = float(np.percentile(z, 95) - np.percentile(z, 5))
    out["elev"] = float(np.median(z))
    out["peak_wl"] = float(peak_wavelength(P, kx, ky))
    out["status"] = "ok"
    return out


def main(limit=None, workers=6):
    POLAR.mkdir(parents=True, exist_ok=True)
    tiles = pd.read_csv(ROOT / "data" / "processed" / "tiles.csv")
    done = pd.read_csv(OUT) if OUT.exists() else pd.DataFrame(columns=["tile"])
    todo = tiles[~tiles.tile.isin(done.tile)]
    if limit:
        todo = todo.head(limit)
    print(f"{len(done)} done, {len(todo)} to do", flush=True)
    rows, n = [], 0
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(process_tile, r): r.tile for r in todo.itertuples()}
        for fu in as_completed(futs):
            try:
                rows.append(fu.result())
            except Exception as err:  # network hiccups: record and retry on next run
                print("tile", futs[fu], "failed:", err, flush=True)
            n += 1
            if len(rows) >= 50 or n == len(futs):
                new = pd.DataFrame(rows)
                done = pd.concat([done, new], ignore_index=True) if len(done) else new
                done.to_csv(OUT, index=False)
                rows = []
                print(f"{n}/{len(futs)}", flush=True)


if __name__ == "__main__":
    main(limit=int(sys.argv[1]) if len(sys.argv) > 1 else None)
