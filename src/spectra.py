"""Crest orientations from the saved polar spectra.

For a wavelength band we average spectral power over the band to get a profile
of power against wave direction (1 deg bins), smooth it circularly, and read:

* the main crest set: the profile peak, with its contrast (peak over median)
  and the doubled-angle concentration R;
* oblique crest sets: other peaks with real prominence, at least 15 deg from
  the main crest. Harmonics of the main dunes share the main crest's direction,
  so they cannot produce an oblique set.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import find_peaks

ROOT = Path(__file__).resolve().parents[1]
POLAR = ROOT / "data" / "processed" / "polar"
WL_EDGES = np.logspace(np.log10(120.0), np.log10(6000.0), 41)
WL_MID = np.sqrt(WL_EDGES[:-1] * WL_EDGES[1:])
MIN_WL = 120.0
SMOOTH = 5.0
WEDGE = 15.0
MIN_PROM = 0.5      # oblique peak prominence, in units of the profile median


def load(tile):
    d = np.load(POLAR / f"{tile}.npz")
    return 10.0 ** d["logp"].astype(np.float64), d["n"].astype(np.float64)


def profile(power, n, lo, hi, smooth=SMOOTH):
    """Smoothed power vs wave azimuth for wavelengths in [lo, hi]."""
    rows = (WL_MID >= lo) & (WL_MID <= hi)
    if rows.sum() == 0:
        return None
    w = np.where(np.isnan(power[rows]), 0.0, power[rows] * n[rows])
    if (n[rows] > 0).sum() < 20:
        return None
    prof = w.sum(0) / np.maximum(n[rows].sum(0), 1)
    x = np.arange(-30, 31)
    ker = np.exp(-0.5 * (x / smooth) ** 2)
    ker /= ker.sum()
    ext = np.r_[prof[-30:], prof, prof[:30]]
    return np.convolve(ext, ker, mode="same")[30:-30]


def wave_to_crest(wave_az):
    return (np.asarray(wave_az) + 0.5 + 90.0) % 180.0


def axial_dist(a, b):
    d = np.abs(np.asarray(a) - np.asarray(b)) % 180.0
    return np.minimum(d, 180.0 - d)


def main_set(prof):
    j = int(np.argmax(prof))
    th = np.radians(np.arange(180) + 0.5)
    R = np.abs(np.sum(prof * np.exp(2j * th))) / prof.sum()
    return float(wave_to_crest(j)), float(prof[j] / np.median(prof)), float(R)


def oblique_sets(prof, main_crest, wedge=WEDGE, min_prom=MIN_PROM):
    """Prominent peaks away from main_crest, strongest first: list of (crest, prominence/median)."""
    med = np.median(prof)
    ext = np.r_[prof, prof, prof]
    peaks, props = find_peaks(ext, prominence=min_prom * med)
    out = {}
    for p, prom in zip(peaks, props["prominences"]):
        if 180 <= p < 360:
            c = float(wave_to_crest(p - 180))
            if axial_dist(c, main_crest) > wedge:
                out[c] = max(out.get(c, 0.0), prom / med)
    return sorted(out.items(), key=lambda kv: -kv[1])


def bands(large_wl, small_wl, two_scale):
    if two_scale and small_wl >= MIN_WL and large_wl > small_wl:
        split = np.sqrt(large_wl * small_wl)
        return (split, min(1.8 * large_wl, 6000.0)), (max(small_wl / 1.8, MIN_WL), split)
    return (max(large_wl / 1.6, MIN_WL), min(1.6 * large_wl, 6000.0)), None


def measure(row):
    power, n = load(row.tile)
    bl, bs = bands(row.large_wl, row.small_wl, bool(row.two_scale))
    out = {"tile": row.tile}
    p = profile(power, n, *bl)
    if p is None:
        return out
    out["large_crest"], out["large_contrast"], out["large_R"] = main_set(p)
    # short band: everything clearly shorter than the large dunes
    hi = 0.6 * (row.small_wl if (bs is not None) else row.large_wl)
    hi = max(hi, 1.5 * MIN_WL)
    q = profile(power, n, MIN_WL, hi)
    if q is not None:
        out["short_main_crest"], out["short_main_contrast"], _ = main_set(q)
        obl = oblique_sets(q, out["large_crest"])
        if obl:
            out["short_obl_crest"], out["short_obl_prom"] = obl[0]
            out["short_n_obl"] = len(obl)
    out["short_hi"] = hi
    if bs is not None:
        s = profile(power, n, *bs)
        if s is not None:
            out["small_crest"], out["small_contrast"], out["small_R"] = main_set(s)
            obl = oblique_sets(s, out["large_crest"])
            if obl:
                out["small_obl_crest"], out["small_obl_prom"] = obl[0]
    return out


def run():
    tiles = pd.read_csv(ROOT / "data" / "processed" / "tiles.csv")
    dem = pd.read_csv(ROOT / "data" / "processed" / "dem_tiles.csv")
    ok = tiles.merge(dem[dem.status == "ok"][["tile", "relief", "elev", "peak_wl"]], on="tile")
    res = pd.DataFrame([measure(r) for r in ok.itertuples()])
    out = ok.merge(res, on="tile", how="left")
    out.to_csv(ROOT / "data" / "processed" / "orientations.csv", index=False)
    return out


if __name__ == "__main__":
    o = run()
    print(len(o), "tiles measured")
    print(o[["large_contrast", "large_R"]].describe().round(2).to_string())
    print("oblique short sets:", o.short_obl_crest.notna().sum(), " oblique small sets:", o.small_obl_crest.notna().sum())
