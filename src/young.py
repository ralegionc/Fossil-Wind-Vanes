"""Test D from docs/PLAN.md: young dunes from Sentinel-2.

Check 1: Sentinel-2 must see the giant crest where the DEM does (median
difference under 10 deg), or imagery directions are not trusted.
Check 2: young oblique crests must follow today's winds (mean misfit to the
bed-instability prediction under 35 deg), or the young-dune arm is uninformative.
Then: (i) within tiles, is the giant crest further from today's winds than the
young crest? (ii) young dunes added to the turnover-time test as the fastest scale.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from addendum import base_tiles
from analysis import RES, boot_weights, load, pred
from flux import axial_diff
import s2
from spectra import main_set, oblique_sets
from turnover import SEC_PER_YEAR, boot, fe_slope

PROC = Path(__file__).resolve().parents[1] / "data" / "processed"
WL_MID = np.sqrt(s2.WL_EDGES[:-1] * s2.WL_EDGES[1:])
YOUNG = (40.0, 150.0)


def profile(tile, season, lo, hi, smooth=5.0):
    f = PROC / "s2polar" / f"{tile}_{season}.npz"
    if not f.exists():
        return None
    d = np.load(f)
    P, n = 10.0 ** d["logp"].astype(float), d["n"].astype(float)
    rows = (WL_MID >= lo) & (WL_MID <= hi)
    if (n[rows] > 0).sum() < 20:
        return None
    w = np.where(np.isnan(P[rows]), 0.0, P[rows] * n[rows])
    p = w.sum(0) / np.maximum(n[rows].sum(0), 1)
    x = np.arange(-30, 31)
    k = np.exp(-0.5 * (x / smooth) ** 2)
    k /= k.sum()
    p = np.convolve(np.r_[p[-30:], p, p[:30]], k, mode="same")[30:-30]
    return p / np.median(p)          # scenes differ in brightness; compare shapes


def both(tile, lo, hi):
    ps = [profile(tile, s, lo, hi) for s in ("jj", "dj")]
    if any(p is None for p in ps):
        return None
    return (ps[0] + ps[1]) / 2.0


def field_mean_ci(x, fields, n=2000):
    x = np.asarray(x, float)
    codes, W = boot_weights(np.asarray(fields), n)
    S = np.bincount(codes, weights=x, minlength=W.shape[1])
    N = np.bincount(codes, minlength=W.shape[1]).astype(float)
    return float(x.mean()), *np.percentile((W @ S) / (W @ N), [2.5, 97.5])


def main():
    b = base_tiles(load()).reset_index(drop=True)
    rows = []
    for r in b.itertuples():
        g = both(r.tile, r.large_wl / 1.6, min(1.6 * r.large_wl, 3000.0))
        y = both(r.tile, *YOUNG)
        if g is None or y is None:
            continue
        gc, gcon, _ = main_set(g)
        yc, ycon, _ = main_set(y)
        obl = oblique_sets(y, r.large_crest)
        rows.append(dict(tile=r.tile, s2_giant=gc, s2_giant_contrast=gcon, young_main=yc,
                         young_main_contrast=ycon,
                         young_obl=obl[0][0] if obl else np.nan, young_obl_prom=obl[0][1] if obl else np.nan))
    y = b.merge(pd.DataFrame(rows), on="tile")
    y.to_csv(PROC / "young.csv", index=False)
    report = {"tiles_with_both_scenes": len(y), "fields": int(y.field.nunique())}

    # check 1: imagery sees the giant crest where the DEM does
    diff = axial_diff(y.s2_giant, y.large_crest)
    report["check1_median_abs_diff"] = float(np.median(diff))
    report["check1_passes"] = bool(np.median(diff) < 10)

    # check 2: oblique young crests follow today's winds
    o = y[y.young_obl.notna()].copy()
    report["share_with_oblique_young"] = float(len(o) / len(y))
    m_bi = axial_diff(o.young_obl, o.mod_bi)
    report["check2_young_misfit_bi"] = float(m_bi.mean())
    report["check2_young_misfit_best"] = float(np.fmin(m_bi, axial_diff(o.young_obl, o.mod_f)).mean())
    report["check2_passes"] = bool(m_bi.mean() < 35)
    # secondary, contaminated by giant-crest edges
    report["secondary_unmasked_young_misfit_bi"] = float(axial_diff(y.young_main, y.mod_bi).mean())

    if report["check1_passes"] and report["check2_passes"]:
        # (i) paired: giant minus young misfit to today's winds, best-of-two for both
        g_m = np.fmin(axial_diff(o.large_crest, o.mod_bi), axial_diff(o.large_crest, o.mod_f))
        y_m = np.fmin(axial_diff(o.young_obl, o.mod_bi), axial_diff(o.young_obl, o.mod_f))
        est, lo, hi = field_mean_ci(g_m - y_m, o.field.values)
        report["paired_giant_minus_young"] = {"n": len(o), "fields": int(o.field.nunique()),
                                              "mean": est, "ci": [lo, hi]}
        # (ii) young dunes as the fastest scale in the turnover test
        aspect = float(np.nanmedian(b.small_height / b.small_wl))
        dp = o.mod_dp.values * SEC_PER_YEAR
        lam = np.sqrt(YOUNG[0] * YOUNG[1])
        T_young = aspect * lam * lam / (2 * dp)
        T_giant = o.large_height.values * o.large_wl.values / (2 * dp)
        logT = np.log10(np.r_[T_giant, T_young])
        mis = np.r_[g_m, y_m]
        fields = np.r_[o.field.values, o.field.values]
        z = [np.tile((o.mod_rdp / o.mod_dp).values, 2)]
        est = fe_slope(mis, logT, z, fields)
        lo, hi = boot(mis, logT, z, fields)
        report["ladder_slope_deg_per_decade"] = {"est": float(est), "ci": [float(lo), float(hi)],
                                                 "young_T_median": float(np.median(T_young)),
                                                 "giant_T_median": float(np.median(T_giant))}
    (RES / "young.json").write_text(json.dumps(report, indent=1))
    print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
