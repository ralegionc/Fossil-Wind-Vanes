"""Test C from docs/PLAN.md: does the misfit to today's winds grow with turnover time?

T = H x wavelength / (2 x DP) is roughly how many years the wind takes to move a
dune's worth of sand. If giant dunes keep old orientations, slow dunes should sit
further from the crest today's winds predict. No ice-age model is involved.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from addendum import base_tiles
from analysis import RES, load, pred
from flux import axial_diff

SEC_PER_YEAR = 3.156e7


def fe_slope(y, x, z, fields):
    """Slope of y on x with field fixed effects and covariates z (columns)."""
    def demean(a):
        s = pd.Series(a)
        return (s - s.groupby(fields).transform("mean")).values
    X = np.column_stack([demean(x)] + [demean(c) for c in z])
    beta, *_ = np.linalg.lstsq(X, demean(y), rcond=None)
    return beta[0]


def boot(y, x, z, fields, n=2000, seed=0):
    rng = np.random.default_rng(seed)
    uniq = np.unique(fields)
    idx = {f: np.flatnonzero(fields == f) for f in uniq}
    vals = []
    for _ in range(n):
        pick = rng.choice(uniq, len(uniq))
        rows = np.concatenate([idx[f] for f in pick])
        # relabel repeated fields so each draw keeps its own fixed effect
        lab = np.concatenate([np.full(len(idx[f]), i) for i, f in enumerate(pick)])
        vals.append(fe_slope(y[rows], x[rows], [c[rows] for c in z], lab))
    return np.percentile(vals, [2.5, 97.5])


def prepare(b, height, wavelength):
    dp = b.mod_dp.values * SEC_PER_YEAR
    T = height * wavelength / (2.0 * dp)
    return np.log10(T)


def run(b, label, logT, misfit, use_dir=True):
    ok = np.isfinite(logT) & np.isfinite(misfit)
    s = b[ok]
    y, x = misfit[ok], logT[ok]
    z = [(s.mod_rdp / s.mod_dp).values] if use_dir else []
    est = fe_slope(y, x, z, s.field.values)
    lo, hi = boot(y, x, z, s.field.values)
    return dict(test=label, n=int(ok.sum()), fields=int(s.field.nunique()),
                slope_deg_per_decade=float(est), lo=float(lo), hi=float(hi))


def main():
    b = base_tiles(load()).reset_index(drop=True)
    misfit = axial_diff(b.large_crest.values, pred(b, "mod", "type"))
    logT = prepare(b, b.large_height.values, b.large_wl.values)

    # known-answer check: the fastest fifth must fit today's winds below 35 deg
    q = pd.qcut(logT, 5, labels=False)
    by_q = pd.DataFrame({"q": q, "T": 10 ** logT, "misfit": misfit}).groupby("q").agg(
        T_median=("T", "median"), misfit_mean=("misfit", "mean"), n=("misfit", "size"))
    check = bool(by_q.misfit_mean.iloc[0] < 35)

    rows = [run(b, "Primary: Gunn H and wavelength, type mode, directionality term", logT, misfit)]
    own = prepare(b, b.relief.values / 2.0, b.peak_wl.values)
    rows.append(run(b, "Own DEM relief/2 and peak wavelength", own, misfit))
    best = np.fmin(axial_diff(b.large_crest, b.mod_bi), axial_diff(b.large_crest, b.mod_f))
    rows.append(run(b, "Best-of-two modes", logT, best))
    rows.append(run(b, "No directionality term", logT, misfit, use_dir=False))
    lin = (b.large_type == "l").values
    rows.append(run(b[lin], "Linear dunes only", logT[lin], misfit[lin]))
    t = pd.DataFrame(rows)
    t.to_csv(RES / "turnover.csv", index=False)
    by_q.to_csv(RES / "turnover_quintiles.csv")
    (RES / "turnover_check.json").write_text(json.dumps({"fastest_fifth_misfit": float(by_q.misfit_mean.iloc[0]),
                                                          "passes": check}, indent=1))
    print(by_q.round(1).to_string())
    print("known-answer check passes:", check)
    print(t.round(2).to_string())


if __name__ == "__main__":
    main()
