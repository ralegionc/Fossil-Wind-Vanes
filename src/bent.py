"""Test F from docs/PLAN.md: young dunes under winds bent by the giant dunes.

Near the ground, giant dunes steer the wind toward their crest line: the turn is
largest when the wind meets the crest at 45 deg and zero along or across it.
Each ERA5 hour is turned by delta = -k sin(2 (theta - crest)). The strength k is
chosen on half the dune fields and scored on the other half (two-fold
cross-fitting), so the young-dune check is made only on held-out tiles.
"""
import json

import numpy as np
import pandas as pd

from analysis import PROC, RES, boot_weights
from flux import axial_diff, flux_rose, orientation_bi, orientation_f

KS = np.arange(0, 50, 5)


def bent_predictions(u, v, crest, k):
    theta = np.degrees(np.arctan2(u, v))                    # downwind azimuth
    delta = -k * np.sin(2 * np.radians(theta - crest[:, None]))
    speed = np.hypot(u, v)
    t2 = np.radians(theta + delta)
    rose = flux_rose(speed * np.sin(t2), speed * np.cos(t2))
    return orientation_bi(rose), orientation_f(rose)


def field_mean_ci(x, fields, n=2000):
    codes, W = boot_weights(np.asarray(fields), n)
    S = np.bincount(codes, weights=x, minlength=W.shape[1])
    N = np.bincount(codes, minlength=W.shape[1]).astype(float)
    return float(np.mean(x)), *np.percentile((W @ S) / (W @ N), [2.5, 97.5])


def main():
    y = pd.read_csv(PROC / "young.csv")
    o = y[y.young_obl.notna()].reset_index(drop=True)
    tiles = pd.read_csv(PROC / "tiles.csv")
    era = np.load(PROC / "era5_winds.npz")
    rows = tiles.set_index("tile").index.get_indexer(o.tile.values)
    u, v = era["u"][rows].astype(float), era["v"][rows].astype(float)
    crest = o.large_crest.values

    bi = {}
    f = {}
    for k in KS:
        bi[k], f[k] = bent_predictions(u, v, crest, k)
    mis = {k: axial_diff(o.young_obl.values, bi[k]) for k in KS}

    fields = o.field.unique()
    rng = np.random.default_rng(0)
    half = set(rng.permutation(fields)[: len(fields) // 2])
    fold = np.where(o.field.isin(half), 0, 1)
    chosen = {}
    held_bi = np.full(len(o), np.nan)
    held_best = np.full(len(o), np.nan)
    for test in (0, 1):
        train = fold != test
        k_star = int(min(KS, key=lambda k: mis[k][train].mean()))
        chosen[test] = k_star
        t = fold == test
        held_bi[t] = mis[k_star][t]
        held_best[t] = np.fmin(mis[k_star][t], axial_diff(o.young_obl.values[t], f[k_star][t]))

    # chance: young crests shuffled between held-out tiles
    chance = np.mean([axial_diff(rng.permutation(o.young_obl.values), np.where(fold == 0, bi[chosen[0]], bi[chosen[1]])).mean()
                      for _ in range(300)])
    report = {"n_tiles": len(o), "fields": int(len(fields)),
              "k_chosen_per_fold": {str(k): int(v) for k, v in chosen.items()},
              "misfit_by_k_all_tiles": {int(k): float(mis[k].mean()) for k in KS},
              "heldout_young_misfit_bi": float(np.nanmean(held_bi)),
              "heldout_chance_shuffled": float(chance),
              "unbent_young_misfit_bi": float(mis[0].mean())}
    report["check_passes"] = bool(report["heldout_young_misfit_bi"] < 35)
    report["k_at_grid_top"] = any(v == KS.max() for v in chosen.values())

    if report["check_passes"]:
        giant = np.fmin(axial_diff(o.large_crest, o.mod_bi), axial_diff(o.large_crest, o.mod_f))
        est, lo, hi = field_mean_ci(giant - held_best, o.field.values)
        report["paired_giant_minus_young"] = {"mean": est, "ci": [lo, hi]}
    (RES / "bent.json").write_text(json.dumps(report, indent=1))
    print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
