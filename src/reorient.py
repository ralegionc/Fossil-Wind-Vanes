"""Test E from docs/PLAN.md: reorientation time for long dunes.

Turning a dune of length L by an angle a moves about H L^2 a / 4 of sand, fed at
DP x L, so turning it by about a radian takes T_rot = H L / (4 DP). For long
linear dunes this can be far longer than rebuilding the cross-section.
"""
import json

import numpy as np
import pandas as pd

from addendum import base_tiles
from analysis import RES, boot_weights, load, pred
from flux import axial_diff
from turnover import SEC_PER_YEAR, run

LGM_YEARS = 21000.0


def t_rot(b):
    dp = b.mod_dp.values * SEC_PER_YEAR
    return b.large_height.values * b.large_width.values / (4.0 * dp)


def main():
    b = base_tiles(load()).reset_index(drop=True)
    T = t_rot(b)
    logT = np.log10(T)
    misfit = axial_diff(b.large_crest.values, pred(b, "mod", "type"))

    q = pd.qcut(logT, 5, labels=False)
    by_q = pd.DataFrame({"q": q, "T": T, "misfit": misfit, "linear": b.large_type.eq("l")}).groupby("q").agg(
        T_median=("T", "median"), misfit_mean=("misfit", "mean"), linear_share=("linear", "mean"),
        n=("misfit", "size"))
    report = {"fastest_fifth_misfit": float(by_q.misfit_mean.iloc[0]),
              "check_passes": bool(by_q.misfit_mean.iloc[0] < 35),
              "T_rot_percentiles": dict(zip(["5", "25", "50", "75", "95"],
                                            np.percentile(T, [5, 25, 50, 75, 95]).round(0).tolist())),
              "share_over_21kyr": float(np.mean(T > LGM_YEARS)),
              "share_over_21kyr_linear": float(np.mean(T[b.large_type.eq("l")] > LGM_YEARS))}

    rows = [run(b, "Primary: T_rot, type mode, directionality term", logT, misfit)]
    best = np.fmin(axial_diff(b.large_crest, b.mod_bi), axial_diff(b.large_crest, b.mod_f))
    rows.append(run(b, "Best-of-two modes", logT, best))
    rows.append(run(b, "No directionality term", logT, misfit, use_dir=False))
    lin = b.large_type.eq("l").values
    rows.append(run(b[lin], "Linear dunes only", logT[lin], misfit[lin]))
    t = pd.DataFrame(rows)

    # slow (> 21 kyr) versus fast tiles, within fields that have both
    slow = T > LGM_YEARS
    has_both = b.groupby("field")["tile"].transform(lambda s: 0) == 0
    f_slow = pd.Series(slow).groupby(b.field.values).transform("any").values
    f_fast = pd.Series(~slow).groupby(b.field.values).transform("any").values
    keep = f_slow & f_fast & has_both.values
    if keep.sum() > 0:
        s = b[keep]
        m = misfit[keep]
        sl = slow[keep]
        # per-field difference (slow minus fast), averaged over fields, field bootstrap
        diffs = (pd.DataFrame({"f": s.field.values, "m": m, "slow": sl})
                 .groupby(["f", "slow"]).m.mean().unstack())
        d = (diffs[True] - diffs[False]).values
        rng = np.random.default_rng(0)
        bs = [rng.choice(d, len(d)).mean() for _ in range(2000)]
        report["slow_minus_fast"] = {"fields": int(len(d)), "tiles_slow": int(sl.sum()),
                                     "tiles_fast": int((~sl).sum()), "mean": float(d.mean()),
                                     "ci": np.percentile(bs, [2.5, 97.5]).tolist()}
    t.to_csv(RES / "reorient.csv", index=False)
    by_q.to_csv(RES / "reorient_quintiles.csv")
    (RES / "reorient.json").write_text(json.dumps(report, indent=1))
    print(by_q.round(2).to_string())
    print(json.dumps(report, indent=1))
    print(t.round(2).to_string())


if __name__ == "__main__":
    main()
