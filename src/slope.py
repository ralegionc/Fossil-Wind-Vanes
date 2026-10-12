"""Exploratory (not pre-registered): memory fraction as a regression slope.

For each tile, signed error = observed giant crest minus the modern prediction,
and predicted turn = LGM prediction minus the modern prediction (both axial, in
(-90, 90]). If crests follow today's winds the slope is 0; with perfect LGM memory
it is 1. Unlike the misfit difference used in the plan, this does not lose
leverage where errors approach 90 deg. Fit through the origin with a field
bootstrap; tiles whose error exceeds 60 deg are dropped as unmatched (the crest
belongs to neither prediction).
"""
import numpy as np
import pandas as pd

from addendum import DAILY_MODELS, MONTHLY, RES, add_daily, base_tiles
from analysis import boot_weights, load, pred
from flux import signed_axial


def slope(err, turn, fields, n=2000, seed=0):
    ok = np.isfinite(err) & np.isfinite(turn)
    err, turn, fields = err[ok], turn[ok], np.asarray(fields)[ok]
    codes, W = boot_weights(fields, n, seed)
    m = W.shape[1]
    sxy = np.bincount(codes, weights=err * turn, minlength=m)
    sxx = np.bincount(codes, weights=turn * turn, minlength=m)
    b = (W @ sxy) / (W @ sxx)
    return float(sxy.sum() / sxx.sum()), *np.percentile(b, [2.5, 97.5])


def split_half_iv(err, turns, fields, n_split=200, n=1000, seed=0):
    """Memory fraction corrected for model noise. The models are split into two
    random halves; one half's mean turn instruments the other's, so noise that is
    independent between models cancels: beta = sum(err * B) / sum(A * B).
    Averaged over random splits, with a field bootstrap of the averaged ratio."""
    rng = np.random.default_rng(seed)
    k = turns.shape[1]
    codes, W = boot_weights(np.asarray(fields), n, seed)
    m = W.shape[1]
    num, den = np.zeros(n), np.zeros(n)
    est_num = est_den = 0.0
    for _ in range(n_split):
        idx = rng.permutation(k)
        A = turns[:, idx[: k // 2]].mean(1)
        B = turns[:, idx[k // 2:]].mean(1)
        sn = np.bincount(codes, weights=err * B, minlength=m)
        sd = np.bincount(codes, weights=A * B, minlength=m)
        est_num += sn.sum()
        est_den += sd.sum()
        num += W @ sn
        den += W @ sd
    return float(est_num / est_den), *np.percentile(num / den, [2.5, 97.5])


def main():
    d, daily = add_daily(load())
    b = base_tiles(d)
    err_all = signed_axial(b.large_crest.values, pred(b, "mod", "type"))
    rows = []
    preds = {"Monthly winds, 10-model mean": "LGM_ens", "Daily winds, 2-model mean": "LGMd4_ens",
             "Daily winds, AWI-ESM": "LGMd4_AWIESM11LR", "Daily winds, MPI-ESM": "LGMd4_MPIESM12LR"}
    preds.update({f"Monthly winds, {m}": f"LGM_{m}" for m in MONTHLY})
    for label, past in preds.items():
        turn = signed_axial(pred(b, past, "type"), pred(b, "mod", "type"))
        for cut in (60, 90):
            keep = np.abs(err_all) <= cut
            est, lo, hi = slope(err_all[keep], turn[keep], b.field.values[keep])
            rows.append(dict(prediction=label, max_error=cut, n=int(keep.sum()),
                             median_turn=float(np.median(np.abs(turn[keep]))),
                             memory_fraction=est, lo=lo, hi=hi))
    turns = np.stack([signed_axial(pred(b, f"LGM_{m}", "type"), pred(b, "mod", "type")) for m in MONTHLY], 1)
    for cut in (60, 90):
        keep = np.abs(err_all) <= cut
        est, lo, hi = split_half_iv(err_all[keep], turns[keep], b.field.values[keep])
        rows.append(dict(prediction="Monthly winds, split-half IV (noise-corrected)", max_error=cut,
                         n=int(keep.sum()), median_turn=float(np.median(np.abs(turns[keep].mean(1)))),
                         memory_fraction=est, lo=lo, hi=hi))
    t = pd.DataFrame(rows)
    t.to_csv(RES / "memory_slope.csv", index=False)
    print(t.round(2).to_string())


if __name__ == "__main__":
    main()
