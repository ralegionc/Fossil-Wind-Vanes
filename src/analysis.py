"""Tests H1 to H3 and the placebo from docs/PLAN.md."""
from pathlib import Path
import json

import numpy as np
import pandas as pd

from flux import axial_diff, signed_axial

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
RES = ROOT / "results"
N_BOOT = 2000
N_PERM = 1000
INFORMATIVE = 10.0


def load():
    o = pd.read_csv(PROC / "orientations.csv")
    p = pd.read_csv(PROC / "predictions.csv")
    return o.merge(p, on="tile")


def pred(d, epoch, mode):
    """Predicted crest per tile for an epoch tag ('mod', 'LGM_ens', ...) and mode."""
    if mode == "type":
        f, b = d[f"{epoch}_f"], d[f"{epoch}_bi"]
        return np.where(d.large_type == "l", f, b)
    return d[f"{epoch}_{mode}"].values


def misfit(obs, d, epoch, mode):
    if mode == "best":
        return np.fmin(axial_diff(obs, d[f"{epoch}_bi"]), axial_diff(obs, d[f"{epoch}_f"]))
    return axial_diff(obs, pred(d, epoch, mode))


def boot_weights(fields, n=N_BOOT, seed=0):
    """Field codes per tile and an (n, n_fields) matrix of how often each field is drawn."""
    codes, uniq = pd.factorize(fields)
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, len(uniq), (n, len(uniq)))
    W = np.stack([np.bincount(r, minlength=len(uniq)) for r in draws]).astype(float)
    return codes, W


def field_boot_mean(x, fields, n=N_BOOT, seed=0):
    """95% CI of the tile mean of x when whole fields are resampled."""
    x = np.asarray(x, float)
    ok = ~np.isnan(x)
    codes, W = boot_weights(np.asarray(fields)[ok], n, seed)
    S = np.bincount(codes, weights=x[ok], minlength=W.shape[1])
    N = np.bincount(codes, minlength=W.shape[1]).astype(float)
    return np.percentile((W @ S) / (W @ N), [2.5, 97.5])


def field_boot_corr(x, y, fields, n=N_BOOT, seed=0):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = ~(np.isnan(x) | np.isnan(y))
    x, y = x[ok], y[ok]
    codes, W = boot_weights(np.asarray(fields)[ok], n, seed)
    m = W.shape[1]
    sums = {k: W @ np.bincount(codes, weights=v, minlength=m)
            for k, v in dict(x=x, y=y, xx=x * x, yy=y * y, xy=x * y, one=np.ones_like(x)).items()}
    n_ = sums["one"]
    cov = sums["xy"] / n_ - sums["x"] * sums["y"] / n_**2
    vx = sums["xx"] / n_ - (sums["x"] / n_) ** 2
    vy = sums["yy"] / n_ - (sums["y"] / n_) ** 2
    return np.nanpercentile(cov / np.sqrt(vx * vy), [2.5, 97.5])


def h1_parts(d, obs_col, past="LGM_ens", mode="type"):
    m_mod = misfit(d[obs_col].values, d, "mod", mode)
    m_past = misfit(d[obs_col].values, d, past, mode)
    return m_mod, m_past


def h1_stat(d, obs_col, past="LGM_ens", mode="type"):
    m_mod, m_past = h1_parts(d, obs_col, past, mode)
    return float(np.nanmean(m_mod - m_past)), float(np.nanmean(m_mod)), float(np.nanmean(m_past))


def informative(d, past="LGM_ens", mode="type"):
    if mode == "best":
        return d
    gap = axial_diff(pred(d, "mod", mode), pred(d, past, mode))
    return d[gap > INFORMATIVE]


def placebo(d, obs_col, mode="type", seed=0):
    """Each field's LGM crests become its modern predictions turned by the median
    predicted LGM rotation of another, randomly chosen field."""
    rng = np.random.default_rng(seed)
    fields = d.field.unique()
    real = h1_stat(d, obs_col, mode=mode)[0]
    null = []
    cols = ["LGM_ens_bi", "LGM_ens_f"]
    # rotation each field's LGM change applies to its own modern prediction
    rot = {c: signed_axial(d[c].values, d[c.replace("LGM_ens", "mod")].values) for c in cols}
    field_rot = {f: {c: np.nanmedian(rot[c][d.field.values == f]) for c in cols} for f in fields}
    dd = d.copy()
    for _ in range(N_PERM):
        perm = dict(zip(fields, rng.permutation(fields)))
        for c in cols:
            r = np.array([field_rot[perm[f]][c] for f in d.field.values])
            dd[c] = (d[c.replace("LGM_ens", "mod")].values + r) % 180.0
        null.append(h1_stat(dd, obs_col, mode=mode)[0])
    null = np.array(null)
    return real, float(np.mean(null >= real)), null


def main():
    RES.mkdir(exist_ok=True)
    d = load()
    base = d[(d.large_contrast >= 2) & (d.large_type != "s") & (d.mod_bi_sd <= 10) & (d.mod_f_sd.fillna(99) <= 10)]
    report = {"n_tiles_measured": int(d.large_crest.notna().sum()), "n_tiles_base": len(base),
              "n_fields_base": int(base.field.nunique())}

    rows = []
    for scale, obs, sub in [("large", "large_crest", base),
                            ("short_oblique", "short_obl_crest", base[base.short_obl_crest.notna()])]:
        for mode in (["type", "bi", "f", "best"] if scale == "large" else ["bi", "f", "best"]):
            for past in ("LGM_ens", "MH_ens"):
                for subset in ("all", "informative"):
                    s = sub if subset == "all" else informative(sub, past, "bi" if scale != "large" and mode == "type" else mode)
                    if len(s) < 10:
                        continue
                    est, mm, mp = h1_stat(s, obs, past, mode)
                    m_mod, m_past = h1_parts(s, obs, past, mode)
                    lo, hi = field_boot_mean(m_mod - m_past, s.field.values)
                    rows.append(dict(scale=scale, mode=mode, past=past, subset=subset, n=len(s),
                                     fields=s.field.nunique(), misfit_modern=mm, misfit_past=mp,
                                     gain=est, lo=lo, hi=hi))
    tab = pd.DataFrame(rows)
    tab.to_csv(RES / "h1_h2.csv", index=False)
    print(tab.round(2).to_string())

    # H3: does the gap between scales match the predicted rotation?
    s = base[base.short_obl_crest.notna()].copy()
    s["gap_obs"] = signed_axial(s.short_obl_crest, s.large_crest)
    s["rot_pred"] = signed_axial(pred(s, "mod", "bi"), pred(s, "LGM_ens", "bi"))
    r = float(np.corrcoef(s.gap_obs, s.rot_pred)[0, 1])
    lo, hi = field_boot_corr(s.gap_obs, s.rot_pred, s.field.values)
    report["h3"] = {"n": len(s), "r": r, "ci": [lo, hi]}
    print("H3 r =", round(r, 3), [round(lo, 3), round(hi, 3)], "n =", len(s))

    real, p, null = placebo(base, "large_crest")
    report["placebo_large_type"] = {"real": real, "p": p, "null_mean": float(null.mean()),
                                    "null_95": float(np.percentile(null, 95))}
    print("placebo: real", round(real, 2), "p =", p, "null mean", round(null.mean(), 2))
    (RES / "summary.json").write_text(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
