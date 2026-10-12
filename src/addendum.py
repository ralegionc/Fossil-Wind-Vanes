"""Addendum 1 tests from docs/PLAN.md: daily-wind predictions (A) and model-agreement subsets (B)."""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis import RES, field_boot_mean, h1_parts, load, pred
from flux import orientation_bi, orientation_f, signed_axial

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
DAILY_MODELS = {"AWI-ESM-1-1-LR": "AWIESM11LR", "MPI-ESM1-2-LR": "MPIESM12LR"}
MONTHLY = ["AWIESM11LR", "INMCM48", "MIROCES2L", "MPIESM12LR", "CNRMCM5",
           "GISSE2R", "IPSLCM5ALR", "MIROCESM", "MPIESMP", "MRICGCM3"]


def base_tiles(d):
    return d[(d.large_contrast >= 2) & (d.large_type != "s") & (d.mod_bi_sd <= 10)
             & (d.mod_f_sd.fillna(99) <= 10)].copy()


def expected_gains(s, past, mode="type", reps=100, seed=7):
    """Mean gain if crests scattered around the modern prediction (pure modern) or
    around the past prediction (perfect memory). Each tile keeps its real residual
    about the modern prediction, with one random sign per dune field: this keeps
    the size of the scatter and its clustering within fields, but removes any
    directional bias that would favour one hypothesis."""
    rng = np.random.default_rng(seed)
    resid = signed_axial(s.large_crest.values, pred(s, "mod", mode))
    codes, uniq = pd.factorize(s.field)
    out = {}
    for name, centre in (("if_modern", pred(s, "mod", mode)), ("if_memory", pred(s, past, mode))):
        g = []
        for _ in range(reps):
            tmp = s.copy()
            r = resid * rng.choice([-1.0, 1.0], len(uniq))[codes]
            tmp["large_crest"] = (centre + r) % 180.0
            a, b = h1_parts(tmp, "large_crest", past, mode)
            g.append(np.nanmean(a - b))
        out[f"expected_{name}"] = float(np.mean(g))
    return out


def h1_row(s, label, past, mode="type"):
    m_mod, m_past = h1_parts(s, "large_crest", past, mode)
    diff = m_mod - m_past
    lo, hi = field_boot_mean(diff, s.field.values) if len(s) >= 10 else (np.nan, np.nan)
    rot = np.abs(signed_axial(pred(s, past, mode), pred(s, "mod", mode)))
    row = dict(test=label, past=past, n=len(s), fields=s.field.nunique(),
               median_rotation=float(np.nanmedian(rot)) if len(s) else np.nan,
               misfit_modern=float(np.nanmean(m_mod)), misfit_past=float(np.nanmean(m_past)),
               gain=float(np.nanmean(diff)), lo=float(lo), hi=float(hi),
               interpreted=bool(s.field.nunique() >= 5))
    if len(s) >= 10:
        row.update(expected_gains(s, past, mode))
        # exploratory (not pre-registered): where the result sits between the two
        # benchmarks, 0 = follows today's winds, 1 = perfect LGM memory
        span = row["expected_if_memory"] - row["expected_if_modern"]
        for k, v in (("memory_fraction", row["gain"]), ("mf_lo", row["lo"]), ("mf_hi", row["hi"])):
            row[k] = (v - row["expected_if_modern"]) / span
    return row


def add_daily(d):
    """Columns LGMd{th}_{model}_{bi,f} and LGMd{th}_ens_{bi,f}: modern prediction turned by
    the crest rotation between each model's daily LGM and pre-industrial flux roses."""
    z = {}
    for f in PROC.glob("daily_roses_*.npz"):
        z.update(dict(np.load(f)))
    avail = {full: short for full, short in DAILY_MODELS.items() if f"{full}|tile_cell" in z}
    tiles = pd.read_csv(PROC / "tiles.csv")
    idx = d.tile.values
    for th in ("3", "4", "5"):
        rots = {"bi": [], "f": []}
        for full, short in avail.items():
            cell = z[f"{full}|tile_cell"][tiles.set_index("tile").index.get_indexer(idx)]
            lg, pi = z[f"{full}|lgm|rose{th}"], z[f"{full}|piControl|rose{th}"]
            for mode, fn in (("bi", orientation_bi), ("f", orientation_f)):
                r = signed_axial(fn(lg), fn(pi))[cell]
                rots[mode].append(r)
                d[f"LGMd{th}_{short}_{mode}"] = (d[f"mod_{mode}"].values + r) % 180.0
        for mode in ("bi", "f"):
            d[f"LGMd{th}_ens_{mode}"] = (d[f"mod_{mode}"].values + np.nanmean(rots[mode], 0)) % 180.0
    return d, list(avail.values())


def signed_rot(s, past):
    return signed_axial(pred(s, past, "type"), pred(s, "mod", "type"))


def main():
    d = load()
    have_daily = any(PROC.glob("daily_roses_*.npz"))
    daily = []
    if have_daily:
        d, daily = add_daily(d)
    base = base_tiles(d)
    rows = []

    # primary tests again, so every row of the summary carries the same columns
    rows.append(h1_row(base, "Primary: all tiles", "LGM_ens"))
    gap = np.abs(signed_rot(base, "LGM_ens"))
    rows.append(h1_row(base[gap > 10], "Primary: informative", "LGM_ens"))

    # B, monthly ensemble: >= 8 of 10 models agree in sign, median rotation >= 10 deg
    R = np.stack([signed_rot(base, f"LGM_{m}") for m in MONTHLY], 1)
    med = np.nanmedian(R, 1)
    agree = (np.sign(R) == np.sign(med)[:, None]).sum(1)
    b_month = base[(np.abs(med) >= 10) & (agree >= 8)]
    rows.append(h1_row(b_month, "B: 10 models agree", "LGM_ens"))

    if have_daily:
        tag = "" if len(daily) == 2 else f" ({', '.join(daily)} only)"
        for th in ("4", "3", "5"):
            past = f"LGMd{th}_ens"
            rows.append(h1_row(base, f"A: daily winds, {th} m/s{tag}", past))
            gap = np.abs(signed_rot(base, past))
            rows.append(h1_row(base[gap > 10], f"A: daily winds, {th} m/s, informative{tag}", past))
            if len(daily) == 2:
                for short in daily:
                    rows.append(h1_row(base, f"A: daily winds, {th} m/s, {short}", f"LGMd{th}_{short}"))
        if len(daily) == 2:
            ra = np.stack([signed_rot(base, f"LGMd4_{s}") for s in daily], 1)
            both = (np.sign(ra[:, 0]) == np.sign(ra[:, 1])) & (np.abs(ra).min(1) >= 10)
            rows.append(h1_row(base[both], "B: both daily models agree", "LGMd4_ens"))

    t = pd.DataFrame(rows)
    t.to_csv(RES / "addendum.csv", index=False)
    print(t.round(2).to_string())
    return d, base


if __name__ == "__main__":
    main()
