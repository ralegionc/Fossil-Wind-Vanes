"""Test G from docs/PLAN.md: weight the climate models by how well they reproduce
today's winds at the dune tiles."""
import json

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import xarray as xr

from addendum import MONTHLY, base_tiles, h1_row, signed_rot
from analysis import PROC, RES, load, pred
from flux import signed_axial
from predict import MODELS, at_tiles, control_for, model_delta, predictions
from slope import slope

SHORT = dict(zip(MODELS, MONTHLY))


def era5_monthly(tiles_idx):
    era = np.load(PROC / "era5_winds.npz")
    month = pd.to_datetime(era["hours"]).month.values - 1
    u, v = era["u"][tiles_idx].astype(float), era["v"][tiles_idx].astype(float)
    U = np.stack([u[:, month == m].mean(1) for m in range(12)], 1)
    V = np.stack([v[:, month == m].mean(1) for m in range(12)], 1)
    return U, V


def model_pi(model, lat, lon):
    ctrl = control_for(model, "lgm")
    out = []
    for var in ("uas", "vas"):
        da = xr.open_dataset(PROC / "pmip" / f"{model}_{ctrl}_{var}.nc")[var]
        out.append(at_tiles(da, lat, lon))
    return out


def rms(du, dv):
    return float(np.sqrt(np.mean(du ** 2 + dv ** 2)))


def main():
    d = load()
    b = base_tiles(d).reset_index(drop=True)
    tiles = pd.read_csv(PROC / "tiles.csv")
    row_of = tiles.set_index("tile").index.get_indexer(b.tile.values)
    U, V = era5_monthly(row_of)
    lat, lon = b.lat.values, b.lon.values
    rng = np.random.default_rng(0)

    skill = []
    for m in MODELS:
        mu, mv = model_pi(m, lat, lon)
        D = rms(mu - U, mv - V)
        null = []
        for _ in range(200):
            p = rng.permutation(len(b))
            null.append(rms(mu[p] - U, mv[p] - V))
        skill.append(dict(model=m, D=D, null_p5=float(np.percentile(null, 5)),
                          beats_null=bool(D < np.percentile(null, 5))))
    sk = pd.DataFrame(skill)
    s = np.median(sk.D)
    sk["weight"] = np.where(sk.beats_null, np.exp(-(sk.D / s) ** 2), 0.0)
    sk["weight"] /= sk.weight.sum()

    # skill-weighted LGM change -> predictions
    era = np.load(PROC / "era5_winds.npz")
    month = pd.to_datetime(era["hours"]).month.values - 1
    u, v = era["u"][row_of].astype(float), era["v"][row_of].astype(float)
    du = np.zeros((len(b), 12))
    dv = np.zeros((len(b), 12))
    for m, w in zip(sk.model, sk.weight):
        if w == 0:
            continue
        a, c = model_delta(m, "lgm", lat, lon)
        du += w * a
        dv += w * c
    for k, val in predictions(u + du[:, month], v + dv[:, month], "LGMw").items():
        b[k] = val

    rows = [h1_row(b, "Skill-weighted: all tiles", "LGMw")]
    gap = np.abs(signed_rot(b, "LGMw"))
    rows.append(h1_row(b[gap > 10], "Skill-weighted: informative", "LGMw"))
    h1 = pd.DataFrame(rows)

    err = signed_axial(b.large_crest.values, pred(b, "mod", "type"))
    turn = signed_axial(pred(b, "LGMw", "type"), pred(b, "mod", "type"))
    keep = np.abs(err) <= 60
    est, lo, hi = slope(err[keep], turn[keep], b.field.values[keep])

    # does skill track how well each model's LGM prediction fits the dunes?
    ms = pd.read_csv(RES / "memory_slope.csv")
    ms = ms[(ms.max_error == 60) & ms.prediction.str.startswith("Monthly winds, ") &
            ~ms.prediction.str.contains("mean|IV")]
    ms["short"] = ms.prediction.str.replace("Monthly winds, ", "", regex=False)
    pm = pd.read_csv(RES / "h1_per_model.csv")
    sk["short"] = sk.model.map(SHORT)
    j = sk.merge(ms[["short", "memory_fraction"]], on="short").merge(pm[["model", "gain"]].rename(
        columns={"model": "short"}), on="short")
    rho_slope = spearmanr(j.D, j.memory_fraction)
    rho_gain = spearmanr(j.D, j.gain)

    report = {"skill": sk.drop(columns="short").to_dict("records"),
              "weighted_memory_fraction": {"est": est, "ci": [lo, hi], "n": int(keep.sum())},
              "weighted_median_turn": float(np.median(np.abs(turn))),
              "spearman_D_vs_memory_fraction": [float(rho_slope.statistic), float(rho_slope.pvalue)],
              "spearman_D_vs_h1_gain": [float(rho_gain.statistic), float(rho_gain.pvalue)]}
    h1.to_csv(RES / "skill_h1.csv", index=False)
    (RES / "skill.json").write_text(json.dumps(report, indent=1))
    print(sk.round(3).to_string())
    print(h1.round(2).to_string())
    print(json.dumps({k: v for k, v in report.items() if k != "skill"}, indent=1))


if __name__ == "__main__":
    main()
