"""Post hoc power check (not in the pre-registered plan).

If giant dunes had perfect LGM memory, observed crests would scatter around
the LGM prediction exactly as much as they actually scatter around the modern
one. We plant that: synthetic crest = LGM prediction + a residual drawn from the
real (observed minus modern) residual of the same tile, with one random sign per
dune field. Keeping residuals in place preserves how errors cluster within fields,
which is what makes the real field-bootstrap intervals wide. Repeated 200 times.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis import RES, field_boot_mean, h1_parts, load, pred
from flux import signed_axial

def main():
    d = load()
    base = d[(d.large_contrast >= 2) & (d.large_type != "s") & (d.mod_bi_sd <= 10) & (d.mod_f_sd.fillna(99) <= 10)].copy()
    resid = signed_axial(base.large_crest.values, pred(base, "mod", "type"))
    lgm = pred(base, "LGM_ens", "type")
    rng = np.random.default_rng(3)
    gains, detected = [], 0
    for i in range(200):
        codes, uniq = pd.factorize(base.field)
        r = resid * rng.choice([-1.0, 1.0], len(uniq))[codes]
        base["synthetic"] = (lgm + r) % 180.0
        m_mod, m_past = h1_parts(base, "synthetic", "LGM_ens", "type")
        g = np.nanmean(m_mod - m_past)
        lo, _ = field_boot_mean(m_mod - m_past, base.field.values, n=500, seed=i)
        gains.append(g)
        detected += lo > 0
    out = {"planted_gain_mean": float(np.mean(gains)), "planted_gain_5_95": np.percentile(gains, [5, 95]).tolist(),
           "power": detected / 200, "median_predicted_rotation": float(np.nanmedian(np.abs(signed_axial(lgm, pred(base, "mod", "type")))))}
    (RES / "power.json").write_text(json.dumps(out, indent=1))
    print(out)

    # per-model H1 (secondary robustness)
    rows = []
    for c in [c for c in d.columns if c.startswith("LGM_") and c.endswith("_bi") and "ens" not in c]:
        tag = c[:-3]
        m_mod, m_past = h1_parts(base, "large_crest", tag, "type")
        lo, hi = field_boot_mean(m_mod - m_past, base.field.values)
        rows.append(dict(model=tag[4:], gain=np.nanmean(m_mod - m_past), lo=lo, hi=hi,
                         rotation=np.nanmedian(np.abs(signed_axial(pred(base, tag, "type"), pred(base, "mod", "type"))))))
    t = pd.DataFrame(rows)
    t.to_csv(RES / "h1_per_model.csv", index=False)
    print(t.round(2).to_string())

if __name__ == "__main__":
    main()
