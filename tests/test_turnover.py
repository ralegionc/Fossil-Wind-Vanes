"""Known-answer check for the within-field turnover regression."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from turnover import fe_slope  # noqa: E402


def test_recovers_planted_slope_despite_field_offsets():
    rng = np.random.default_rng(0)
    fields = np.repeat(np.arange(20), 50)
    offset = rng.normal(0, 20, 20)[fields]          # big differences between fields
    x = rng.normal(3, 0.5, fields.size) + 0.1 * offset  # x correlated with field offsets
    z = rng.uniform(0, 1, fields.size)
    y = offset + 4.0 * x - 6.0 * z + rng.normal(0, 3, fields.size)
    assert abs(fe_slope(y, x, [z], fields) - 4.0) < 0.6


def test_flat_within_fields_gives_zero_even_if_pooled_trend():
    fields = np.repeat(np.arange(10), 30)
    x = np.tile(np.linspace(2, 4, 30), 10) + fields * 0.3
    y = 10.0 * fields + np.random.default_rng(1).normal(0, 1, x.size)  # trend only between fields
    assert abs(fe_slope(y, x, [], fields)) < 1.0
