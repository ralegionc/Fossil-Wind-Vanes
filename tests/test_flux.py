"""Known-answer checks for flux roses and the two orientation modes."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from flux import (axial_diff, flux_rose, orientation_bi, orientation_f,  # noqa: E402
                  resultant, sand_flux, signed_axial)


def winds(towards_deg, speed, n):
    a = np.radians(towards_deg)
    return np.full(n, speed * np.sin(a)), np.full(n, speed * np.cos(a))


def mix(*parts):
    u = np.concatenate([p[0] for p in parts])
    v = np.concatenate([p[1] for p in parts])
    return u, v


def test_threshold_and_direction():
    q, az = sand_flux(np.array([0.0, 3.0, 10.0]), np.array([10.0, 0.0, 0.0]))
    assert q[1] == 0.0 and q[0] > 0 and q[2] > 0      # 3 m/s is below threshold
    assert abs(az[0] - 0.0) < 1e-6 and abs(az[2] - 90.0) < 1e-6


def test_unidirectional_wind():
    rose = flux_rose(*winds(90.0, 10.0, 100))  # sand moves east
    assert axial_diff(orientation_bi(rose), 0.0) < 2.0    # transverse crests run N-S
    assert axial_diff(orientation_f(rose), 90.0) < 2.0    # elongation along the wind


def test_rubin_hunter_obtuse_gives_longitudinal_crests():
    # two equal winds 120 deg apart: crests parallel to the resultant (linear dunes)
    rose = flux_rose(*mix(winds(60.0, 10.0, 100), winds(180.0, 10.0, 100)))
    rdd, _, _ = resultant(rose)
    assert abs(rdd - 120.0) < 2.0
    assert axial_diff(orientation_bi(rose), 120.0) < 3.0


def test_rubin_hunter_acute_gives_transverse_crests():
    rose = flux_rose(*mix(winds(60.0, 10.0, 100), winds(120.0, 10.0, 100)))
    assert axial_diff(orientation_bi(rose), 0.0) < 3.0


def test_fingering_follows_resultant_for_bimodal_wind():
    # unequal bimodal wind: elongation sits between the two, nearer the stronger
    rose = flux_rose(*mix(winds(30.0, 12.0, 100), winds(150.0, 9.0, 100)))
    f = orientation_f(rose)
    assert 30.0 < f < 90.0


def test_signed_axial():
    assert signed_axial(10.0, 170.0) == 20.0
    assert signed_axial(170.0, 10.0) == -20.0
