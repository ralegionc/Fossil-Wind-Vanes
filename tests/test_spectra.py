"""The polar-spectrum route must give the same answers as the full 2D spectrum."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from orient import polar_spectrum, power_spectrum  # noqa: E402
from spectra import axial_dist, main_set, oblique_sets, profile  # noqa: E402
from test_orient import synthetic_field  # noqa: E402


def polar(z):
    P, kx, ky = power_spectrum(z, 26.0, 30.87)
    logp, n = polar_spectrum(P, kx, ky)
    return 10.0 ** logp.astype(float), n.astype(float)


def test_main_set_from_polar():
    for crest in (5.0, 47.0, 133.0):
        power, n = polar(synthetic_field(crest, 1000.0))
        c, contrast, R = main_set(profile(power, n, 700.0, 1500.0))
        assert axial_dist(c, crest) < 2.0 and contrast > 5 and R > 0.5


def test_oblique_set_found_and_harmonics_ignored():
    big = synthetic_field(30.0, 1200.0, noise=0.1)
    harm = 0.5 * synthetic_field(30.0, 600.0, noise=0.0) + 0.3 * synthetic_field(30.0, 400.0, noise=0.0)
    obl = 0.25 * synthetic_field(95.0, 300.0, noise=0.0)
    power, n = polar(big + harm + obl)
    sets = oblique_sets(profile(power, n, 150.0, 700.0), main_crest=30.0)
    assert sets and axial_dist(sets[0][0], 95.0) < 3.0


def test_no_oblique_set_when_absent():
    z = synthetic_field(30.0, 1200.0, noise=0.3) + 0.5 * synthetic_field(30.0, 600.0, noise=0.0)
    power, n = polar(z)
    assert oblique_sets(profile(power, n, 150.0, 700.0), main_crest=30.0) == []
