"""Known-answer checks for the spectral crest-orientation method."""
import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from orient import band_orientation, peak_wavelength, power_spectrum, secondary_orientation  # noqa: E402


def synthetic_field(crest_az, wl, dx=26.0, dy=30.87, n=1000, noise=0.5, seed=0):
    """Sinusoidal dunes whose crests run along compass azimuth crest_az."""
    rng = np.random.default_rng(seed)
    y = -np.arange(n)[:, None] * dy  # rows run south, so northing decreases
    x = np.arange(n)[None, :] * dx
    wave_az = math.radians(crest_az + 90.0)  # wave vector is perpendicular to crests
    phase = 2 * math.pi * (x * math.sin(wave_az) + y * math.cos(wave_az)) / wl
    return 10 * np.sin(phase) + noise * 10 * rng.standard_normal((n, n)) + 0.01 * x


def axial_diff(a, b):
    d = abs(a - b) % 180.0
    return min(d, 180.0 - d)


@pytest.mark.parametrize("crest", [0.0, 25.0, 60.0, 90.0, 135.0, 170.0])
def test_recovers_crest_azimuth(crest):
    z = synthetic_field(crest, 900.0)
    P, kx, ky = power_spectrum(z, 26.0, 30.87)
    got, R, share = band_orientation(P, kx, ky, 600.0, 1400.0)
    assert axial_diff(got, crest) < 2.0
    assert R > 0.8


def test_two_scales_are_separated():
    z = synthetic_field(20.0, 1500.0) + 0.6 * synthetic_field(110.0, 400.0, seed=1)
    P, kx, ky = power_spectrum(z, 26.0, 30.87)
    big, _, _ = band_orientation(P, kx, ky, 800.0, 2500.0)
    small, _, _ = band_orientation(P, kx, ky, 200.0, 800.0)
    assert axial_diff(big, 20.0) < 3.0
    assert axial_diff(small, 110.0) < 3.0


def test_secondary_ignores_harmonics_and_finds_oblique_set():
    # sawtooth-like large dunes carry strong harmonics in the small band
    big = synthetic_field(30.0, 1200.0, noise=0.1)
    harmonic = 0.5 * synthetic_field(30.0, 600.0, noise=0.0) + 0.3 * synthetic_field(30.0, 400.0, noise=0.0)
    oblique = 0.25 * synthetic_field(95.0, 450.0, noise=0.0)
    P, kx, ky = power_spectrum(big + harmonic + oblique, 26.0, 30.87)
    plain, _, _ = band_orientation(P, kx, ky, 300.0, 800.0)
    assert axial_diff(plain, 30.0) < 3.0  # the naive estimate is fooled by harmonics
    sec, contrast = secondary_orientation(P, kx, ky, 300.0, 800.0, exclude_crest=30.0)
    assert axial_diff(sec, 95.0) < 3.0
    assert contrast > 3.0


def test_secondary_returns_nothing_without_oblique_set():
    z = synthetic_field(30.0, 1200.0, noise=0.3) + 0.5 * synthetic_field(30.0, 600.0, noise=0.0)
    P, kx, ky = power_spectrum(z, 26.0, 30.87)
    sec, contrast = secondary_orientation(P, kx, ky, 300.0, 800.0, exclude_crest=30.0)
    assert np.isnan(sec) or contrast < 2.0


def test_peak_wavelength():
    z = synthetic_field(45.0, 1200.0, noise=0.2)
    P, kx, ky = power_spectrum(z, 26.0, 30.87)
    assert abs(peak_wavelength(P, kx, ky) - 1200.0) / 1200.0 < 0.1
