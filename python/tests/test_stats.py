import numpy as np
import pytest

from gsua_csb import band_depth, median_ci


def test_median_ci_bounds_are_ordered_and_sane():
    rng = np.random.default_rng(0)
    x = rng.normal(loc=10, scale=2, size=(3, 500))
    lo, hi = median_ci(x, alpha=0.05)
    assert np.all(lo < hi)
    # 500 samples should give a tight interval around the true median (10)
    assert np.all(np.abs(lo - 10) < 1.0)
    assert np.all(np.abs(hi - 10) < 1.0)


def test_median_ci_empirical_coverage_near_nominal():
    # A single draw's containment is not a valid test of a 95% CI (~14% chance of a false fail
    # across 3 independent rows); check the empirical coverage rate over many trials instead.
    rng = np.random.default_rng(0)
    hits = 0
    trials = 500
    for _ in range(trials):
        x = rng.normal(loc=10, scale=2, size=(1, 100))
        lo, hi = median_ci(x, alpha=0.05)
        hits += int(lo[0] < 10 < hi[0])
    coverage = hits / trials
    assert 0.90 <= coverage <= 0.99


def test_median_ci_shrinks_with_more_samples():
    rng = np.random.default_rng(1)
    small = rng.normal(size=(1, 20))
    large = rng.normal(size=(1, 2000))
    lo_s, hi_s = median_ci(small)
    lo_l, hi_l = median_ci(large)
    assert (hi_l - lo_l) < (hi_s - lo_s)


def test_band_depth_ranks_central_curve_first():
    t = np.linspace(0, 1, 30)
    rng = np.random.default_rng(2)
    center = np.sin(2 * np.pi * t)
    curves = center[None, :] + 0.1 * rng.standard_normal((40, 30))
    outlier = center + 5.0
    x = np.vstack([curves, outlier[None, :]])  # outlier is the last row

    depth, idx = band_depth(x)
    assert idx[-1] == x.shape[0] - 1  # outlier should rank as the shallowest


def test_band_depth_multi_output_shape():
    t = np.linspace(0, 1, 20)
    rng = np.random.default_rng(3)
    x = rng.standard_normal((15, 20, 2))
    depth, idx = band_depth(x)
    assert depth.shape == (15,)
    assert idx.shape == (15,)
