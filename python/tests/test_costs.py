import numpy as np
import pytest

from gsua_csb import coverage_metric, likecost, costf, costf_multi, rcostf
from gsua_csb._costs import _pearson_rows


def _pre_nan_tolerant_costf(ydata, yfunction, regulator, alpha=2.0):
    """Reimplementation of costf's pre-NaN-tolerant formula (plain sum/length, no masking), for
    byte-identity comparison against the current implementation on NaN-free data."""
    ydata = np.atleast_2d(np.asarray(ydata, dtype=np.float64))
    yfunction = np.atleast_2d(np.asarray(yfunction, dtype=np.float64))
    regulator = np.atleast_1d(np.asarray(regulator, dtype=np.float64))
    n_outputs, length = ydata.shape
    cost = np.sum((ydata - yfunction) ** 2, axis=1) / length / regulator
    r = _pearson_rows(ydata, yfunction)
    r = np.where(np.isnan(r), 1.0, r)
    cost = ((2 - r) * cost) ** alpha
    return float(cost.sum() / n_outputs)


def test_costf_perfect_fit_is_zero():
    y = np.sin(np.linspace(0, 2 * np.pi, 50))[None, :]
    assert costf(y, y, regulator=[1.0]) == pytest.approx(0.0, abs=1e-10)


def test_rcostf_perfect_fit_is_zero():
    y = np.sin(np.linspace(0, 2 * np.pi, 50))[None, :]
    assert rcostf(y, y) == pytest.approx(0.0, abs=1e-10)


def test_rcostf_worse_fit_costs_more():
    t = np.linspace(0, 2 * np.pi, 50)
    ydata = np.sin(t)[None, :]
    close = ydata + 0.01 * np.cos(t)[None, :]
    far = ydata + 2.0 * np.cos(t)[None, :]
    assert rcostf(ydata, close) < rcostf(ydata, far)


def test_costf_multi_ranks_runs_by_fit_quality():
    t = np.linspace(0, 1, 40)
    ref = (np.sin(2 * np.pi * t) * 5 + 2)[None, :]
    rng = np.random.default_rng(0)
    good_runs = ref[0] + 0.02 * rng.standard_normal((50, 40))
    bad_runs = ref[0] * 0.3 + 3 + 0.02 * rng.standard_normal((50, 40))
    ydata = np.stack([good_runs, bad_runs], axis=0).reshape(100, 40, 1)
    cost = costf_multi(ydata, ref)
    assert cost[:50].mean() < cost[50:].mean()


def test_likecost_zero_when_matching_margin_baseline():
    t = np.linspace(0, 1, 30)
    ynom = np.sin(2 * np.pi * t) * 5 + 10
    margin = 0.1
    ymargin = ynom * (1 + abs(margin))
    cost = likecost(ymargin[None, :], ynom, margin)
    assert cost[0] == pytest.approx(0.0, abs=1e-8)


class TestCoverageMetric:
    """Mirrors the MATLAB gsua_covmetric validation scenarios."""

    def setup_method(self):
        rng = np.random.default_rng(1)
        self.rng = rng
        t = np.linspace(0, 1, 50)
        self.ydata = np.sin(2 * np.pi * t) * 5 + 2

    def test_good_fit_thin_band_both_below_one(self):
        y = self.ydata[None, :] + 0.02 * self.rng.standard_normal((200, 50))
        cost_data, cost_band, *_ = coverage_metric(y, self.ydata)
        assert cost_data < 1
        assert cost_band < 1

    def test_poor_fit_thin_band_distinguishes_accuracy_from_precision(self):
        y = (self.ydata * 0.3 + 3)[None, :] + 0.02 * self.rng.standard_normal((200, 50))
        cost_data, cost_band, *_ = coverage_metric(y, self.ydata)
        assert cost_data >= 1  # poor fit
        assert cost_band < 1  # but converged/thin

    def test_good_median_wide_band(self):
        y = self.ydata[None, :] + 3 * self.rng.standard_normal((200, 50))
        cost_data, cost_band, *_ = coverage_metric(y, self.ydata)
        assert cost_band >= 1  # not converged

    def test_zero_variance_band_no_nan(self):
        y = np.tile(self.ydata, (200, 1))
        cost_data, cost_band, p5, p50, p95 = coverage_metric(y, self.ydata)
        assert not np.isnan(cost_data)
        assert not np.isnan(cost_band)
        assert cost_band == pytest.approx(0.0, abs=1e-10)

    def test_excludes_nonfinite_runs_with_warning(self):
        y = self.ydata[None, :] + 0.02 * self.rng.standard_normal((200, 50))
        y[:5] = np.inf
        with pytest.warns(UserWarning, match="non-finite"):
            cost_data, cost_band, *_ = coverage_metric(y, self.ydata)
        assert np.isfinite(cost_data)
        assert np.isfinite(cost_band)

    def test_multi_output(self):
        ydata2 = np.stack([self.ydata, self.ydata * 0.5 + 1], axis=0)
        y = np.stack(
            [ydata2[0][None, :] + 0.02 * self.rng.standard_normal((150, 50)),
             ydata2[1][None, :] + 0.02 * self.rng.standard_normal((150, 50))],
            axis=2,
        ).reshape(150, 50, 2)
        cost_data, cost_band, p5, p50, p95 = coverage_metric(y, ydata2)
        assert p50.shape == (2, 50)
        assert cost_data < 1


def test_pearson_rows_matches_corrcoef_when_no_nan():
    rng = np.random.default_rng(2)
    a = rng.standard_normal((3, 20))
    b = a * 2 + 1 + 0.1 * rng.standard_normal((3, 20))
    r = _pearson_rows(a, b)
    expected = np.array([np.corrcoef(a[i], b[i])[0, 1] for i in range(3)])
    np.testing.assert_allclose(r, expected, atol=1e-10)


def test_pearson_rows_below_two_valid_points_is_zero():
    a = np.array([[1.0, np.nan, np.nan], [1.0, 2.0, 3.0]])
    b = np.array([[2.0, 3.0, 4.0], [2.0, 4.0, 6.0]])
    r = _pearson_rows(a, b)
    assert r[0] == 0.0  # row 0 has only 1 jointly-valid entry
    assert r[1] == pytest.approx(1.0)  # row 1 fully valid, perfectly correlated


def test_costf_no_nan_propagates_to_cost():
    ydata = np.array([[1.0, 2.0, np.nan, np.nan, np.nan, 6.0], [2.0, 4.0, 6.0, 8.0, 10.0, 12.0]])
    yfunction = ydata + 0.5
    cost = costf(ydata, yfunction, regulator=[1.0, 1.0])
    assert np.isfinite(cost)


def test_costf_nan_row_matches_manual_restriction():
    ydata = np.array([[1.0, 2.0, np.nan, np.nan, np.nan, 6.0]])
    yfunction = ydata + 0.5
    mask = ~np.isnan(ydata[0])
    restricted = costf(ydata[:, mask], yfunction[:, mask], regulator=[1.0])
    full = costf(ydata, yfunction, regulator=[1.0])
    assert full == pytest.approx(restricted, abs=1e-10)


def test_costf_byte_identical_to_pre_nan_tolerant_formula_when_no_nan():
    rng = np.random.default_rng(3)
    ydata = rng.uniform(1, 10, (2, 15))
    yfunction = ydata + rng.normal(0, 0.5, (2, 15))
    regulator = np.array([2.0, 3.0])
    assert costf(ydata, yfunction, regulator) == pytest.approx(
        _pre_nan_tolerant_costf(ydata, yfunction, regulator), abs=1e-12
    )


def test_rcostf_no_nan_propagates_to_cost():
    ydata = np.array([[5.0, 6.0, np.nan, np.nan, 9.0, 10.0], [2.0, 4.0, 6.0, 8.0, 10.0, 12.0]])
    yfunction = ydata + 0.5
    cost = rcostf(ydata, yfunction)
    assert np.isfinite(cost)
