import numpy as np
import pytest

from gsua_csb import coverage_metric, likecost, costf, costf_multi, rcostf


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
