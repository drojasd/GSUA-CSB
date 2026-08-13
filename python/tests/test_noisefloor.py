import warnings

import numpy as np
import pytest

from gsua_csb import UserFunctionModel, costf, noise_floor, parameter_estimation

XDATA = np.linspace(0, 4, 15)
TRUE_PARS = np.array([4.0, 0.4, 0.8])


def _fixture_func(pars, d):
    """2-output toy model: row 0 incident-like, row 1 cumulative-like (linear in pars[0], so
    scaling pars[0] and ydata together reproduces a joint (ydata,yfunction) rescale)."""
    row0 = pars[0] * np.exp(pars[1] * d)
    row1 = pars[2] * np.cumsum(row0)
    return np.vstack([row0, row1])


def _model():
    return UserFunctionModel(
        func=_fixture_func, names=["a", "b", "c"],
        range=np.array([[0.1, 20.0], [0.01, 2.0], [0.01, 5.0]]), domain=XDATA,
    )


def _fit_pool(ydata, n=1, margin=0.1, solver="minimize", initial_points=None):
    model = _model()
    pe = parameter_estimation(
        model, XDATA, ydata, n=n, solver=solver, margin=margin,
        initial_points=TRUE_PARS[None, :] if initial_points is None else initial_points,
    )
    return model, pe


class TestNb2Parameterization:
    """Verifies noise_floor's documented NB2 reparameterization (R=k, P=k/(k+mu)) directly
    against numpy's own negative_binomial, not by reaching into noise_floor's private bootstrap
    generator."""

    @pytest.mark.parametrize(
        "mu,k", [(50, 8), (2, 8), (100, 200)], ids=["moderate", "small_mu", "large_k"]
    )
    def test_nb2_mean_and_variance_match_target(self, mu, k):
        rng = np.random.default_rng(0)
        sample = rng.negative_binomial(k, k / (k + mu), size=200000).astype(float)
        target_var = mu + mu**2 / k
        assert sample.mean() == pytest.approx(mu, rel=0.02)
        assert sample.var() == pytest.approx(target_var, rel=0.05)


def test_quasipoisson_fallback_on_noiseless_data():
    # Noiseless data -> residuals identically zero at theta* -> pooled Pearson dispersion phi<=1
    # by construction -- quasi-Poisson's k=mu/(phi-1) reparameterization is degenerate there, so
    # noise_floor must fall back to plain Poisson draws with a warning, not error.
    ydata = _fixture_func(TRUE_PARS, XDATA)
    model, pe = _fit_pool(ydata)

    with pytest.warns(UserWarning, match="dispersion phi"):
        out = noise_floor(model, XDATA, ydata, pe.x, pe.cost, margin=pe.margin, alpha=pe.alpha,
                           b=2000, seed=1)

    assert out.quasi_fallback is True
    # Different point in the shared RNG stream than the separately-run 'poisson' entry, so not
    # bit-identical -- but both are literally Poisson(mu) draws, so summary stats should agree.
    assert out.by_model["quasipoisson"].cost.mean() == pytest.approx(
        out.by_model["poisson"].cost.mean(), rel=0.15
    )


def test_cumulative_synthetic_nan_mask_matches_real_data():
    ydata = _fixture_func(TRUE_PARS, XDATA)
    ydata[1, 5:8] = np.nan   # gap in the cumulative row
    model, pe = _fit_pool(ydata)

    out = noise_floor(model, XDATA, ydata, pe.x, np.array([0.01]), margin=pe.margin,
                       alpha=pe.alpha, b=50, seed=2, cumulative=[False, True])

    synth_row1 = out.example_synthetic["quasipoisson"][1]
    assert np.all(np.isnan(synth_row1[5:8]))
    assert np.all(np.isfinite(synth_row1[np.r_[0:5, 8:15]]))
    assert np.all(np.isfinite(out.example_synthetic["quasipoisson"][0]))
    assert not np.any(np.isnan(out.cost))


def test_cumulative_row_synthetic_sample_is_nondecreasing():
    ydata = _fixture_func(TRUE_PARS, XDATA)
    model, pe = _fit_pool(ydata)

    out = noise_floor(model, XDATA, ydata, pe.x, np.array([0.01]), margin=pe.margin,
                       alpha=pe.alpha, b=50, seed=3, cumulative=[False, True])

    synth_row1 = out.example_synthetic["quasipoisson"][1]
    assert np.all(np.diff(synth_row1) >= 0)


def test_best_cost_above_threshold_warns():
    ydata = _fixture_func(TRUE_PARS, XDATA)
    model, pe = _fit_pool(ydata)

    with pytest.warns(UserWarning, match="exceeds the noise-floor threshold"):
        noise_floor(model, XDATA, ydata, pe.x, np.array([1e12]), margin=pe.margin,
                    alpha=pe.alpha, b=50, seed=4)


def test_best_cost_below_threshold_does_not_warn_about_rejection():
    ydata = _fixture_func(TRUE_PARS, XDATA)
    model, pe = _fit_pool(ydata)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        noise_floor(model, XDATA, ydata, pe.x, np.array([1e-9]), margin=pe.margin,
                    alpha=pe.alpha, b=50, seed=4)
    assert not any("exceeds the noise-floor threshold" in str(w.message) for w in caught)


def test_seed_reproducibility():
    ydata = _fixture_func(TRUE_PARS, XDATA)
    model, pe = _fit_pool(ydata)
    cost_pool = np.array([0.05])

    out1 = noise_floor(model, XDATA, ydata, pe.x, cost_pool, margin=pe.margin, alpha=pe.alpha,
                        b=100, seed=7)
    out2 = noise_floor(model, XDATA, ydata, pe.x, cost_pool, margin=pe.margin, alpha=pe.alpha,
                        b=100, seed=7)
    np.testing.assert_array_equal(out1.cost, out2.cost)
    assert out1.threshold == out2.threshold


def test_different_seeds_give_different_draws():
    ydata = _fixture_func(TRUE_PARS, XDATA)
    model, pe = _fit_pool(ydata)
    cost_pool = np.array([0.05])

    out1 = noise_floor(model, XDATA, ydata, pe.x, cost_pool, margin=pe.margin, alpha=pe.alpha,
                        b=100, seed=7)
    out2 = noise_floor(model, XDATA, ydata, pe.x, cost_pool, margin=pe.margin, alpha=pe.alpha,
                        b=100, seed=8)
    assert not np.array_equal(out1.cost, out2.cost)


def test_degenerate_margin_raises():
    ydata = _fixture_func(TRUE_PARS, XDATA)
    model, pe = _fit_pool(ydata)
    with pytest.raises(ValueError, match="must not be 0"):
        noise_floor(model, XDATA, ydata, pe.x, pe.cost, margin=0.0, alpha=2.0, b=10)


def test_cumulative_wrong_length_raises():
    ydata = _fixture_func(TRUE_PARS, XDATA)
    model, pe = _fit_pool(ydata)
    with pytest.raises(ValueError, match="one entry per fitted output"):
        noise_floor(model, XDATA, ydata, pe.x, pe.cost, margin=pe.margin, alpha=pe.alpha, b=10,
                    cumulative=[False])


def test_costf_is_scale_invariant_under_joint_data_rescale():
    # gsua_costf/costf's own regulator normalization is provably scale-invariant under a joint
    # (ydata,yfunction) rescale by a positive constant for any margin/alpha (MSE numerator and
    # regulator both scale as scale**2 and cancel; the Pearson correlation term is scale-invariant
    # by construction) -- verified directly here via costf itself, not assumed. This is what makes
    # cost values comparable/rankable across a model's absolute cost units; it is a DIFFERENT
    # claim from noise_floor's own threshold, which is intentionally NOT scale-invariant (see the
    # module docstring) because it is calibrated against a physically meaningful count-noise model.
    model = _model()
    ydata = _fixture_func(TRUE_PARS, XDATA)
    pars = np.array([4.2, 0.35, 0.9])
    yf = model.evaluate(pars, XDATA)
    margin_internal = 1.1
    regulator = np.sum((ydata - ydata * margin_internal) ** 2, axis=1) / ydata.shape[1]
    res_original = costf(ydata, yf, regulator, alpha=2.0)

    scale = 37.0
    pars_scaled = pars.copy()
    pars_scaled[0] *= scale
    ydata_scaled = ydata * scale
    yf_scaled = model.evaluate(pars_scaled, XDATA)
    regulator_scaled = np.sum((ydata_scaled - ydata_scaled * margin_internal) ** 2, axis=1) / ydata_scaled.shape[1]
    res_scaled = costf(ydata_scaled, yf_scaled, regulator_scaled, alpha=2.0)

    assert res_scaled == pytest.approx(res_original, abs=1e-9)


def test_calibrated_band_covers_at_least_as_well_as_naive_rule():
    rng = np.random.default_rng(55)
    ydata0 = _fixture_func(TRUE_PARS, XDATA)
    ydata = ydata0 + np.vstack([
        rng.poisson(np.maximum(ydata0[0] * 0.3, 0)),
        rng.poisson(np.maximum(ydata0[1] * 0.02, 0)),
    ])

    model = _model()
    initial = TRUE_PARS[None, :] * rng.uniform(0.7, 1.3, size=(20, 3))
    initial = np.clip(initial, model.range[:, 0], model.range[:, 1])
    pe = parameter_estimation(model, XDATA, ydata, n=20, solver="minimize", margin=0.1,
                               initial_points=initial)

    out = noise_floor(model, XDATA, ydata, pe.x, pe.cost, margin=pe.margin, alpha=pe.alpha,
                       b=500, seed=30)

    def band_coverage(accepted_mask):
        curves = np.stack([model.evaluate(p, XDATA) for p in pe.x[accepted_mask]], axis=-1)
        lo = curves.min(axis=-1)
        hi = curves.max(axis=-1)
        in_band = (ydata >= lo) & (ydata <= hi)
        return in_band.mean()

    coverage_cal = band_coverage(out.accepted)

    lim_naive = int(np.sum(pe.cost < 1.5 * pe.cost.min()))
    order = np.argsort(pe.cost)
    accepted_naive = np.zeros_like(out.accepted)
    accepted_naive[order[:lim_naive]] = True
    coverage_naive = band_coverage(accepted_naive)

    assert coverage_cal >= coverage_naive
