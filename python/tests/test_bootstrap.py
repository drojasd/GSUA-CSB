"""Parametric-bootstrap ensemble (parametric_bootstrap / gsua_bootstrap).

Refitting the same data from many starts collapses to one optimum; the parametric bootstrap must
instead reveal the wide S0 uncertainty the single best fit hides, on the discrete-SIR case, and be
reproducible by seed.
"""
import numpy as np
import pytest

import gsua_csb as g


def _sir(factores, t):
    t = np.asarray(t, dtype=float)
    S0, I0, R0, beta, gamma, alpha = factores
    tmin, tmax = 1, int(np.max(t))
    nsteps = tmax - tmin + 1
    S = np.empty(nsteps); I = np.empty(nsteps); R = np.empty(nsteps)
    S[0], I[0], R[0] = S0, I0, R0
    for k in range(nsteps - 1):
        Nt = S[k] + I[k] + R[k]
        P = 1.0 - np.exp(-beta * I[k] / Nt)
        S[k + 1] = S[k] * (1 - P) + alpha * R[k]
        I[k + 1] = S[k] * P + I[k] * (1 - gamma)
        R[k + 1] = I[k] * gamma + R[k] * (1 - alpha)
    grid = np.arange(tmin, tmax + 1)
    return np.vstack([np.interp(t, grid, S), np.interp(t, grid, I), np.interp(t, grid, R)])


def _sir_setup():
    names = ["S0", "I0", "R0", "beta", "gamma", "alpha"]
    ranges = np.array([[1e3, 1e4], [0, 10], [10, 100], [0, 1], [0, 1], [0, 1]], float)
    truth = np.array([5000, 5, 50, 0.5, 0.2, 0.01])
    xdata = np.arange(1, 51)
    model = g.UserFunctionModel(_sir, names=names, range=ranges, domain=[1, 200],
                                output_names=["S", "I", "R"])
    model.set_output("I")
    yI = np.asarray(model.evaluate(truth, xdata), float).ravel()
    rng = np.random.default_rng(1)
    k = 10
    ydata = rng.negative_binomial(k, k / (k + yI)).astype(float)
    return model, xdata, ydata, truth


def test_bootstrap_shape_and_reproducibility():
    model, xdata, ydata, _ = _sir_setup()
    pe = g.parameter_estimation(model, xdata, ydata, n=10, solver="minimize", margin=0.3, seed=2)
    kw = dict(b=6, n=4, margin=pe.margin, alpha=pe.alpha, solver="minimize", seed=7, progress=False)
    b1 = g.parametric_bootstrap(model, xdata, ydata, pe.x, pe.cost, **kw)
    b2 = g.parametric_bootstrap(model, xdata, ydata, pe.x, pe.cost, **kw)
    assert b1.x.shape == (6, 6)
    assert b1.cost.shape == (6,)
    np.testing.assert_array_equal(b1.x, b2.x)   # a fixed seed reproduces the ensemble


def test_bootstrap_reveals_s0_nonidentifiability():
    model, xdata, ydata, _ = _sir_setup()
    pe = g.parameter_estimation(model, xdata, ydata, n=12, solver="minimize", margin=0.3, seed=2)
    bs = g.parametric_bootstrap(model, xdata, ydata, pe.x, pe.cost, b=10, n=6,
                                margin=pe.margin, alpha=pe.alpha, solver="minimize", seed=1,
                                progress=False)
    assert np.all(np.isfinite(bs.cost))
    s0 = bs.x[:, 0]
    assert s0.max() - s0.min() > 1000   # bootstrap exposes the S0 spread multistart hides


def test_bootstrap_bad_cumulative_raises():
    model, xdata, ydata, _ = _sir_setup()
    est = np.tile(model.nominal, (3, 1))
    cost = np.array([1.0, 2.0, 3.0])
    with pytest.raises(ValueError, match="one entry per fitted output"):
        g.parametric_bootstrap(model, xdata, ydata, est, cost, b=3, n=2, margin=0.3, alpha=2.0,
                               cumulative=[False, False], progress=False)
