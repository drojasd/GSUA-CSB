import numpy as np
import pytest

sp = pytest.importorskip("sympy")

from gsua_csb import SymbolicODEModel


def test_exponential_decay_matches_analytic_solution():
    # solve_ivp's default tolerance (rtol=1e-3) only guarantees ~0.1% accuracy -- tighten it here
    # so the test actually checks correctness rather than default numerical-integration slop.
    t, y, k = sp.symbols("t y k")
    model = SymbolicODEModel(
        odes=[-k * y],
        state_vars=[y],
        t=t,
        params=[k],
        domain=np.linspace(0, 5, 50),
        names=["y0", "k"],
        range=np.array([[1.0, 1.0], [0.5, 0.5]]),
        solver_kwargs={"rtol": 1e-10, "atol": 1e-12},
    )
    out = model.evaluate(np.array([1.0, 0.5]))
    expected = np.exp(-0.5 * model.domain)
    np.testing.assert_allclose(out[0], expected, rtol=1e-8)


def test_lotka_volterra_stays_positive_and_correct_shape():
    t = sp.Symbol("t")
    x, y = sp.symbols("x y")
    a, b, c, d = sp.symbols("a b c d")
    model = SymbolicODEModel(
        odes=[a * x - b * x * y, -c * y + d * x * y],
        state_vars=[x, y],
        t=t,
        params=[a, b, c, d],
        domain=np.linspace(0, 10, 50),
        names=["x0", "y0", "a", "b", "c", "d"],
        range=np.tile([[1.0, 1.0]], (6, 1)),
    )
    params = np.array([1.0, 1.0, 1.1, 0.4, 0.4, 0.1])
    out = model.evaluate(params)
    assert out.shape == (2, 50)
    assert np.all(out >= 0)  # predator-prey populations should stay non-negative here


def test_evaluate_batch_default_loop_over_independent_ivps():
    t, y, k = sp.symbols("t y k")
    model = SymbolicODEModel(
        odes=[-k * y],
        state_vars=[y],
        t=t,
        params=[k],
        domain=np.linspace(0, 5, 20),
        solver_kwargs={"rtol": 1e-10, "atol": 1e-12},
    )
    batch = np.array([[1.0, 0.5], [2.0, 0.5], [1.0, 1.0]])
    out = model.evaluate_batch(batch)
    assert out.shape == (3, 1, 20)
    np.testing.assert_allclose(out[0, 0], np.exp(-0.5 * model.domain), rtol=1e-8)
    np.testing.assert_allclose(out[1, 0], 2 * np.exp(-0.5 * model.domain), rtol=1e-8)


def test_from_bounds_percent_method():
    t, y, k = sp.symbols("t y k")
    model = SymbolicODEModel.from_bounds(
        odes=[-k * y],
        state_vars=[y],
        t=t,
        params=[k],
        domain=np.linspace(0, 5, 10),
        values=[1.0, 0.5],
        spread=[10.0, 20.0],
        method="percent",
    )
    np.testing.assert_allclose(model.range, [[0.9, 1.1], [0.4, 0.6]])


def test_cross_validated_against_matlab_gsua_dpmat(matlab_lotka_volterra_reference):
    """Same ODE system, same parameters, evaluated in both MATLAB (gsua_dpmat/gsua_deval) and
    here -- compares against a reference captured directly from the real MATLAB toolbox.

    The reference was captured at tight solver tolerance (RelTol=1e-10) deliberately: at default
    tolerance, MATLAB's ode45 and scipy's RK45 take different numerical paths through the
    Lotka-Volterra system's sensitive peak region and disagree by ~5%, even though both are
    correct implementations converging to the same true solution -- verified separately by
    tightening both sides and confirming they then agree to ~9-10 significant figures. Using a
    tight reference (and tight solver_kwargs here) tests actual correctness instead of being
    sensitive to each library's default step-size heuristics.
    """
    t = sp.Symbol("t")
    x, y = sp.symbols("x y")
    a, b, c, d = sp.symbols("a b c d")
    model = SymbolicODEModel(
        odes=[a * x - b * x * y, -c * y + d * x * y],
        state_vars=[x, y],
        t=t,
        params=[a, b, c, d],
        domain=matlab_lotka_volterra_reference["t"],
        names=["x0", "y0", "a", "b", "c", "d"],
        range=np.tile([[1.0, 1.0]], (6, 1)),
        solver_kwargs={"rtol": 1e-10, "atol": 1e-12},
    )
    params = matlab_lotka_volterra_reference["params"]
    out = model.evaluate(params)
    np.testing.assert_allclose(out[0], matlab_lotka_volterra_reference["x"], rtol=1e-6, atol=1e-6)
    np.testing.assert_allclose(out[1], matlab_lotka_volterra_reference["y"], rtol=1e-6, atol=1e-6)
