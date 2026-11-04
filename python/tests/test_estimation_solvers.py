"""Smoke test: every parameter_estimation solver runs and returns a finite cost.

Python parity for the MATLAB gsua_pe solver-dispatch fix (#12). scipy's ``dual_annealing`` and
``differential_evolution`` take ``bounds`` (not a mandatory ``x0``), so the MATLAB
``simulannealbnd``/``patternsearch`` x0 bug has no analogue here -- this just guards the four
dispatch paths against regressions.
"""
import numpy as np
import pytest

import gsua_csb as g


def _model():
    def decay(p, x):
        return p[0] * np.exp(-p[1] * np.asarray(x, dtype=float))

    return g.UserFunctionModel(
        decay, names=["a", "b"], range=np.array([[0.1, 6.0], [0.05, 3.0]]),
        domain=np.linspace(0, 2, 12), output_names=["y"],
    )


@pytest.mark.parametrize(
    "solver", ["least_squares", "minimize", "differential_evolution", "dual_annealing"]
)
def test_solver_runs_without_error(solver):
    model = _model()
    x = np.linspace(0, 2, 12)
    y = 3.0 * np.exp(-1.0 * x)
    pe = g.parameter_estimation(model, x, y, n=1, solver=solver, seed=0)
    assert pe.x.shape == (1, 2)
    assert np.isfinite(pe.cost[0])


def test_unknown_solver_raises():
    model = _model()
    x = np.linspace(0, 2, 12)
    y = 3.0 * np.exp(-1.0 * x)
    with pytest.raises(ValueError, match="Unknown solver"):
        g.parameter_estimation(model, x, y, n=1, solver="nope")
