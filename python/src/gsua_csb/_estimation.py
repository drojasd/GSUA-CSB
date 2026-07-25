"""Parameter estimation: the Python replacement for ``gsua_pe``.

MATLAB's ``gsua_pe`` dispatches across seven MATLAB optimizers (``lsqcurvefit``, ``lsqnonlin``,
``ga``, ``particleswarm``, ``patternsearch``, ``surrogateopt``, ``simulannealbnd``, ``fmincon``)
behind one ``'solver'`` string, repeats the estimation ``N`` times from different starting points
(by default drawn the same way as :func:`gsua_csb.design_matrix`), and returns every attempt sorted
by fit quality -- the standard defense against a nonlinear optimizer converging to a local rather
than global minimum. This module ports that shape onto four ``scipy.optimize`` backends chosen to
cover the same three problem classes MATLAB's seven span:

- ``"least_squares"`` (default) -- bounded Levenberg-Marquardt-style local refinement via
  ``scipy.optimize.least_squares`` (``method="trf"``). Replaces ``lsqcurvefit``/``lsqnonlin``.
- ``"minimize"`` -- bounded gradient-based local refinement via ``scipy.optimize.minimize``
  (default ``method="L-BFGS-B"``). Replaces ``fmincon``.
- ``"differential_evolution"`` -- population-based global search via
  ``scipy.optimize.differential_evolution``. Replaces ``ga``.
- ``"dual_annealing"`` -- simulated-annealing-style global search via
  ``scipy.optimize.dual_annealing``. Replaces ``simulannealbnd``.

Not ported: ``particleswarm``, ``patternsearch``, ``surrogateopt``, and ``MultiStart`` have no
close `scipy.optimize` equivalent (the closest analogues live in separate, less-maintained
third-party packages) -- flagged as a lower-priority follow-up rather than approximated poorly.
MATLAB's retry-on-exception loop (``err_counter``/``err_limit``, silently resampling the start
point up to 100 times after any error) is also not ported: a well-posed bounded scipy call raising
an exception is a real bug to surface, not routine noise to paper over by retrying blindly.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import differential_evolution, dual_annealing, least_squares, minimize

from ._costs import rcostf
from ._model import Model
from ._sampling import design_matrix

PESolver = Literal["least_squares", "minimize", "differential_evolution", "dual_annealing"]
_VALID_SOLVERS = {"least_squares", "minimize", "differential_evolution", "dual_annealing"}


@dataclass
class PEResult:
    """Result of :func:`parameter_estimation`.

    Attributes:
        names: Parameter names, length Np (``model.names``).
        x: (N, Np) estimated parameter sets, one row per attempt, sorted by ``cost`` ascending
            (``x[0]`` is the best fit found). Fixed parameters are carried through at
            ``model.nominal`` in every row, matching every other free-parameter-only routine in
            this package.
        cost: (N,) fit cost for each row of ``x``, sorted ascending. For ``solver="least_squares"``
            this is the summed squared residual; for the other solvers it is the objective
            actually minimized (MSE, or the margin-normalized :func:`gsua_csb.rcostf` cost if
            ``margin`` was given).
        solver: The solver used.
    """

    names: list[str]
    x: NDArray[np.float64]
    cost: NDArray[np.float64]
    solver: str


def parameter_estimation(
    model: Model,
    xdata: ArrayLike,
    ydata: ArrayLike,
    n: int = 1,
    solver: PESolver = "least_squares",
    initial_points: ArrayLike | None = None,
    margin: float = 0.0,
    alpha: float = 2.0,
    seed: int | np.random.Generator | None = None,
    solver_kwargs: dict | None = None,
) -> PEResult:
    """Fit ``model``'s free parameters to ``(xdata, ydata)``, from ``n`` starting points.

    Args:
        model: The model to fit. ``model.range`` supplies optimizer bounds; fixed parameters
            (``model.fixed``) are excluded from the search and held at ``model.nominal``.
        xdata: Points to evaluate the model at (passed to ``model.evaluate``).
        ydata: Reference data to fit against, same shape as ``model.evaluate(params, xdata)``.
        n: Number of independent estimation attempts (multistart). Each attempt gets its own
            starting point; results are returned sorted by fit quality so a converged-to-local-
            minimum attempt doesn't silently win.
        solver: Which `scipy.optimize` backend to use -- see the module docstring for the mapping
            to MATLAB's solvers.
        initial_points: (n, Np) starting points, one row per attempt. Defaults to
            :func:`gsua_csb.design_matrix(model, n, seed=seed)`, matching MATLAB's default
            behavior of generating a fresh design matrix when no ``ipoint`` is given.
        margin: For ``solver`` in ``{"minimize", "differential_evolution", "dual_annealing"}``:
            if ``0`` (default), the objective is plain MSE (mirrors MATLAB's ``immse`` default
            path). If nonzero, the objective is :func:`gsua_csb.rcostf` with this margin (mirrors
            MATLAB's ``gsua_costf``-based path) -- below 1 means within tolerance. Ignored for
            ``solver="least_squares"``, which always minimizes the raw residual.
        alpha: Exponent passed to :func:`gsua_csb.rcostf` when ``margin != 0``.
        seed: Seed for the default ``initial_points`` and for the stochastic global solvers
            (``differential_evolution``, ``dual_annealing``).
        solver_kwargs: Extra keyword arguments passed through to the underlying `scipy.optimize`
            call.

    Returns:
        A :class:`PEResult`.

    Raises:
        ValueError: If ``solver`` is not one of the four supported values.
    """
    if solver not in _VALID_SOLVERS:
        raise ValueError(f"Unknown solver: {solver!r}, expected one of {sorted(_VALID_SOLVERS)}")

    ydata = np.asarray(ydata, dtype=np.float64)
    solver_kwargs = dict(solver_kwargs or {})

    free_idx = np.where(~model.fixed)[0]
    lb = model.range[free_idx, 0]
    ub = model.range[free_idx, 1]

    if initial_points is None:
        initial_points = design_matrix(model, n, seed=seed)
    else:
        initial_points = np.asarray(initial_points, dtype=np.float64)
        if initial_points.shape[0] != n:
            raise ValueError(f"initial_points must have {n} rows (n={n}), got {initial_points.shape[0]}")

    def full_params(p_free: NDArray[np.float64]) -> NDArray[np.float64]:
        full = model.nominal.copy()
        full[free_idx] = p_free
        return full

    def residual(p_free: NDArray[np.float64]) -> NDArray[np.float64]:
        y = model.evaluate(full_params(p_free), xdata)
        return np.ravel(np.asarray(y, dtype=np.float64) - ydata)

    def cost(p_free: NDArray[np.float64]) -> float:
        y = model.evaluate(full_params(p_free), xdata)
        if margin == 0:
            return float(np.mean((ydata - np.asarray(y, dtype=np.float64)) ** 2))
        return rcostf(ydata, y, margin=margin, alpha=alpha)

    results_x = np.zeros((n, model.n_params))
    results_cost = np.full(n, np.inf)

    for i in range(n):
        x0_free = initial_points[i, free_idx]

        if solver == "least_squares":
            sol = least_squares(residual, x0_free, bounds=(lb, ub), **solver_kwargs)
            results_x[i] = full_params(sol.x)
            results_cost[i] = float(np.sum(sol.fun**2))
        elif solver == "minimize":
            sol = minimize(cost, x0_free, bounds=list(zip(lb, ub)), **solver_kwargs)
            results_x[i] = full_params(sol.x)
            results_cost[i] = float(sol.fun)
        elif solver == "differential_evolution":
            sol = differential_evolution(cost, bounds=list(zip(lb, ub)), seed=seed, **solver_kwargs)
            results_x[i] = full_params(sol.x)
            results_cost[i] = float(sol.fun)
        elif solver == "dual_annealing":
            sol = dual_annealing(cost, bounds=list(zip(lb, ub)), seed=seed, **solver_kwargs)
            results_x[i] = full_params(sol.x)
            results_cost[i] = float(sol.fun)

    order = np.argsort(results_cost)
    return PEResult(
        names=list(model.names), x=results_x[order], cost=results_cost[order], solver=solver
    )
