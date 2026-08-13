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
        margin: The ``margin`` this call was made with -- recorded so a downstream noise-floor
            calibration (:func:`gsua_csb.noise_floor`) can recompute the identical ``rcostf``
            regulator without the caller having to re-specify (and risk mismatching) it. Only
            meaningful (and only comparable to a fresh ``rcostf`` bootstrap) when ``cost`` actually
            came from the ``rcostf``-scored path -- i.e. ``solver`` in ``{"minimize",
            "differential_evolution", "dual_annealing"}`` and ``margin != 0``; ``cost`` is raw
            squared-residual for ``solver="least_squares"`` regardless of this value.
        alpha: The ``alpha`` this call was made with -- same recovery purpose as ``margin``.
    """

    names: list[str]
    x: NDArray[np.float64]
    cost: NDArray[np.float64]
    solver: str
    margin: float
    alpha: float


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
            (``model.fixed``) are excluded from the search and held at ``model.nominal``. A free
            parameter with ``model.log_scale`` set is searched in log10 space internally (see
            :attr:`Model.log_scale`) -- transparent to the caller: ``initial_points`` (if given) and
            every entry of the returned ``PEResult.x`` are always in natural/linear units.
        xdata: Points to evaluate the model at (passed to ``model.evaluate``).
        ydata: Reference data to fit against, same shape as ``model.evaluate(params, xdata)``.
        n: Number of independent estimation attempts (multistart). Each attempt gets its own
            starting point; results are returned sorted by fit quality so a converged-to-local-
            minimum attempt doesn't silently win.
        solver: Which `scipy.optimize` backend to use -- see the module docstring for the mapping
            to MATLAB's solvers.
        initial_points: (n, Np) starting points (natural/linear units, regardless of
            ``model.log_scale``), one row per attempt. Defaults to
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
        ValueError: If ``solver`` is not one of the four supported values, or a free,
            ``log_scale`` parameter has a lower bound that isn't strictly positive.
    """
    if solver not in _VALID_SOLVERS:
        raise ValueError(f"Unknown solver: {solver!r}, expected one of {sorted(_VALID_SOLVERS)}")

    ydata = np.asarray(ydata, dtype=np.float64)
    solver_kwargs = dict(solver_kwargs or {})

    free_idx = np.where(~model.fixed)[0]
    lb_natural = model.range[free_idx, 0]
    ub_natural = model.range[free_idx, 1]

    # log_scale (see Model.log_scale) means this free parameter is searched in log10 space, not
    # linear -- essential for a parameter spanning many orders of magnitude (e.g. a PEtab
    # parameterScale="log10" rate constant): a linear-space local optimizer cannot traverse that
    # many orders of magnitude in a reasonable number of steps, regardless of starting point.
    log_scale = np.asarray(getattr(model, "log_scale", np.zeros(model.n_params, dtype=bool)), dtype=bool)
    free_log = log_scale[free_idx]
    if np.any(free_log):
        bad = free_idx[free_log][lb_natural[free_log] <= 0]
        if bad.size:
            raise ValueError(
                "log_scale parameters must have a strictly positive lower bound; failed for "
                f"{[model.names[i] for i in bad]}"
            )
    with np.errstate(divide="ignore", invalid="ignore"):
        # np.where evaluates both branches eagerly, so log10 of a non-log-scale bound (which may
        # legitimately be <= 0) still runs -- its result is simply never selected.
        lb = np.where(free_log, np.log10(lb_natural), lb_natural)
        ub = np.where(free_log, np.log10(ub_natural), ub_natural)

    if initial_points is None:
        initial_points = design_matrix(model, n, seed=seed)
    else:
        initial_points = np.asarray(initial_points, dtype=np.float64)
        if initial_points.shape[0] != n:
            raise ValueError(f"initial_points must have {n} rows (n={n}), got {initial_points.shape[0]}")

    def to_search_space(x0_free_natural: NDArray[np.float64]) -> NDArray[np.float64]:
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(free_log, np.log10(x0_free_natural), x0_free_natural)

    def full_params(p_search: NDArray[np.float64]) -> NDArray[np.float64]:
        p_free_natural = np.where(free_log, 10.0**p_search, p_search)
        full = model.nominal.copy()
        full[free_idx] = p_free_natural
        return full

    def residual(p_search: NDArray[np.float64]) -> NDArray[np.float64]:
        y = model.evaluate(full_params(p_search), xdata)
        return np.ravel(np.asarray(y, dtype=np.float64) - ydata)

    def cost(p_search: NDArray[np.float64]) -> float:
        y = model.evaluate(full_params(p_search), xdata)
        if margin == 0:
            return float(np.mean((ydata - np.asarray(y, dtype=np.float64)) ** 2))
        return rcostf(ydata, y, margin=margin, alpha=alpha)

    results_x = np.zeros((n, model.n_params))
    results_cost = np.full(n, np.inf)

    for i in range(n):
        x0_free = to_search_space(initial_points[i, free_idx])

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
        names=list(model.names), x=results_x[order], cost=results_cost[order], solver=solver,
        margin=margin, alpha=alpha,
    )
