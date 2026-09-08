"""Profile-likelihood confidence intervals: the Python replacement for ``gsua_likelihood``.

For each parameter, steps away from ``model.nominal`` (independently, up and down) while every
*other* free parameter is re-optimized to best explain the data with that parameter held fixed at
the trial value, until the resulting likelihood-ratio statistic crosses the chi-squared threshold
implied by ``alpha``. This is an alternative to the sampling-based ranges from
:func:`gsua_csb.range_refinement`/:func:`gsua_csb.confidence_subcontour_box`, grounded directly in
the likelihood-ratio test rather than Monte Carlo coverage.

Two deliberate deviations from the MATLAB source:

- ``margin`` is offset by one relative to MATLAB's. **MATLAB's ``margin=1.1`` is this module's
  ``margin=0.1``**; both mean a 10% relative standard deviation.

  MATLAB is internally consistent about this, contrary to an earlier note here which claimed a
  unit mismatch and was wrong. Its threshold uses ``desv = (ydata*(margin-1))^2``
  (``gsua_likelihood.m:62``) and its inner re-optimization passes ``'margin', -margin`` to
  ``gsua_pe`` (``:117``), which adds one (``gsua_pe.m:81``) to get ``-(margin-1)``; being
  negative, that selects the Gaussian-NLL branch with ``desv = abs((ydata*margin)^2)``
  (``gsua_pe.m:177``) -- the same ``(ydata*(margin-1))^2``. Both ends agree.

  This module drops the offset so that ``margin`` is simply the relative standard deviation,
  matching how ``margin`` reads everywhere else in this package. Carrying a value across
  languages unchanged is the thing to avoid: MATLAB's ``margin=0.1`` implies a 90% relative std
  *and* makes ``gsua_pe``'s offset positive, which silently switches its inner refit to plain
  MSE instead of the likelihood -- so the threshold and the refit stop matching. Pass ``1+m`` in
  MATLAB where this module takes ``m``. With that correction the two implementations agree
  closely (on a 3-parameter pharmacokinetic fit, ``ka`` intervals of width 0.555 and 0.606).
- The inner re-optimization at each trial point uses a single local refit
  (``scipy.optimize.minimize``) from the current nominal point rather than MATLAB's
  multistart-and-take-the-best (``gsua_pe`` with ``reps`` random starts). Profile likelihood
  already requires many refits per parameter (one per bisection step, two directions, for every
  profiled parameter), so multistarting each one would multiply an already expensive computation;
  a single local refit from a good starting point (the previous trial's optimum) is the standard
  simplification used in most profile-likelihood implementations.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import minimize
from scipy.stats import chi2

from ._model import Model


def _gaussian_nll(ydata: NDArray[np.float64], y: NDArray[np.float64], margin: float) -> float:
    """Gaussian negative log-likelihood of ``y`` given ``ydata`` with std ``|ydata| * margin``."""
    desv = (ydata * margin) ** 2
    return float(np.nansum(np.log(2 * np.pi * desv) + (ydata - y) ** 2 / desv) / 2)


def _refit_nll(
    model: Model,
    xdata: NDArray[np.float64] | None,
    ydata: NDArray[np.float64],
    margin: float,
    fixed_idx: int,
    fixed_value: float,
    x0_free: NDArray[np.float64],
) -> tuple[float, NDArray[np.float64]]:
    """Re-optimize every free parameter except ``fixed_idx`` (held at ``fixed_value``); return
    ``(nll, best_free_params)``."""
    free_idx = np.where(~model.fixed)[0]
    free_idx = free_idx[free_idx != fixed_idx]

    def full_params(p_free: NDArray[np.float64]) -> NDArray[np.float64]:
        full = model.nominal.copy()
        full[fixed_idx] = fixed_value
        full[free_idx] = p_free
        return full

    if free_idx.size == 0:
        y = model.evaluate(full_params(np.array([])), xdata)
        return _gaussian_nll(ydata, np.asarray(y, dtype=np.float64), margin), x0_free

    lb = model.range[free_idx, 0]
    ub = model.range[free_idx, 1]

    def obj(p_free: NDArray[np.float64]) -> float:
        y = model.evaluate(full_params(p_free), xdata)
        return _gaussian_nll(ydata, np.asarray(y, dtype=np.float64), margin)

    sol = minimize(obj, x0_free, bounds=list(zip(lb, ub)))
    return float(sol.fun), sol.x


@dataclass
class ProfileLikelihoodResult:
    """Result of :func:`profile_likelihood`.

    Attributes:
        names: Parameter names, length Np (``model.names``).
        range: (Np, 2) profile-likelihood confidence interval per profiled parameter. Parameters
            not in ``params`` keep their original ``model.range`` row unchanged.
        threshold: The Gaussian NLL acceptance threshold (``nominal_cost + chi2.ppf(alpha, 1) / 2``)
            a trial point's refit cost must stay under.
        nominal_cost: The Gaussian NLL at ``model.nominal``.
    """

    names: list[str]
    range: NDArray[np.float64]
    threshold: float
    nominal_cost: float


def profile_likelihood(
    model: Model,
    xdata: ArrayLike,
    ydata: ArrayLike | None = None,
    alpha: float = 0.95,
    step: float = 0.1,
    margin: float = 0.1,
    tol_obj: float = 1e-3,
    tol_step: float = 1e-3,
    max_iter: int = 50,
    params: Sequence[int] | None = None,
) -> ProfileLikelihoodResult:
    """Compute a profile-likelihood confidence interval for each of ``model``'s free parameters.

    Args:
        model: The model to profile. ``model.nominal`` is the reference point;
            ``model.range`` bounds both the search and the inner refits.
        xdata: Points to evaluate the model at.
        ydata: Reference data. Defaults to ``model.evaluate(model.nominal, xdata)`` (profiling
            around a perfect fit to the nominal model, useful for a purely structural
            identifiability check with no real data).
        alpha: Confidence level for the chi-squared threshold (e.g. 0.95).
        step: Initial relative step size (as a fraction of the nominal-to-bound distance) used to
            perturb each parameter away from its nominal value.
        margin: Relative standard deviation for the Gaussian likelihood weighting (see module
            docstring for why this uses one consistent basis, unlike MATLAB).
        tol_obj: Convergence tolerance on the likelihood-ratio statistic (bisection stop).
        tol_step: Convergence tolerance on the step size, relative to the initial nominal-to-bound
            distance.
        max_iter: Maximum number of bisection iterations per parameter bound.
        params: Indices of the parameters to profile. Defaults to every free (non-fixed) parameter.

    Returns:
        A :class:`ProfileLikelihoodResult`.
    """
    d = np.asarray(xdata, dtype=np.float64)
    y_nom = np.asarray(model.evaluate(model.nominal, d), dtype=np.float64)
    ydata_arr = y_nom if ydata is None else np.asarray(ydata, dtype=np.float64)

    lam = chi2.ppf(alpha, df=1) / 2
    nominal_cost = _gaussian_nll(ydata_arr, y_nom, margin)
    threshold = nominal_cost + lam

    profile_idx = np.where(~model.fixed)[0] if params is None else np.asarray(params, dtype=int)
    new_range = model.range.copy()

    for i in profile_idx:
        for direction in (-1, +1):
            j = 0 if direction < 0 else 1
            new_range[i, j] = _bisect_profile_bound(
                model, d, ydata_arr, margin, i, direction, step, tol_obj, tol_step, max_iter, threshold
            )

    return ProfileLikelihoodResult(
        names=list(model.names), range=new_range, threshold=threshold, nominal_cost=nominal_cost
    )


def _bisect_profile_bound(
    model: Model,
    xdata: NDArray[np.float64],
    ydata: NDArray[np.float64],
    margin: float,
    i: int,
    direction: int,
    step: float,
    tol_obj: float,
    tol_step: float,
    max_iter: int,
    threshold: float,
) -> float:
    nom_i = float(model.nominal[i])
    orig_bound = float(model.range[i, 0 if direction < 0 else 1])

    step_par = abs(nom_i - orig_bound) * step
    if step_par == 0:
        step_par = abs(nom_i * step) or 1.0

    free_idx = np.where(~model.fixed)[0]
    free_idx = free_idx[free_idx != i]
    x0_free = model.nominal[free_idx].copy()

    trial = nom_i + direction * step_par
    par_hist = [nom_i, trial]
    obj_hist: list[float] = []
    bisecting = False
    bracket: list[float] | None = None
    counter = 0
    bound = orig_bound

    while True:
        if (direction < 0 and trial < model.range[i, 0]) or (direction > 0 and trial > model.range[i, 1]):
            bound = orig_bound
            break

        cost, x0_free = _refit_nll(model, xdata, ydata, margin, i, trial, x0_free)
        obj = threshold - cost
        obj_hist.append(obj)

        if abs(obj) < tol_obj:
            bound = trial
            break

        if not bisecting:
            if obj > 0:
                trial = trial + direction * step_par
                par_hist.append(trial)
            else:
                trial = (par_hist[-1] + par_hist[-2]) / 2
                par_hist.append(trial)
                bisecting = True
                bracket = [par_hist[-3], par_hist[-2]]
        else:
            counter += 1
            if obj_hist[-1] * obj_hist[-2] > 0:
                bracket = [bracket[0], trial]
            else:
                bracket = [par_hist[-2], trial]
            trial = (bracket[0] + bracket[1]) / 2
            par_hist.append(trial)

        if counter > max_iter:
            bound = trial
            break
        if abs(par_hist[-1] - par_hist[-2]) < tol_step * abs(nom_i - orig_bound):
            bound = trial
            break

    return bound
