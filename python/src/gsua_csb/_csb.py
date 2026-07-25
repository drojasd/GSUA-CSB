"""Range refinement and the Confidence Sub-contour Box: the Python replacement for
``gsua_oatr``/``gsua_oatr2``/``gsua_csb`` -- the toolbox's namesake algorithm.

:func:`range_refinement` (MATLAB: ``gsua_oatr``/``gsua_oatr2``) independently searches, for each
free parameter, how far it can move from ``model.nominal`` before the model output diverges from
the nominal curve by more than a fractional tolerance -- a one-at-a-time range expansion/reduction.
MATLAB implements the boundary search as a hand-rolled step-doubling/halving loop with an internal
retry counter; this is a textbook bracket-and-refine root-find, so it is reimplemented here with
``scipy.optimize.brentq`` instead of translating that loop literally. ``gsua_oatr`` (additive step)
and ``gsua_oatr2`` (multiplicative step) are two step-size strategies for the same underlying
search and are collapsed into the one robust implementation here.

:func:`confidence_subcontour_box` (MATLAB: ``gsua_csb``) iteratively narrows the *whole* parameter
box (not one parameter at a time): sample a design matrix from the current box, score every run's
fit against a reference output, and if too few runs are "behavioral" (fit within tolerance), trim
each free parameter's marginal distribution among the behavioral/best-fit runs back toward
uniformity (via repeated Kolmogorov-Smirnov-vs-uniform testing) and shrink the box to that trimmed
range. Repeat until enough runs are behavioral or a repetition budget is spent.

Not ported (documented gap): MATLAB's ``gsua_csb`` has a recursive escape hatch -- when the
per-iteration ``select`` fraction shrinks below a threshold without the box narrowing further, it
recursively re-invokes itself with a larger sample size ``N``, up to a ``stretch`` recursion depth.
This module instead simply continues iterating with the current sample size and reports whether the
stop condition was reached (``CSBResult.converged``) -- for a problem that needs the recursive
resample-larger escape hatch to converge, increase ``n`` directly rather than relying on an
automatic escalation.
"""

from __future__ import annotations

import copy
import warnings
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import brentq
from scipy.stats import kstest
from scipy.stats import kurtosis as _kurtosis
from scipy.stats import skew as _skew

from ._costs import costf_multi, rcostf
from ._evalutils import nominal_output
from ._model import Model
from ._sampling import design_matrix


@dataclass
class RangeRefinementResult:
    """Result of :func:`range_refinement`.

    Attributes:
        names: Parameter names, length Np (``model.names``).
        range: (Np, 2) refined [lower, upper] bounds.
    """

    names: list[str]
    range: NDArray[np.float64]


def _search_boundary(
    cost_at, nom: float, sign: int, step0: float, max_expansions: int
) -> float:
    """Bracket-then-refine the point where ``cost_at`` crosses zero, walking outward from ``nom``."""
    a, fa = nom, cost_at(nom)
    step = step0
    b = nom + sign * step
    fb = cost_at(b)
    n = 0
    while fb < 0 and n < max_expansions:
        step *= 2
        a, fa = b, fb
        b = nom + sign * step
        fb = cost_at(b)
        n += 1
    if fb < 0:
        return b  # tolerance boundary never reached within budget -- unconstrained in this direction
    lo, hi = (a, b) if sign > 0 else (b, a)
    try:
        return brentq(cost_at, lo, hi, xtol=1e-8)
    except ValueError:
        return b


def range_refinement(
    model: Model,
    xdata: ArrayLike | None = None,
    lim: float = 0.3,
    correct: bool = False,
    max_expansions: int = 30,
) -> RangeRefinementResult:
    """Independently expand/reduce each free parameter's range around ``model.nominal``.

    For each free parameter, searches outward in both directions from ``model.nominal`` for the
    point where :func:`gsua_csb.rcostf` (margin ``1 + lim``) crosses 1 -- the boundary past which
    the model output has diverged from the nominal curve by more than the ``lim`` fractional
    tolerance.

    Args:
        model: The model to refine. ``model.nominal`` is the reference point; ``model.range``
            supplies the initial step size and (if ``correct=True``) the outer clipping bound.
        xdata: Points to evaluate at. Defaults to ``model.domain``.
        lim: Fractional tolerance defining the acceptance boundary. Default 0.3 (30%).
        correct: If True, clip the refined range to ``model.range`` where it would otherwise
            extend beyond it.
        max_expansions: Maximum number of step-doubling expansions per direction before giving up
            and treating that direction as unconstrained (returning the farthest point tried).

    Returns:
        A :class:`RangeRefinementResult`.
    """
    d = model.domain if xdata is None else np.asarray(xdata, dtype=np.float64)
    y_nom = nominal_output(model, d)
    free_idx = np.where(~model.fixed)[0]
    new_range = model.range.copy()

    for i in free_idx:
        nom = float(model.nominal[i])
        orig_lo, orig_hi = model.range[i]

        def cost_at(value: float, i: int = i) -> float:
            params = model.nominal.copy()
            params[i] = value
            try:
                y = model.evaluate(params, d)
            except Exception:
                return np.inf
            return rcostf(y_nom, y, margin=1 + lim) - 1.0

        step_hi = abs(nom - orig_hi) * 0.1 or abs(nom) * 0.1 or 1.0
        step_lo = abs(nom - orig_lo) * 0.1 or abs(nom) * 0.1 or 1.0
        hi = _search_boundary(cost_at, nom, +1, step_hi, max_expansions)
        lo = _search_boundary(cost_at, nom, -1, step_lo, max_expansions)
        if lo > hi:
            lo, hi = hi, lo
        if correct:
            lo = max(lo, orig_lo)
            hi = min(hi, orig_hi)
        new_range[i] = [lo, hi]

    return RangeRefinementResult(names=list(model.names), range=new_range)


def _trim_to_uniform(
    data: NDArray[np.float64], alpha: float = 0.1, step: float = 0.05, max_iter: int = 40
) -> tuple[float, float]:
    """Iteratively trim the tails of ``data`` (already normalized to [0, 1]) toward uniformity.

    Mirrors ``gsua_csb``'s inner loop: while a KS test rejects uniformity and the (Pearson,
    non-excess) kurtosis exceeds 1.79 (a uniform distribution's kurtosis is exactly 1.8), trim 5%
    off whichever tail(s) the skewness indicates, stopping early if the KS p-value stops improving.

    Returns:
        ``(lo, hi)`` bounds in [0, 1] relative to the original (pre-trim) range.
    """
    lo, hi = 0.0, 1.0
    if data.size < 4:
        return lo, hi

    def _norm(x: NDArray[np.float64]) -> NDArray[np.float64]:
        rng = x.max() - x.min()
        return (x - x.min()) / rng if rng > 0 else np.zeros_like(x)

    mask = (data >= lo) & (data <= hi)
    x = data[mask]
    _, pval = kstest(_norm(x), "uniform")
    kurt = _kurtosis(x, fisher=False)
    p_hist: list[float] = []

    it = 0
    while pval < alpha and kurt > 1.79 and it < max_iter:
        sk = _skew(x)
        width = hi - lo
        if abs(sk) < 0.1:
            lo, hi = lo + width * step, hi - width * step
        elif sk > 0:
            hi = hi - width * step
        else:
            lo = lo + width * step

        mask = (data >= lo) & (data <= hi)
        x = data[mask]
        if x.size < 4:
            break
        _, pval = kstest(_norm(x), "uniform")
        kurt = _kurtosis(x, fisher=False)
        p_hist.append(pval)
        if len(p_hist) > 2:
            prev = p_hist[-2]
            if prev != 0 and abs((prev - p_hist[-1]) / prev) < 2:
                break
            if prev > p_hist[-1]:
                break
        it += 1

    return lo, hi


@dataclass
class CSBResult:
    """Result of :func:`confidence_subcontour_box`.

    Attributes:
        names: Parameter names, length Np (``model.names``).
        range: (Np, 2) final refined [lower, upper] bounds.
        range_history: (n_iterations + 1, Np, 2) range at each iteration, ``range_history[0]`` is
            the starting range.
        behavioral_fraction: (n_iterations,) fraction of the ``n`` sampled runs within tolerance
            (cost < 1) at each iteration.
        converged: True if ``behavioral_fraction`` reached ``stop`` within ``reps`` iterations.
    """

    names: list[str]
    range: NDArray[np.float64]
    range_history: NDArray[np.float64]
    behavioral_fraction: NDArray[np.float64]
    converged: bool


def confidence_subcontour_box(
    model: Model,
    n: int | None = None,
    xdata: ArrayLike | None = None,
    y_exp: ArrayLike | None = None,
    reps: int = 100,
    select: float = 0.5,
    lim: float = 0.3,
    stop: float = 0.95,
    switch_fraction: float = 0.75,
    protect: bool = False,
    correct: bool = True,
    seed: int | None = None,
) -> CSBResult:
    """Iteratively narrow ``model``'s parameter box to the region consistent with a reference output.

    Args:
        model: The model to refine. ``model.range`` is the starting box; ``model.nominal`` is used
            as the known-good point when ``protect=True``.
        n: Sample size per iteration. Defaults to ``max(10 * Np, 10 * Np)`` (MATLAB's minimum);
            raised to that minimum if given smaller.
        xdata: Points to evaluate at. Defaults to ``model.domain``.
        y_exp: Reference output to score against. Defaults to ``model.evaluate(model.nominal, xdata)``.
        reps: Maximum number of narrowing iterations.
        select: Fraction of the sample to treat as "best fit" when fewer than
            ``switch_fraction * n`` runs are behavioral (cost < 1).
        lim: Margin passed to :func:`gsua_csb.costf_multi` (fractional tolerance).
        stop: Target behavioral fraction; iteration stops early once reached.
        switch_fraction: Behavioral-count threshold (as a fraction of ``n``) above which every
            behavioral run is used instead of just the top ``select`` fraction.
        protect: If True, keep ``model.nominal`` inside the narrowed range for every free
            parameter, shifting the range by 10% of its width if narrowing would exclude it.
        correct: If True (default), clip the narrowed range to never exceed ``model.range``.
        seed: Seed for the per-iteration design matrix sampling.

    Returns:
        A :class:`CSBResult`.
    """
    d = model.domain if xdata is None else np.asarray(xdata, dtype=np.float64)
    y_nom = nominal_output(model, d) if y_exp is None else np.atleast_1d(np.asarray(y_exp, dtype=np.float64))

    free_idx = np.where(~model.fixed)[0]
    Np = model.n_params
    if n is None or n < Np * 10:
        n = Np * 10

    original_range = model.range.copy()
    current_range = model.range.copy()
    history = [current_range.copy()]
    frac_history = []
    converged = False

    rng = np.random.default_rng(seed)

    for _ in range(reps):
        temp_model = copy.copy(model)
        temp_model.range = current_range
        M = design_matrix(temp_model, n, method="latin_hypercube", seed=rng)
        Y = model.evaluate_batch(M, d)
        J = costf_multi(Y, y_nom, margin=lim, alpha=1.0)

        behavioral = J < 1
        n_behavioral = int(behavioral.sum())
        frac = n_behavioral / n
        frac_history.append(frac)

        if frac >= stop:
            converged = True
            break

        order = np.argsort(J)
        if n_behavioral > int(np.ceil(switch_fraction * n)):
            best = M[order[:n_behavioral]]
        else:
            best = M[order[: int(np.ceil(select * n))]]

        new_range = current_range.copy()
        for i in free_idx:
            lo0, hi0 = current_range[i]
            width0 = hi0 - lo0
            if width0 <= 0:
                continue
            norm_col = (best[:, i] - lo0) / width0
            lo_n, hi_n = _trim_to_uniform(norm_col)
            new_lo = lo_n * width0 + lo0
            new_hi = hi_n * width0 + lo0

            if protect:
                nom_i = model.nominal[i]
                w = new_hi - new_lo
                if new_lo >= nom_i:
                    new_lo = nom_i - w * 0.1
                    new_hi = new_hi - w * 0.1
                elif new_hi <= nom_i:
                    new_hi = nom_i + w * 0.1
                    new_lo = new_lo + w * 0.1

            if correct:
                new_lo = max(new_lo, original_range[i, 0])
                new_hi = min(new_hi, original_range[i, 1])

            new_range[i] = [new_lo, new_hi]

        current_range = new_range
        history.append(current_range.copy())

    if not converged:
        warnings.warn(
            f"confidence_subcontour_box did not reach the stop={stop} behavioral fraction "
            f"within reps={reps} iterations (last fraction: {frac_history[-1] if frac_history else 0.0})",
            stacklevel=2,
        )

    return CSBResult(
        names=list(model.names),
        range=current_range,
        range_history=np.stack(history),
        behavioral_fraction=np.array(frac_history),
        converged=converged,
    )
