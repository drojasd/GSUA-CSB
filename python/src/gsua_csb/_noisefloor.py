"""Noise-calibrated fit-acceptance threshold: the Python replacement for
``res < 1.5*res(1)``/``lims = sum(res < 1.5*res(1))``.

MATLAB equivalent: ``gsua_noisefloor`` (no prior MATLAB toolbox equivalent existed either -- this
is new capability, ported from Python back to MATLAB and Python together in the same pass).

That naive rule asks "how close is this to the best fit found", not "is this statistically
distinguishable from a perfect model given how noisy the data actually is" -- the former gets
*narrower* as the fit improves (backwards: a calibration band should not shrink because the
optimizer got lucky), and has no mechanism to guarantee the accepted band actually contains the
data it was fit to.

Method (parametric bootstrap of the noise floor):
    1. Take the best pool member (by cost) as theta*, simulate it to get fitted curves f_i(t).
    2. Estimate a pooled observation-noise model from theta*'s own residuals.
    3. Generate ``b`` synthetic datasets from ``f`` under that noise model.
    4. Score the TRUE model ``f`` against each synthetic dataset with :func:`gsua_csb.costf`,
       recomputing the regulator from the synthetic data exactly as :func:`gsua_csb.rcostf` does
       inside :func:`gsua_csb.parameter_estimation` -- so the resulting cost distribution is
       directly comparable to ``cost``.
    5. That distribution is the cost a perfect model incurs from noise alone; its upper quantile
       (default 95%) is the largest cost still consistent with being correct.
    6. Accept every pool member with ``cost < threshold``.

Positioning: this is a goodness-of-fit BAND-SIZING CALIBRATION, not a formal confidence region --
it does not replace :func:`gsua_csb.profile_likelihood` or :func:`gsua_csb.confidence_subcontour_box`,
it is complementary and cheap (no refitting). It is a PARAMETRIC bootstrap: the noise floor is
conditional on theta* being approximately correct -- under model misspecification the floor is
optimistic. The quantile is a user choice with real consequences; 95% is a convention, not a
derivation.

A genuine nuance worth knowing before relying on this: the calibrated *threshold* is intentionally
**not** invariant under an arbitrary rescale of the data's units, because it is calibrated against
a physically meaningful count-noise model (``var ~ mu``, tied to the data's actual units -- a real
observation of 500 counts has real Poisson-like noise of about sqrt(500), which is not preserved by
relabeling it "50000 sub-units"). What *is* scale-invariant is :func:`gsua_csb.costf`'s own
regulator normalization (verified directly, not assumed) -- that is what makes ``cost`` values
comparable/rankable across pool members regardless of a model's absolute cost units, independent of
this module. The noise floor exists to be tied to something *real* (the data's own noise), unlike
``1.5*cost.min()``, which is tied to nothing but an arbitrary multiplier on whatever the optimizer's
raw cost happened to be.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ._costs import costf
from ._model import Model

NoiseModel = Literal["poisson", "quasipoisson", "nbglobal"]
_MODELS: tuple[NoiseModel, ...] = ("poisson", "quasipoisson", "nbglobal")


@dataclass
class NoiseFloorByModel:
    """One noise model's bootstrap result, part of :class:`NoiseFloorResult.by_model`."""

    threshold: float
    cost: NDArray[np.float64]


@dataclass
class NoiseFloorResult:
    """Result of :func:`noise_floor`.

    Attributes:
        threshold: Scalar cutoff on the :func:`gsua_csb.costf`/``cost`` scale, for the chosen
            ``noise`` model.
        cost: (B,) bootstrap cost distribution, for the chosen ``noise`` model.
        accepted: (N,) bool, ``cost < threshold`` -- replaces ``cost < 1.5*cost.min()``.
        phi: Pooled Pearson dispersion estimate.
        k: Pooled method-of-moments global-NB dispersion estimate (``inf`` if degenerate).
        by_model: ``{"poisson": ..., "quasipoisson": ..., "nbglobal": ...}``, each a
            :class:`NoiseFloorByModel`, for the sensitivity comparison the module docstring
            describes -- quasi-Poisson and global-NB agreeing closely is reassuring; a large
            disagreement is itself diagnostic.
        model: Which model ``threshold``/``cost`` came from.
        quasi_fallback: True if quasi-Poisson degenerated to Poisson draws (``phi <= 1``).
        example_synthetic: ``{"poisson": ..., "quasipoisson": ..., "nbglobal": ...}``, each one
            example ``(n_outputs, len(xdata))`` synthetic dataset (the first bootstrap draw) --
            for plotting a representative noise-floor sample against real data, or for verifying
            e.g. that a cumulative row's NaN mask matches ``ydata``'s own.
    """

    threshold: float
    cost: NDArray[np.float64]
    accepted: NDArray[np.bool_]
    phi: float
    k: float
    by_model: dict[str, NoiseFloorByModel]
    model: NoiseModel
    quasi_fallback: bool
    example_synthetic: dict[str, NDArray[np.float64]]


def noise_floor(
    model: Model,
    xdata: ArrayLike,
    ydata: ArrayLike,
    x: ArrayLike,
    cost: ArrayLike,
    *,
    margin: float,
    alpha: float,
    b: int = 2000,
    noise: NoiseModel = "quasipoisson",
    quantile: float = 0.95,
    cumulative: ArrayLike | None = None,
    seed: int | np.random.Generator | None = None,
) -> NoiseFloorResult:
    """Noise-calibrated fit-acceptance threshold via parametric bootstrap.

    Replaces the scale-dependent ``lims = sum(cost < 1.5*cost.min())`` idiom. See the module
    docstring for the method and its positioning/caveats.

    Args:
        model: The model ``x``/``cost`` were estimated for (e.g. from
            :func:`gsua_csb.parameter_estimation`). ``model.fixed`` supplies the free-parameter
            count used in the pooled dispersion estimate's degrees of freedom.
        xdata: Domain the pool was fit against.
        ydata: (n_outputs, len(xdata)) reference data the pool was fit against, same orientation
            as :func:`gsua_csb.costf`/:func:`gsua_csb.parameter_estimation`. May contain NaN for
            missing samples.
        x: (N, Np) pool of estimated parameter vectors (e.g. ``PEResult.x``).
        cost: (N,) fit cost per pool member (e.g. ``PEResult.cost``), on the same scale
            :func:`gsua_csb.rcostf` produces.
        margin: The **raw** margin the pool was estimated with (e.g. ``PEResult.margin`` -- the
            same convention :func:`gsua_csb.parameter_estimation`'s own ``margin`` argument uses,
            not pre-offset). Required, not auto-recovered -- pass it through explicitly from
            whatever produced ``x``/``cost`` (typically a ``PEResult``).
        alpha: The cost exponent the pool was estimated with (e.g. ``PEResult.alpha``). Required
            for the same reason as ``margin``.
        b: Bootstrap sample size. Default 2000.
        noise: Noise model used for the primary ``threshold``/``cost``: ``"poisson"``,
            ``"quasipoisson"`` (default, recommended -- see below), or ``"nbglobal"``.
        quantile: Upper quantile of the bootstrap cost distribution defining the threshold.
            Default 0.95.
        cumulative: (n_outputs,) bool, one flag per fitted output row (default: all False). For a
            row whose ``ydata`` is a CUMULATIVE series, noise is estimated and generated on the
            INCIDENT (first-differenced) series and re-accumulated -- applying independent noise
            directly to a cumulative series would imply independent errors on a running sum,
            which is wrong.
        seed: Seed (or `numpy.random.Generator`) for the bootstrap draws, for reproducibility.

    Noise models (do not default to Poisson -- it is almost always far too tight for real
    surveillance-style data, since real observation processes are typically overdispersed
    relative to a Poisson count model):
        ``"poisson"``: ``var = mu``. Included for completeness/comparison, never the default.
        ``"quasipoisson"``: ``var = phi*mu``, ``phi`` = pooled Pearson dispersion
            ``sum(r**2/mu)/(n-p)``. RECOMMENDED DEFAULT. Implemented as NB2 with a per-point
            ``k = mu/(phi-1)`` (yields ``var = mu + mu**2/k = phi*mu`` exactly). Falls back to
            Poisson draws with a warning if ``phi <= 1`` (data equi/under-dispersed relative to
            Poisson, for which this reparameterization is degenerate).
        ``"nbglobal"``: ``var = mu + mu**2/k``, a single pooled ``k`` fitted by method of moments
            (``k = sum(mu**2)/(sum(r**2)-sum(mu))``). Report as a bound only -- a single global
            ``k`` is fitted almost entirely by the largest-mu points (e.g. epidemic peaks) and can
            predict implausible spread elsewhere. Falls back to Poisson draws with a warning if
            the method-of-moments estimate is degenerate (``k<=0``).
        ``by_model`` reports the threshold/distribution under all three regardless of ``noise``, as
        a sanity check.

    Returns:
        A :class:`NoiseFloorResult`. If the best pool member's cost exceeds the threshold, a
        warning fires: the model itself is rejected as a description of the data at the chosen
        quantile, which is diagnostically important and must not pass silently.

    Idiom replacing ``lims = sum(cost < 1.5*cost.min()); accepted_x = x[:lims]``::

        out = noise_floor(model, xdata, ydata, pe.x, pe.cost, margin=pe.margin, alpha=pe.alpha)
        ia = identifiability_analysis(model, pe.x[out.accepted], cost=pe.cost[out.accepted])

    Raises:
        ValueError: If ``margin`` (after the same raw-to-internal transform
            :func:`gsua_csb.rcostf` applies) is degenerate (``== 1``), or ``cumulative`` has the
            wrong length.
    """
    xdata = np.asarray(xdata, dtype=np.float64)
    ydata = np.atleast_2d(np.asarray(ydata, dtype=np.float64))
    x = np.atleast_2d(np.asarray(x, dtype=np.float64))
    cost = np.atleast_1d(np.asarray(cost, dtype=np.float64))

    n_outputs, length = ydata.shape

    if cumulative is None:
        cumulative_mask = np.zeros(n_outputs, dtype=bool)
    else:
        cumulative_mask = np.asarray(cumulative, dtype=bool)
        if cumulative_mask.shape[0] != n_outputs:
            raise ValueError(
                f"cumulative must have one entry per fitted output row ({n_outputs}), "
                f"got {cumulative_mask.shape[0]}"
            )

    # Same raw-to-internal transform rcostf applies internally, so the regulator formula below is
    # identical to what parameter_estimation actually scored with.
    margin_internal = abs(margin)
    if margin_internal < 1:
        margin_internal += 1
    if margin_internal == 1:
        raise ValueError(
            f"margin (internal value {margin_internal:.6g}) must not be 0 (or 1, raw) -- it "
            "degenerates the regulator to zero. The noise floor is only meaningful for a pool "
            "whose cost came from the rcostf-scored path."
        )

    rng = np.random.default_rng(seed)

    # theta*: best pool member by cost, not assumed to be the first row -- robust to caller sort
    # order.
    ibest = int(np.argmin(cost))
    best_cost = float(cost[ibest])
    theta_star = x[ibest]
    f = np.atleast_2d(np.asarray(model.evaluate(theta_star, xdata), dtype=np.float64))

    p_free = int((~model.fixed).sum())

    # Per-output incident/cumulative transform, built once and reused for both dispersion
    # estimation and synthetic generation. For a cumulative row, diffing lets NaN
    # propagate/smear naturally across any real gap in ydata -- mathematically correct (a
    # cumulative gap genuinely does not tell you the incident breakdown across it), and only ever
    # feeds a pooled scalar (phi/k), not a positional array, so no gap-bridging logic is needed.
    mu_grid = np.zeros((n_outputs, length))
    r_grid = np.full((n_outputs, length), np.nan)
    for i in range(n_outputs):
        if cumulative_mask[i]:
            f_inc = np.diff(np.concatenate([[0.0], f[i]]))
            y_inc = np.diff(np.concatenate([[0.0], ydata[i]]))
            mu_grid[i] = f_inc
            r_grid[i] = y_inc - f_inc
        else:
            mu_grid[i] = f[i]
            r_grid[i] = ydata[i] - f[i]
    mu_clamped = np.maximum(mu_grid, 0.0)

    # Pooled dispersion across every fitted output (not per-row): a single phi/k, matching how
    # this calibration is reported and used downstream.
    valid = np.isfinite(r_grid) & np.isfinite(mu_grid) & (mu_grid > 0)
    r_valid = r_grid[valid]
    mu_valid = mu_grid[valid]
    n_valid = r_valid.size
    if n_valid - p_free <= 0:
        warnings.warn(
            f"Only {n_valid} valid (non-NaN, mu>0) pooled residuals for {p_free} free "
            "parameters; dispersion estimates may be unreliable.",
            stacklevel=2,
        )
    denom = max(n_valid - p_free, 1)
    phi = float(np.sum(r_valid**2 / mu_valid) / denom)

    k_denom = float(np.sum(r_valid**2) - np.sum(mu_valid))
    if k_denom <= 0:
        k_global = float("inf")
        warnings.warn(
            "Method-of-moments global-NB k is degenerate (data at or below Poisson dispersion); "
            "nbglobal falls back to Poisson.",
            stacklevel=2,
        )
    else:
        k_global = float(np.sum(mu_valid**2) / k_denom)

    quasi_fallback = phi <= 1
    if quasi_fallback:
        warnings.warn(
            f"Pooled Pearson dispersion phi={phi:.6g} <= 1; quasipoisson falls back to Poisson draws.",
            stacklevel=2,
        )

    by_model: dict[str, NoiseFloorByModel] = {}
    example_synthetic: dict[str, NDArray[np.float64]] = {}
    for model_name in _MODELS:
        yb = _generate_bootstrap(
            rng, model_name, mu_clamped, cumulative_mask, ydata, phi, k_global, quasi_fallback, b
        )
        example_synthetic[model_name] = yb[:, :, 0].copy()
        cost_dist = np.zeros(b)
        for j in range(b):
            yb_slice = yb[:, :, j]
            n_valid_row = (~np.isnan(yb_slice)).sum(axis=1)
            regulator_j = np.nansum((yb_slice - yb_slice * margin_internal) ** 2, axis=1) / n_valid_row
            cost_dist[j] = costf(yb_slice, f, regulator_j, alpha)
        threshold_j = float(np.quantile(cost_dist, quantile))
        by_model[model_name] = NoiseFloorByModel(threshold=threshold_j, cost=cost_dist)

    chosen = by_model[noise]
    result = NoiseFloorResult(
        threshold=chosen.threshold,
        cost=chosen.cost,
        accepted=cost < chosen.threshold,
        phi=phi,
        k=k_global,
        by_model=by_model,
        model=noise,
        quasi_fallback=quasi_fallback,
        example_synthetic=example_synthetic,
    )

    if best_cost > result.threshold:
        warnings.warn(
            f"Best pool cost ({best_cost:.6g}) exceeds the noise-floor threshold "
            f"({result.threshold:.6g}) at the {quantile*100:.0f}% level: the model is rejected "
            "as a description of the data.",
            stacklevel=2,
        )

    return result


def _generate_bootstrap(
    rng: np.random.Generator,
    model_name: str,
    mu_clamped: NDArray[np.float64],
    cumulative_mask: NDArray[np.bool_],
    ydata: NDArray[np.float64],
    phi: float,
    k_global: float,
    quasi_fallback: bool,
    b: int,
) -> NDArray[np.float64]:
    """Draw ``b`` synthetic datasets (n_outputs, len(xdata), b) under the given noise model."""
    n_outputs, length = mu_clamped.shape
    mu_rep = np.broadcast_to(mu_clamped[:, :, None], (n_outputs, length, b))

    if model_name == "poisson":
        s = rng.poisson(mu_rep).astype(np.float64)
    elif model_name == "quasipoisson":
        if quasi_fallback:
            s = rng.poisson(mu_rep).astype(np.float64)
        else:
            k_point = np.maximum(mu_clamped / (phi - 1), np.finfo(np.float64).tiny)
            k_rep = np.broadcast_to(k_point[:, :, None], (n_outputs, length, b))
            p_rep = k_rep / (k_rep + mu_rep)
            s = rng.negative_binomial(k_rep, p_rep).astype(np.float64)
    elif model_name == "nbglobal":
        if np.isinf(k_global):
            s = rng.poisson(mu_rep).astype(np.float64)
        else:
            k_rep = np.full((n_outputs, length, b), max(k_global, np.finfo(np.float64).tiny))
            p_rep = k_rep / (k_rep + mu_rep)
            s = rng.negative_binomial(k_rep, p_rep).astype(np.float64)
    else:  # pragma: no cover - guarded by the Literal type / caller validation
        raise ValueError(f"Unknown noise model: {model_name!r}")

    s = np.where(np.isnan(s), 0.0, s)

    yb = np.zeros((n_outputs, length, b))
    for i in range(n_outputs):
        if cumulative_mask[i]:
            yb[i] = np.cumsum(s[i], axis=0)
        else:
            yb[i] = s[i]
        yb[i, np.isnan(ydata[i]), :] = np.nan
    return yb
