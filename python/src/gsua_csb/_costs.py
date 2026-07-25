"""Correlation-penalized, margin-normalized cost functions.

Python port of ``gsua_costf``/``gsua_rcostf``/``gsua_costfMulti``/``gsua_likecost``/
``gsua_covmetric``. All share the same shape: a shape-normalized MSE penalized by how well the
*shape* (correlation) of the candidate matches the reference, so a flat/rescaled fit still incurs
cost. A value below 1 means "within tolerance"; the toolbox's optimizers and stopping checks all
rely on that convention.

Two deliberate deviations from the MATLAB originals:

* MATLAB's ``inputs``/``len`` arguments (the shape of ``ydata``) are dropped here -- they are
  inferred from the array shape instead of passed redundantly by the caller.
* MATLAB's ``gsua_costf`` is missing the ``isnan(corr) -> 1`` guard that ``gsua_rcostf`` and
  ``gsua_costfMulti`` already have (a near-constant reference/candidate makes the correlation
  undefined). This is a fresh implementation, not a literal translation, so the guard is applied
  uniformly rather than reproducing that inconsistency.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def _pearson_rows(a: NDArray[np.float64], b: NDArray[np.float64]) -> NDArray[np.float64]:
    """Row-wise Pearson correlation between two (inputs, len) arrays, NaN where undefined."""
    a_c = a - a.mean(axis=1, keepdims=True)
    b_c = b - b.mean(axis=1, keepdims=True)
    num = (a_c * b_c).sum(axis=1)
    den = np.sqrt((a_c**2).sum(axis=1) * (b_c**2).sum(axis=1))
    with np.errstate(invalid="ignore", divide="ignore"):
        r = num / den
    return r


def costf(ydata: ArrayLike, yfunction: ArrayLike, regulator: ArrayLike, alpha: float = 2.0) -> float:
    """Scalar correlation-penalized MSE cost, externally-supplied regulator.

    MATLAB equivalent: ``gsua_costf``. ``regulator`` is typically built once by the caller (e.g. a
    margin-perturbed baseline) and reused across many calls inside an optimizer loop -- kept as an
    explicit parameter here for exactly that reason, unlike ``rcostf`` which builds its own.

    Args:
        ydata: (n_outputs, n_samples) or (n_samples,) reference data.
        yfunction: Same shape as ``ydata`` -- candidate/model output.
        regulator: (n_outputs,) or scalar normalization, one per output row.
        alpha: Exponent applied to each per-output cost term. Default 2 (matches ``gsua_pe``'s
            default).

    Returns:
        Scalar cost, mean of the per-output penalized costs. Below 1 means within tolerance.
    """
    ydata = np.atleast_2d(np.asarray(ydata, dtype=np.float64))
    yfunction = np.atleast_2d(np.asarray(yfunction, dtype=np.float64))
    regulator = np.atleast_1d(np.asarray(regulator, dtype=np.float64))
    n_outputs, length = ydata.shape

    cost = np.sum((ydata - yfunction) ** 2, axis=1) / length / regulator
    r = _pearson_rows(ydata, yfunction)
    r = np.where(np.isnan(r), 1.0, r)
    cost = ((2 - r) * cost) ** alpha
    return float(cost.sum() / n_outputs)


def rcostf(ydata: ArrayLike, yfunction: ArrayLike, margin: float = 1.1, alpha: float = 2.0) -> float:
    """Scalar correlation-penalized MSE cost with a self-derived regulator.

    MATLAB equivalent: ``gsua_rcostf``. Standalone alternative to :func:`costf` for use directly as
    a scalar objective with an optimizer -- the normalization is derived internally from ``ydata``
    (the MSE between ``ydata`` and a margin-perturbed copy of itself), so no external regulator is
    needed.

    Args:
        ydata: (n_outputs, n_samples) or (n_samples,) reference data.
        yfunction: Same shape as ``ydata`` -- candidate/model output.
        margin: Relative perturbation used to build the internal regulator. Values ``< 1`` are
            treated as a fraction and converted to ``1 + margin``. Default 1.1 (= 10%).
        alpha: Exponent applied to each per-output cost term. Default 2.

    Returns:
        Scalar cost. Below 1 means within tolerance.
    """
    ydata = np.atleast_2d(np.asarray(ydata, dtype=np.float64))
    yfunction = np.atleast_2d(np.asarray(yfunction, dtype=np.float64))
    margin = abs(margin)
    if margin < 1:
        margin += 1
    n_outputs, length = ydata.shape

    regulator = np.sum((ydata - ydata * margin) ** 2, axis=1) / length
    cost = np.sum((ydata - yfunction) ** 2, axis=1) / length / regulator
    r = _pearson_rows(ydata, yfunction)
    r = np.where(np.isnan(r), 1.0, r)
    cost = ((2 - r) * cost) ** alpha
    return float(cost.sum() / n_outputs)


def costf_multi(
    ydata: ArrayLike,
    yfunction: ArrayLike,
    margin: float = 0.1,
    alpha: float = 1.0,
) -> NDArray[np.float64]:
    """Vectorized correlation-penalized MSE cost, many runs against one reference.

    MATLAB equivalent: ``gsua_costfMulti``. Scores every run in a batch of model outputs against a
    single reference in one call -- used by sensitivity analysis (multi-objective branch) and CSB
    range refinement to reduce a whole population of simulations to one scalar cost per run.

    Fully vectorized over runs via numpy broadcasting (looping only over the small ``n_outputs``
    dimension), replacing MATLAB's ``parfor``-over-runs -- with typically hundreds to thousands of
    runs and a handful of outputs, looping the small axis and vectorizing the large one is the
    faster arrangement regardless of parallelism.

    Args:
        ydata: (n_runs, n_samples, n_outputs) array, one page of output per run (e.g. from
            :meth:`gsua_csb.Model.evaluate_batch`).
        yfunction: (n_outputs, n_samples) reference/nominal output.
        margin: Relative perturbation for the internal regulator. Values ``<= 1`` are treated as a
            fraction and converted to ``1 + margin``. Default 0.1.
        alpha: Exponent applied to each per-output cost term. Default 1.

    Returns:
        (n_runs,) array of scalar costs, one per run.
    """
    ydata = np.asarray(ydata, dtype=np.float64)
    if ydata.ndim == 2:
        ydata = ydata[:, :, None]
    yfunction = np.atleast_2d(np.asarray(yfunction, dtype=np.float64))
    n_runs, length, n_outputs = ydata.shape
    if margin <= 1:
        margin = 1 + abs(margin)

    regulator = np.sum((yfunction - yfunction * margin) ** 2, axis=1) / length  # (n_outputs,)

    cost = np.zeros((n_runs, n_outputs))
    for i in range(n_outputs):
        y = ydata[:, :, i]  # (n_runs, length)
        ref = yfunction[i]  # (length,)
        y_c = y - y.mean(axis=1, keepdims=True)
        ref_c = ref - ref.mean()
        num = y_c @ ref_c
        den = np.sqrt((y_c**2).sum(axis=1) * (ref_c**2).sum())
        with np.errstate(invalid="ignore", divide="ignore"):
            r = num / den
        r = np.where(np.isnan(r), 1.0, r)
        mse = np.sum((y - ref) ** 2, axis=1) / length / regulator[i]
        cost[:, i] = ((2 - r) * mse) ** alpha

    return cost.sum(axis=1) / n_outputs


def likecost(ydata: ArrayLike, ynom: ArrayLike, margin: float) -> NDArray[np.float64]:
    """Gaussian negative-log-likelihood cost relative to a margin baseline.

    MATLAB equivalent: ``gsua_likecost``. Computes the Gaussian negative log-likelihood of each run
    in ``ydata`` against ``ynom`` (standard deviation proportional to ``ynom * margin``), minus the
    negative log-likelihood of the margin-perturbed baseline itself -- zero when a run matches
    ``ynom`` as closely as the margin baseline does, growing as it deviates further.

    Args:
        ydata: (n_runs, n_samples) array, one row per run.
        ynom: (n_samples,) reference/nominal output.
        margin: Relative standard deviation for the likelihood weighting (absolute value used).

    Returns:
        (n_runs,) array of relative negative log-likelihood costs.
    """
    ydata = np.atleast_2d(np.asarray(ydata, dtype=np.float64))
    ynom = np.asarray(ynom, dtype=np.float64)
    margin = 1 + abs(margin)
    ymargin = ynom * margin
    desv = np.abs(ymargin**2)

    costmargin = np.nansum(np.log(2 * np.pi * desv) + (ymargin - ynom) ** 2 / desv) / 2
    cost = np.nansum(np.log(2 * np.pi * desv) + (ydata - ynom) ** 2 / desv, axis=1) / 2
    return cost - costmargin


def coverage_metric(
    y: ArrayLike,
    ydata: ArrayLike,
    margin: float = 0.1,
    alpha: float = 2.0,
    plevels: tuple[float, float] = (5, 95),
) -> tuple[float, float, NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Percentile-band coverage/tightness metrics for uncertainty-analysis results.

    MATLAB equivalent: ``gsua_covmetric``. Reduces a Monte-Carlo uncertainty-analysis output to a
    low/median/high percentile band per output and sample point, then scores that band with two
    costs on the same margin-normalized scale as :func:`costf`: below 1 means within tolerance.

    ``cost_data`` (accuracy) is the distance from the band's median curve to ``ydata`` -- high means
    the identified region's central tendency does not track the observed data. ``cost_band``
    (precision) is the distance between the low and high percentile curves themselves, independent
    of ``ydata`` -- low means repeated estimations/samples converged to a tight region of parameter
    space, *regardless* of how well that region fits the data. These can genuinely disagree: a model
    can converge to a thin, well-identified region that still does not perfectly track noisy or
    structurally-imperfect data. That is a model-fit limitation, not an identifiability problem --
    see the semi-automation identification routine's Phase 4/5 for how to act on the distinction.

    Args:
        y: (n_runs, n_samples) or (n_runs, n_samples, n_outputs) Monte-Carlo output (e.g. from
            :meth:`gsua_csb.Model.evaluate_batch`). Runs containing any non-finite value are
            excluded from the percentile band.
        ydata: (n_samples,) or (n_outputs, n_samples) reference data.
        margin: Relative tolerance for the shared normalization regulator. Values ``< 1`` are
            treated as a fraction. Default 0.1.
        alpha: Exponent applied to each per-output cost term. Default 2.
        plevels: (low, high) percentile levels defining the band. Default (5, 95).

    Returns:
        ``(cost_data, cost_band, p_low, p_median, p_high)`` -- the two scalar costs, and the
        (n_outputs, n_samples) low/median/high percentile curves actually used.
    """
    y = np.asarray(y, dtype=np.float64)
    if y.ndim == 2:
        y = y[:, :, None]
    ydata = np.atleast_2d(np.asarray(ydata, dtype=np.float64))

    valid = np.isfinite(y).all(axis=(1, 2))
    n_bad = int((~valid).sum())
    if n_bad:
        import warnings

        warnings.warn(
            f"{n_bad} of {y.shape[0]} runs contain non-finite values and were excluded "
            "from the percentile band",
            stacklevel=2,
        )
    y = y[valid]
    if y.shape[0] < 2:
        raise ValueError("Not enough finite runs left to compute a percentile band")

    p_low = np.percentile(y, plevels[0], axis=0).T   # (n_outputs, n_samples)
    p_med = np.percentile(y, 50, axis=0).T
    p_high = np.percentile(y, plevels[1], axis=0).T

    cost_data = rcostf(ydata, p_med, margin=margin, alpha=alpha)
    cost_band = rcostf(p_low, p_high, margin=margin, alpha=alpha)
    return cost_data, cost_band, p_low, p_med, p_high
