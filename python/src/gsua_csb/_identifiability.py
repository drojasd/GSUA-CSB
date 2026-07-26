"""Practical identifiability analysis: the Python replacement for ``gsua_ia``/``gsua_dia``.

Ports the shared data/statistics core of MATLAB's ``gsua_ia`` (plotting) and ``gsua_dia``
(headless) into a single headless function -- this package is headless-first throughout (plotting
is a separate, later module), so there is no need for two parallel variants here the way MATLAB
has one for "with plots" and one for "without."

One deliberate deviation: MATLAB's two functions actually compute the confidence interval on the
median differently -- ``gsua_ia`` uses the distribution-free order-statistic construction
(``gsua_medianCI``, ported here as :func:`gsua_csb.median_ci`), while ``gsua_dia`` uses a cruder
parametric normal approximation (``med +/- 1.96*sqrt(pi/2)*std/sqrt(n)``) that assumes the
estimates are asymptotically normal around the median. This module always uses the distribution-
free method (``gsua_ia``'s, the more defensible one for a set of repeated point-estimates whose
actual distribution is unknown) rather than reproducing both.

Also not ported: MATLAB's outlier removal calls a bundled third-party utility
(``Add Funcs/AntonSemechko-Multivariate-Outliers-.../DetectMultVarOutliers.m``), not part of the
core toolbox. This module implements the same idea (classic squared Mahalanobis distance vs. a
chi-squared critical value) directly instead of vendoring a third-party MATLAB file's exact
percentile-table convention.

One capability neither ``gsua_ia`` nor ``gsua_dia`` has: filtering repeated estimates by fit
quality before computing statistics. A multistart run that converged to a bad local optimum still
contributes its (essentially arbitrary) parameter values to the correlation matrix, confidence
interval, and clustering -- which can make a genuinely well-identified parameter look poorly
identified purely because an optimizer run failed, not because of real non-identifiability.
Passing ``cost=`` (e.g. ``PEResult.cost``) screens those runs out *before* anything else is
computed, via one of two methods (``cost_method``):

- ``"rtol"`` (default) -- keep runs within a fixed relative tolerance (``cost_rtol``) of the best
  cost. Simple and predictable, but a fixed percentage is an arbitrary answer to what is genuinely
  a data-dependent question.
- ``"gap"`` -- keep runs up to the largest relative jump in the *sorted* costs. This is the
  automatic version of how a "waterfall plot" (e.g. pyPESTO's) is read by eye: a cluster of
  similar, small costs (converged runs), a jump, then a scattered tail (failed runs). No tolerance
  to pick, but it can keep everything if costs form a smooth continuum with no clear jump.

Regardless of method, ``min_keep`` puts a floor under how aggressive the filter can be -- with a
small or noisy multistart sample, a strict cutoff can otherwise leave a single point (or none),
which silently produces degenerate or ``NaN`` correlation/clustering/confidence-interval results
downstream instead of a clear signal that something needs attention.

Multiple global minima: repeated estimations from a non-convex problem can land in distinct basins
of attraction instead of scattering around one point. When ``cluster=True``, this runs
``sklearn.cluster.SpectralClustering`` on the (fit-quality-filtered) normalized estimates for
candidate counts ``2..max_k`` and keeps the split with the best mean silhouette score
(``sklearn.metrics.silhouette_score``), the same model-selection idea as MATLAB's
``evalclusters``. If the best mean silhouette is below ``sil_threshold``, the estimates are
treated as a single basin -- a low silhouette means the best split found is not actually well
separated, so the sample is more likely unimodal than confidently multimodal.

Clustering runs *before* multivariate outlier removal, on purpose. Global Mahalanobis-distance
outlier detection assumes one unimodal population; run it before clustering and a real second
basin looks exactly like "outliers" relative to the pooled mean/covariance across both basins --
so outlier removal would delete the very thing ``cluster=True`` exists to find. When a genuine
multi-cluster split is found, every downstream summary statistic (``range``, ``correlation``,
``index``, ``nominal``) is instead computed from the *dominant* (largest) cluster's points only --
pooling a median or a correlation across separated basins isn't a meaningful summary either -- and
``outlier=True``, if also requested, is scoped to that single cluster, where the unimodal
assumption behind Mahalanobis distance is actually valid.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.stats import chi2
from sklearn.cluster import SpectralClustering
from sklearn.metrics import silhouette_score

from ._model import Model
from ._stats import median_ci

CostMethod = Literal["rtol", "gap"]


def _remove_multivariate_outliers(
    X: NDArray[np.float64], alpha: float
) -> tuple[NDArray[np.float64], NDArray[np.bool_]]:
    """Drop rows whose squared Mahalanobis distance exceeds the chi-squared ``1 - alpha`` critical value."""
    mean = X.mean(axis=0)
    cov = np.cov(X, rowvar=False)
    diff = X - mean
    inv_cov = np.linalg.pinv(np.atleast_2d(cov))
    md2 = np.einsum("ij,jk,ik->i", diff, inv_cov, diff)
    crit = chi2.ppf(1 - alpha, df=X.shape[1])
    keep = md2 < crit
    return X[keep], keep


def _gap_cutoff(cost: NDArray[np.float64], gap_ratio: float = 3.0) -> NDArray[np.bool_]:
    """Keep costs up to the largest relative jump in the sorted (log) costs, if any jump is salient.

    The automatic version of reading a waterfall plot: find where the sorted costs jump from a
    cluster of small, similar values to a scattered tail, and cut there. A smooth continuum of
    costs (e.g. evenly spaced) always has *some* largest gap between consecutive log-costs, so the
    largest gap alone isn't evidence of a real break -- it must also be at least ``gap_ratio``
    times the *median* gap to count as one, otherwise every point is treated as part of one
    continuum and nothing is dropped.
    """
    n = cost.shape[0]
    if n < 3:
        return np.ones(n, dtype=bool)
    order = np.argsort(cost)
    sorted_cost = cost[order]
    log_cost = np.log(np.maximum(sorted_cost, np.finfo(np.float64).tiny))
    gaps = np.diff(log_cost)
    max_gap = np.max(gaps)
    median_gap = np.median(gaps)
    if max_gap <= 0 or median_gap <= 0 or max_gap < gap_ratio * median_gap:
        return np.ones(n, dtype=bool)
    cut = int(np.argmax(gaps))  # keep sorted_cost[: cut + 1]
    keep_sorted = np.zeros(n, dtype=bool)
    keep_sorted[: cut + 1] = True
    keep = np.zeros(n, dtype=bool)
    keep[order] = keep_sorted
    return keep


def _filter_by_cost(
    X: NDArray[np.float64],
    cost: NDArray[np.float64],
    method: CostMethod,
    rtol: float,
    atol: float,
    gap_ratio: float,
    min_keep: int,
) -> tuple[NDArray[np.float64], int]:
    if method == "rtol":
        threshold = np.min(cost) * (1 + rtol) + atol
        keep = cost <= threshold
    elif method == "gap":
        keep = _gap_cutoff(cost, gap_ratio=gap_ratio)
    else:
        raise ValueError(f"Unknown cost_method: {method!r}, expected 'rtol' or 'gap'")

    min_keep = min(min_keep, X.shape[0])
    if keep.sum() < min_keep:
        # Floor: always keep at least the min_keep best-cost runs, even if the chosen method's
        # cutoff would otherwise leave fewer -- a degenerate 0- or 1-point result silently breaks
        # every downstream statistic instead of signaling that the filter was too aggressive.
        rank = np.argsort(np.argsort(cost))  # rank[i] = cost[i]'s position among all costs, 0-based
        keep = rank < min_keep

    n_removed = int((~keep).sum())
    return X[keep], n_removed


@dataclass
class ClusterInfo:
    """Multiple-global-minima clustering result (part of :class:`IdentifiabilityResult`).

    Attributes:
        num_clusters: Number of candidate global minima detected. ``1`` when clustering was not
            requested, or no split beat ``sil_threshold``.
        labels: (N,) cluster label (0-indexed) for each estimation run.
        centers: (num_clusters, Np) candidate global minima, one row per cluster -- the per-cluster
            median of the original (unnormalized) estimates.
        cluster_sizes: (num_clusters,) number of runs in each cluster.
        silhouette: Mean silhouette score of the chosen split, or ``None`` if clustering was not
            requested or run.
        dominant_cluster: Index (into ``centers``/``cluster_sizes``) of the largest cluster -- the
            one :class:`IdentifiabilityResult`'s ``range``/``correlation``/``index``/``nominal``
            are computed from when ``num_clusters > 1``. Always ``0`` when ``num_clusters == 1``.
    """

    num_clusters: int
    labels: NDArray[np.intp]
    centers: NDArray[np.float64]
    cluster_sizes: NDArray[np.intp]
    silhouette: float | None
    dominant_cluster: int = 0


@dataclass
class IdentifiabilityResult:
    """Result of :func:`identifiability_analysis`.

    Attributes:
        names: Parameter names, length Np (``model.names``).
        range: (Np, 2) updated confidence interval for each parameter, replacing ``model.range``.
            Computed from the dominant cluster's points only when ``cluster.num_clusters > 1``
            (see :class:`ClusterInfo`) -- pooling separated basins into one interval isn't
            meaningful.
        nominal: (Np,) median of the (possibly filtered) estimates used for ``range``.
        index: (Np,) identifiability index per parameter, in ``[0, 1]`` -- ``0`` well identified,
            ``1`` poorly identified. Combines interval width, mean correlation with other
            parameters, and count of strong (``|r| > 0.5``) correlations. Always ``0`` for a fixed
            parameter (``model.fixed``) -- "fixed" means "not being estimated," not "poorly
            estimated."
        correlation: (Np, Np) correlation matrix among the estimates used for ``range``. Entries
            involving a fixed parameter are ``NaN`` (undefined: a fixed parameter has zero variance
            across repeated estimations), except the diagonal, which is always 1.
        n_bad_fit_removed: Number of estimation runs dropped by the fit-quality filter (0 if
            ``cost=None``).
        n_outliers_removed: Number of estimation runs dropped by outlier removal (0 if
            ``outlier=False``).
        cluster: Multiple-global-minima clustering result.
    """

    names: list[str]
    range: NDArray[np.float64]
    nominal: NDArray[np.float64]
    index: NDArray[np.float64]
    correlation: NDArray[np.float64]
    n_bad_fit_removed: int
    n_outliers_removed: int
    cluster: ClusterInfo


def _summary_stats(
    X: NDArray[np.float64],
    model: Model,
    correction: bool,
    ci_alpha: float,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Compute (range, nominal, index, correlation) from a single (presumed-unimodal) point set."""
    lb0 = model.range[:, 0]
    ub0 = model.range[:, 1]
    width0 = ub0 - lb0
    Np = X.shape[1]
    free_mask = width0 > 0

    # A fixed parameter has zero range width and (if it stays fixed across every repeated
    # estimation, as gsua_pe/parameter_estimation guarantee) zero variance -- correlation against
    # it is mathematically undefined (0/0), not just numerically unstable. Leave fixed parameters
    # at well-defined defaults (NaN correlation, 0 identifiability index -- "fixed" means "not
    # being estimated", not "poorly estimated").
    correlation = np.full((Np, Np), np.nan)
    np.fill_diagonal(correlation, 1.0)
    if free_mask.sum() > 1:
        free_idx = np.where(free_mask)[0]
        correlation[np.ix_(free_idx, free_idx)] = np.corrcoef(X[:, free_idx], rowvar=False)

    lb, ub = median_ci(X.T, alpha=ci_alpha)
    if correction:
        lb = np.where(lb < lb0, lb0, lb)
        ub = np.where(ub > ub0, ub0, ub)

    boxin = np.zeros(Np)
    boxin[free_mask] = (ub[free_mask] - lb[free_mask]) / width0[free_mask]
    corrin = np.zeros(Np)
    extrin = np.zeros(Np)
    with np.errstate(invalid="ignore"):
        corrin[free_mask] = np.nanmean(np.abs(correlation[:, free_mask]), axis=0)
        extrin[free_mask] = np.nansum(np.abs(correlation[:, free_mask]) > 0.5, axis=0) / Np
    index = (2 * boxin + corrin + extrin) / 4

    median_est = np.median(X, axis=0)
    return np.column_stack([lb, ub]), median_est, index, correlation


def identifiability_analysis(
    model: Model,
    estimates: ArrayLike,
    cost: ArrayLike | None = None,
    cost_method: CostMethod = "rtol",
    cost_rtol: float = 0.1,
    cost_atol: float = 1e-8,
    cost_gap_ratio: float = 3.0,
    min_keep: int = 3,
    correction: bool = False,
    outlier: bool = False,
    outlier_alpha: float = 0.025,
    cluster: bool = False,
    max_k: int = 5,
    sil_threshold: float = 0.5,
    ci_alpha: float = 0.05,
    seed: int | None = None,
) -> IdentifiabilityResult:
    """Analyze practical parameter identifiability from a set of repeated estimations.

    Args:
        model: Source of the original parameter range (``model.range``) that the new confidence
            interval is measured against/optionally clipped to.
        estimates: (N, Np) repeated parameter estimates, e.g. ``PEResult.x`` from
            :func:`gsua_csb.parameter_estimation` called with ``n > 1``.
        cost: (N,) fit cost for each row of ``estimates``, e.g. ``PEResult.cost``. If given, runs
            that converged to a much worse cost than the best one are dropped before any other
            statistic is computed -- see ``cost_method``. ``None`` (default) skips this filter.
        cost_method: How to decide the fit-quality cutoff. ``"rtol"`` (default) keeps runs within
            ``cost_rtol`` of the best cost -- simple and predictable. ``"gap"`` keeps runs up to
            the largest relative jump in the sorted costs -- data-driven, no tolerance to pick, but
            can keep everything if costs form a smooth continuum. See the module docstring.
        cost_rtol: Relative tolerance above the best cost for ``cost_method="rtol"``. Default 0.1
            (10%, matching this package's other margin-normalized-cost conventions). Ignored for
            ``cost_method="gap"``.
        cost_atol: Absolute tolerance added to the ``"rtol"`` threshold, so the filter doesn't
            become vacuously strict when the best cost is at or near zero (e.g. a synthetic,
            noiseless fit). Ignored for ``cost_method="gap"``.
        cost_gap_ratio: For ``cost_method="gap"``, how many times larger than the *median*
            consecutive gap (in sorted log-cost) the largest gap must be to count as a real break
            rather than ordinary spacing in a smooth continuum. Default 3.0. Ignored for
            ``cost_method="rtol"``.
        min_keep: Minimum number of runs the fit-quality filter is allowed to leave (falls back to
            keeping the ``min_keep`` best-cost runs if the chosen ``cost_method`` would leave
            fewer). Guards against a degenerate 0- or 1-point result silently propagating into
            ``NaN``/meaningless correlation, confidence-interval, and clustering results. Ignored
            if ``cost=None``.
        correction: If True, clip the new confidence interval to ``model.range`` where it would
            otherwise extend beyond it.
        outlier: If True, remove multivariate outliers (via Mahalanobis distance) before computing
            statistics -- scoped to the dominant cluster's points if ``cluster=True`` found a
            genuine multi-cluster split (see the module docstring for why outlier removal runs
            after, not before, clustering).
        outlier_alpha: Significance level for the outlier chi-squared cutoff. Default 0.025.
        cluster: If True, run spectral clustering to detect multiple candidate global minima.
        max_k: Largest cluster count to test.
        sil_threshold: Minimum mean silhouette score required to accept a multi-cluster split over
            a single basin.
        ci_alpha: Significance level for the median confidence interval. Default 0.05 (95% CI).
        seed: Seed for ``SpectralClustering``'s stochastic eigensolver, for reproducibility.

    Returns:
        An :class:`IdentifiabilityResult`.

    Raises:
        ValueError: If ``cost`` is given and its length doesn't match ``estimates``, or
            ``cost_method`` is not one of the two supported values.
    """
    X = np.atleast_2d(np.asarray(estimates, dtype=np.float64))

    n_bad_fit_removed = 0
    if cost is not None:
        cost = np.atleast_1d(np.asarray(cost, dtype=np.float64))
        if cost.shape[0] != X.shape[0]:
            raise ValueError(f"cost must have {X.shape[0]} entries (one per row of estimates), got {cost.shape[0]}")
        X, n_bad_fit_removed = _filter_by_cost(
            X, cost, cost_method, cost_rtol, cost_atol, cost_gap_ratio, min_keep
        )

    lb0 = model.range[:, 0]
    width0 = model.range[:, 1] - lb0
    free_mask = width0 > 0

    def normalize(pts: NDArray[np.float64]) -> NDArray[np.float64]:
        norm = np.zeros_like(pts)
        norm[:, free_mask] = (pts[:, free_mask] - lb0[free_mask]) / width0[free_mask]
        return norm

    all_median = np.median(X, axis=0)
    cluster_info = ClusterInfo(
        num_clusters=1,
        labels=np.zeros(X.shape[0], dtype=np.intp),
        centers=all_median[None, :],
        cluster_sizes=np.array([X.shape[0]]),
        silhouette=None,
        dominant_cluster=0,
    )

    # Clustering runs on the fit-quality-filtered points, BEFORE any outlier removal -- see the
    # module docstring for why: outlier removal assumes one unimodal population, and running it
    # first would treat a genuine second basin as outliers relative to the pooled distribution.
    stats_source = X
    if cluster:
        n_runs = X.shape[0]
        kmax = min(max_k, n_runs - 1)
        if kmax >= 2:
            normalized = normalize(X)
            best_k = 1
            best_sil = 0.0
            for k in range(2, kmax + 1):
                n_neighbors = min(10, n_runs - 1)
                labels = SpectralClustering(
                    n_clusters=k,
                    affinity="nearest_neighbors",
                    n_neighbors=n_neighbors,
                    random_state=seed,
                ).fit_predict(normalized)
                if len(np.unique(labels)) < 2:
                    continue
                sil = silhouette_score(normalized, labels)
                if sil > best_sil:
                    best_sil = sil
                    best_k = k

            if best_k > 1 and best_sil >= sil_threshold:
                labels = SpectralClustering(
                    n_clusters=best_k,
                    affinity="nearest_neighbors",
                    n_neighbors=min(10, n_runs - 1),
                    random_state=seed,
                ).fit_predict(normalized)
                centers = np.array([np.median(X[labels == k], axis=0) for k in range(best_k)])
                sizes = np.array([int(np.sum(labels == k)) for k in range(best_k)])
                dominant = int(np.argmax(sizes))
                cluster_info = ClusterInfo(
                    num_clusters=best_k,
                    labels=labels,
                    centers=centers,
                    cluster_sizes=sizes,
                    silhouette=best_sil,
                    dominant_cluster=dominant,
                )
                # Pooling separated basins into one median/correlation/interval isn't a meaningful
                # summary -- report the dominant (largest) basin's statistics instead. The other
                # basins are still fully available via cluster.centers.
                stats_source = X[labels == dominant]

    n_outliers_removed = 0
    if outlier:
        stats_source, keep = _remove_multivariate_outliers(stats_source, alpha=outlier_alpha)
        n_outliers_removed = int((~keep).sum())

    rng, nominal, index, correlation = _summary_stats(stats_source, model, correction, ci_alpha)

    return IdentifiabilityResult(
        names=list(model.names),
        range=rng,
        nominal=nominal,
        index=index,
        correlation=correlation,
        n_bad_fit_removed=n_bad_fit_removed,
        n_outliers_removed=n_outliers_removed,
        cluster=cluster_info,
    )
