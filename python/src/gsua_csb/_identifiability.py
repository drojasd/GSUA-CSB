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

Multiple global minima: repeated estimations from a non-convex problem can land in distinct basins
of attraction instead of scattering around one point. When ``cluster=True``, this runs
``sklearn.cluster.SpectralClustering`` on the normalized estimates for candidate counts
``2..max_k`` and keeps the split with the best mean silhouette score (``sklearn.metrics.
silhouette_score``), the same model-selection idea as MATLAB's ``evalclusters``. If the best mean
silhouette is below ``sil_threshold``, the estimates are treated as a single basin -- a low
silhouette means the best split found is not actually well separated, so the sample is more likely
unimodal than confidently multimodal.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.stats import chi2
from sklearn.cluster import SpectralClustering
from sklearn.metrics import silhouette_score

from ._model import Model
from ._stats import median_ci


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
    """

    num_clusters: int
    labels: NDArray[np.intp]
    centers: NDArray[np.float64]
    cluster_sizes: NDArray[np.intp]
    silhouette: float | None


@dataclass
class IdentifiabilityResult:
    """Result of :func:`identifiability_analysis`.

    Attributes:
        names: Parameter names, length Np (``model.names``).
        range: (Np, 2) updated confidence interval for each parameter, replacing ``model.range``.
        nominal: (Np,) median of the (possibly outlier-filtered) estimates.
        index: (Np,) identifiability index per parameter, in ``[0, 1]`` -- ``0`` well identified,
            ``1`` poorly identified. Combines interval width, mean correlation with other
            parameters, and count of strong (``|r| > 0.5``) correlations.
        correlation: (Np, Np) correlation matrix among the (possibly outlier-filtered) estimates.
        n_outliers_removed: Number of estimation runs dropped by outlier removal (0 if
            ``outlier=False``).
        cluster: Multiple-global-minima clustering result.
    """

    names: list[str]
    range: NDArray[np.float64]
    nominal: NDArray[np.float64]
    index: NDArray[np.float64]
    correlation: NDArray[np.float64]
    n_outliers_removed: int
    cluster: ClusterInfo


def identifiability_analysis(
    model: Model,
    estimates: ArrayLike,
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
        correction: If True, clip the new confidence interval to ``model.range`` where it would
            otherwise extend beyond it.
        outlier: If True, remove multivariate outliers from ``estimates`` (via Mahalanobis
            distance) before computing statistics.
        outlier_alpha: Significance level for the outlier chi-squared cutoff. Default 0.025.
        cluster: If True, run spectral clustering to detect multiple candidate global minima.
        max_k: Largest cluster count to test.
        sil_threshold: Minimum mean silhouette score required to accept a multi-cluster split over
            a single basin.
        ci_alpha: Significance level for the median confidence interval. Default 0.05 (95% CI).
        seed: Seed for ``SpectralClustering``'s stochastic eigensolver, for reproducibility.

    Returns:
        An :class:`IdentifiabilityResult`.
    """
    X = np.atleast_2d(np.asarray(estimates, dtype=np.float64))
    n_outliers_removed = 0
    if outlier:
        X, keep = _remove_multivariate_outliers(X, alpha=outlier_alpha)
        n_outliers_removed = int((~keep).sum())

    lb0 = model.range[:, 0]
    ub0 = model.range[:, 1]
    width0 = ub0 - lb0

    normalized = (X - lb0) / width0  # (N, Np)
    correlation = np.corrcoef(X, rowvar=False)
    correlation = np.atleast_2d(correlation)

    lb, ub = median_ci(X.T, alpha=ci_alpha)
    if correction:
        lb = np.where(lb < lb0, lb0, lb)
        ub = np.where(ub > ub0, ub0, ub)

    Np = X.shape[1]
    boxin = (ub - lb) / width0
    corrin = np.mean(np.abs(correlation), axis=0)
    extrin = np.sum(np.abs(correlation) > 0.5, axis=0) / Np
    index = (2 * boxin + corrin + extrin) / 4

    median_est = np.median(X, axis=0)

    cluster_info = ClusterInfo(
        num_clusters=1,
        labels=np.zeros(X.shape[0], dtype=np.intp),
        centers=median_est[None, :],
        cluster_sizes=np.array([X.shape[0]]),
        silhouette=None,
    )
    if cluster:
        n_runs = normalized.shape[0]
        kmax = min(max_k, n_runs - 1)
        if kmax >= 2:
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
                cluster_info = ClusterInfo(
                    num_clusters=best_k,
                    labels=labels,
                    centers=centers,
                    cluster_sizes=sizes,
                    silhouette=best_sil,
                )

    return IdentifiabilityResult(
        names=list(model.names),
        range=np.column_stack([lb, ub]),
        nominal=median_est,
        index=index,
        correlation=correlation,
        n_outliers_removed=n_outliers_removed,
        cluster=cluster_info,
    )
