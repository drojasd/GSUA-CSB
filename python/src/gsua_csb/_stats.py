"""Distribution-free statistics used by identifiability analysis.

Python port of ``gsua_medianCI`` and ``gsua_depth``.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.stats import binom, rankdata


def median_ci(x: ArrayLike, alpha: float = 0.05) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Distribution-free confidence interval for the median of each row of ``x``.

    MATLAB equivalent: ``gsua_medianCI``. Uses the order-statistics/binomial construction
    ``CI = [X_(r), X_(s)]``, where ``X_(k)`` is the k-th order statistic and r, s are chosen by
    inverting a Binomial(n, 1/2) distribution to achieve coverage approximately ``1 - alpha``. Makes
    no normality assumption.

    Args:
        x: (n_vars, n_samples) array. Each row is treated as an independent set of samples for one
            variable (e.g. one parameter's repeated estimates).
        alpha: Significance level. Default 0.05 for a 95% CI.

    Returns:
        ``(ci_lower, ci_upper)``, each (n_vars,).
    """
    x = np.atleast_2d(np.asarray(x, dtype=np.float64))
    n = x.shape[1]

    r = 1
    while r <= n and binom.cdf(r - 1, n, 0.5) < alpha / 2:
        r += 1
    r = min(r, n)

    s = n
    while s >= 1 and binom.cdf(s - 1, n, 0.5) > 1 - alpha / 2:
        s -= 1
    s = max(s, 1)

    sorted_x = np.sort(x, axis=1)
    return sorted_x[:, r - 1], sorted_x[:, s - 1]


def band_depth(x: ArrayLike) -> tuple[NDArray[np.float64], NDArray[np.intp]]:
    """Generalized band depth of a set of curves, ranked deepest-first.

    MATLAB equivalent: ``gsua_depth``. Computes the generalized band depth (Lopez-Pintado & Romo)
    of every curve relative to the full set -- a measure of how "central" each curve is among the
    others. Used to pick a representative reference output automatically when no experimental/
    nominal output is supplied for multi-objective sensitivity analysis.

    Fully vectorized (loops only over the small ``n_outputs`` axis via rank computation per sample
    point); MATLAB's ``parallel`` flag has no equivalent here since there is no per-curve Python
    loop to parallelize.

    Args:
        x: (n_curves, n_samples) or (n_curves, n_samples, n_outputs) array. Depth is averaged
            across outputs when 3-D.

    Returns:
        ``(depth, idx)`` -- (n_curves,) average band depth per curve, and the curve indices sorted
        from deepest (most central, ``idx[0]``) to shallowest.
    """
    x = np.asarray(x, dtype=np.float64)
    if x.ndim == 2:
        x = x[:, :, None]
    n, d, n_outputs = x.shape

    depth = np.zeros((n, n_outputs))
    denom = n * (n + 1) / 2
    for h in range(n_outputs):
        posi = np.apply_along_axis(lambda col: rankdata(col, method="min"), axis=0, arr=x[:, :, h])
        depth[:, h] = ((posi - 1) * (n - posi) + n).sum(axis=1) / denom

    depth = depth.mean(axis=1)
    idx = np.argsort(-depth)
    return depth, idx
