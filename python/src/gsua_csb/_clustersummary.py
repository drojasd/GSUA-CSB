"""Per-cluster, per-parameter summary with an identifiability verdict: the Python ``gsua_clustersummary``.

For each parameter, in each cluster, reports where the group CENTER sits and how WIDE the plausible
range is, plus a suggested verdict -- the summary used in the Semillero epidemiology hands-on.

Two kinds of interval, answering different questions:
    * ``median_ci_*`` / ``mean_ci_*`` describe the CENTER of the group (how well pinned the central
      value is): the distribution-free CI of the median (:func:`gsua_csb.median_ci`) and the
      Student-t CI of the mean. They shrink as the cluster grows.
    * The 2.5-97.5 PERCENTILE interval ``[p2_5, p97_5]`` describes the SPREAD of plausible values.
      Its width relative to the parameter's ORIGINAL range (``pct_width_rel``) is what the verdict
      uses, and it does NOT shrink just because more estimates were collected.

Verdict on ``pct_width_rel``: ``< 0.25`` -> ``"identifiable"``; ``> 0.50`` -> ``"not identifiable"``;
in between -> ``"dubious"``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.stats import t as _student_t

from ._identifiability import IdentifiabilityResult
from ._stats import median_ci


@dataclass
class ClusterSummaryRow:
    """One (cluster, parameter) row of :func:`cluster_summary`. Tabulate with
    ``pandas.DataFrame([vars(r) for r in rows])`` if desired."""

    cluster: int
    parameter: str
    n: int
    median: float
    median_ci_low: float
    median_ci_high: float
    mean: float
    mean_ci_low: float
    mean_ci_high: float
    p2_5: float
    p97_5: float
    pct_width_rel: float
    verdict: str


def cluster_summary(result: IdentifiabilityResult) -> list[ClusterSummaryRow]:
    """Per-cluster, per-parameter statistical summary of an identifiability result.

    MATLAB equivalent: ``gsua_clustersummary``. See the module docstring for what each interval
    means and how the verdict is derived.

    Args:
        result: Result from :func:`gsua_csb.identifiability_analysis` (uses ``result.names`` and
            ``result.cluster``'s per-cluster ``estimates`` and ``orig_range``).

    Returns:
        A list of :class:`ClusterSummaryRow`, one per (cluster, parameter).
    """
    names = result.names
    ci = result.cluster
    orig = np.asarray(ci.orig_range, dtype=np.float64)
    range_w = orig[:, 1] - orig[:, 0]

    rows: list[ClusterSummaryRow] = []
    for c in range(ci.num_clusters):
        est = ci.estimates[c]            # (n_c, Np)
        n_c = est.shape[0]
        mlo, mhi = median_ci(est.T, alpha=0.05)
        for i, name in enumerate(names):
            e = est[:, i]
            med = float(np.median(e))
            mean = float(np.mean(e))
            sd = float(np.std(e, ddof=1)) if n_c > 1 else 0.0
            moe = float(_student_t.ppf(0.975, n_c - 1) * sd / np.sqrt(n_c)) if n_c > 1 else 0.0
            p2_5 = float(np.percentile(e, 2.5))
            p97_5 = float(np.percentile(e, 97.5))
            pw = float((p97_5 - p2_5) / range_w[i]) if range_w[i] > 0 else 0.0
            verdict = "identifiable" if pw < 0.25 else ("not identifiable" if pw > 0.5 else "dubious")
            rows.append(ClusterSummaryRow(
                cluster=c, parameter=name, n=n_c, median=med,
                median_ci_low=float(mlo[i]), median_ci_high=float(mhi[i]),
                mean=mean, mean_ci_low=mean - moe, mean_ci_high=mean + moe,
                p2_5=p2_5, p97_5=p97_5, pct_width_rel=pw, verdict=verdict,
            ))
    return rows
