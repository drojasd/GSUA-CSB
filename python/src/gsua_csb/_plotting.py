"""Plotting: the Python replacement for ``gsua_plot``.

MATLAB's ``gsua_plot`` is a single function dispatched on a ``plot_type`` string with a variable,
positional argument list whose meaning changes per case (``gsua_plot('Bar', T, S)`` vs.
``gsua_plot('Bar', T, S, t, tref)`` are different plots entirely). This module instead exposes one
small function per plot type, each taking the relevant result dataclass directly
(:class:`gsua_csb.UncertaintyResult`, :class:`gsua_csb.SensitivityResult`,
:class:`gsua_csb.IdentifiabilityResult`, :class:`gsua_csb.MCFResult`) -- real keyword arguments
instead of positional-count-dependent dispatch, and every function accepts an optional ``ax`` (or
``axes``) so callers can compose these into their own figure layouts, following the standard
``matplotlib``/``pandas`` convention.

Not ported (deliberate simplifications, not oversights):

- ``'Pie'`` -- superseded by :func:`plot_sensitivity_bar`/:func:`plot_identifiability_index`. Pie
  charts encode magnitude via angle, which is harder to compare accurately than a bar's length;
  modern visualization guidance generally prefers bars for exactly this comparison.
- ``'FractionalSensitivityPlots'``/``'TotalSensitivityPlots'`` (one subplot per parameter) --
  superseded by :func:`plot_sensitivity_area`'s single stacked-area view, which scales to more
  parameters without a growing subplot grid and shows relative contribution directly.
- ``'ScatterOutput'``/``'ScatterParameter'`` -- exploratory pairwise/output scatter grids with no
  dedicated result dataclass behind them; a caller with a design matrix or output array can build
  these directly with ``matplotlib`` (``ax.scatter(M[:, i], Y[:, t])``) with no toolbox-specific
  logic to wrap.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import ArrayLike, NDArray

if TYPE_CHECKING:
    from matplotlib.axes import Axes

    from ._identifiability import IdentifiabilityResult
    from ._sensitivity import SensitivityResult
    from ._uncertainty import MCFResult, UncertaintyResult


def plot_uncertainty(
    ua: "UncertaintyResult",
    ax: "Axes | None" = None,
    color: str = "tab:blue",
    nominal_color: str = "tab:red",
    alpha: float = 0.3,
) -> "Axes":
    """Plot a Monte Carlo uncertainty-analysis ensemble with the reference output overlaid.

    MATLAB equivalent: ``gsua_plot('UncertaintyAnalysis', ...)``. Every row of ``ua.Y`` is drawn as
    a thin translucent line ("spaghetti plot"), with ``ua.y_nom`` overlaid in a contrasting color
    and full opacity.

    Args:
        ua: Result from :func:`gsua_csb.uncertainty_analysis`.
        ax: Axes to draw on. A new figure/axes is created if not given.
        color: Line color for the ensemble.
        nominal_color: Line color for the reference/nominal curve.
        alpha: Opacity for each ensemble line (low, so overlapping density is visible).

    Returns:
        The Axes drawn on.
    """
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots()

    if ua.xdata is None:
        ax.boxplot(ua.Y.ravel())
        ax.axhline(float(np.ravel(ua.y_nom)[0]), color=nominal_color, label="Nominal/reference")
    else:
        for row in ua.Y:
            ax.plot(ua.xdata, row, color=color, alpha=alpha, linewidth=0.7)
        ax.plot(ua.xdata, ua.y_nom, color=nominal_color, linewidth=2, label="Nominal/reference")
        ax.set_xlabel("Domain")
    ax.set_ylabel("Output")
    ax.set_title(f"Uncertainty analysis (N={ua.Y.shape[0]})")
    ax.legend()
    return ax


def plot_sensitivity_bar(sens: "SensitivityResult", index: str = "Si", ax: "Axes | None" = None) -> "Axes":
    """Horizontal bar chart of a sensitivity result's scalar indices, sorted ascending.

    MATLAB equivalent: ``gsua_plot('Bar', T, S)``.

    Args:
        sens: Result from :func:`gsua_csb.sensitivity_analysis`.
        index: Which scalar index to plot -- ``"Si"`` (first-order) or ``"STi"`` (total-order).
        ax: Axes to draw on. A new figure/axes is created if not given.

    Returns:
        The Axes drawn on.
    """
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots()

    values = getattr(sens, index)
    order = np.argsort(values)
    ax.barh(range(len(values)), values[order])
    ax.set_yticks(range(len(values)))
    ax.set_yticklabels([sens.names[i] for i in order])
    ax.set_xlabel(index)
    ax.set_title(f"Sensitivity indices ({sens.method})")
    return ax


def plot_sensitivity_area(
    sens: "SensitivityResult", xdata: ArrayLike, index: str = "Si_vec", ax: "Axes | None" = None
) -> "Axes":
    """Stacked-area plot of a sensitivity result's time-dependent indices.

    MATLAB equivalent: ``gsua_plot('FractionalSensitivityArea', ...)``/
    ``gsua_plot('TotalSensitivityArea', ...)``. Parameters are stacked largest-mean-contribution
    first (bottom of the stack), so the most influential parameter is easiest to read off.

    Args:
        sens: Result from :func:`gsua_csb.sensitivity_analysis`.
        xdata: Domain points matching the columns of ``Si_vec``/``STi_vec`` (e.g. ``model.domain``).
        index: Which time-dependent index to plot -- ``"Si_vec"`` or ``"STi_vec"``.
        ax: Axes to draw on. A new figure/axes is created if not given.

    Returns:
        The Axes drawn on.

    Raises:
        ValueError: If ``sens`` has no time-dependent indices (``method="oat"``).
    """
    import matplotlib.pyplot as plt

    S = getattr(sens, index)
    if S is None:
        raise ValueError(f"{index} is None for method={sens.method!r} (e.g. 'oat' has no time-dependent indices)")

    if ax is None:
        _, ax = plt.subplots()

    order = np.argsort(-np.nanmean(S, axis=1))
    ax.stackplot(np.asarray(xdata), S[order], labels=[sens.names[i] for i in order])
    ax.set_xlabel("Domain")
    ax.set_ylabel(index)
    ax.set_title(f"Time-dependent sensitivity indices ({sens.method})")
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0))
    return ax


def plot_identifiability_correlation(
    ia: "IdentifiabilityResult", ax: "Axes | None" = None, cmap: str = "coolwarm"
) -> "Axes":
    """Heatmap of the parameter correlation matrix from a practical identifiability analysis.

    Entries involving a fixed parameter (``NaN`` in ``ia.correlation``, see
    :func:`gsua_csb.identifiability_analysis`) render as blank/masked cells rather than a color.

    Args:
        ia: Result from :func:`gsua_csb.identifiability_analysis`.
        ax: Axes to draw on. A new figure/axes is created if not given.
        cmap: Colormap name, centered at zero. Default ``"coolwarm"``.

    Returns:
        The Axes drawn on.
    """
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots()

    corr = np.ma.masked_invalid(ia.correlation)
    im = ax.imshow(corr, vmin=-1, vmax=1, cmap=cmap)
    n = len(ia.names)
    ax.set_xticks(range(n))
    ax.set_xticklabels(ia.names, rotation=45, ha="right")
    ax.set_yticks(range(n))
    ax.set_yticklabels(ia.names)
    plt.colorbar(im, ax=ax, label="Correlation")
    ax.set_title("Parameter correlation")
    return ax


def plot_identifiability_graph(
    ia: "IdentifiabilityResult", ax: "Axes | None" = None, threshold: float = 0.5, cmap: str = "jet"
) -> "Axes":
    """Network view of strong pairwise correlations from a practical identifiability analysis.

    MATLAB equivalent: ``gsua_ia``'s "Identifiability graph" figure. An edge is drawn between two
    parameters when ``|correlation| > threshold``; node color ("heat") is the identifiability
    index. A cluster of parameters all connected to each other reads as "these are entangled --
    the data constrains some combination of them, not each individually" (the toolbox's namesake
    example: a two-parameter product `a*b` in a model where only the product is identifiable).

    If ``ia.correlation_reliable`` is ``False`` (too few points fed ``ia.correlation`` -- see
    :func:`gsua_csb.identifiability_analysis`), every entry is ``NaN`` and no edges are drawn: the
    graph renders as isolated nodes rather than the complete graph a naive reading of ``NaN``-free
    but degenerate 2-point correlations (always exactly +-1) would otherwise produce.

    Args:
        ia: Result from :func:`gsua_csb.identifiability_analysis`.
        ax: Axes to draw on. A new figure/axes is created if not given.
        threshold: Minimum ``|correlation|`` to draw an edge. Default 0.5, matching
            ``identifiability_analysis``'s own "strong correlation" convention (``extrin``).
        cmap: Colormap name for node color. Default ``"jet"``, matching the MATLAB original.

    Returns:
        The Axes drawn on.
    """
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    if ax is None:
        _, ax = plt.subplots()

    n = len(ia.names)
    angles = 2 * np.pi * np.arange(n) / max(n, 1)
    pos = np.column_stack([np.cos(angles), np.sin(angles)])

    with np.errstate(invalid="ignore"):
        strong = np.abs(ia.correlation) > threshold
    np.fill_diagonal(strong, False)
    edges = [(pos[i], pos[j]) for i, j in zip(*np.where(np.triu(strong, k=1)))]
    if edges:
        ax.add_collection(LineCollection(edges, colors="gray", linewidths=1.0, zorder=1))

    sc = ax.scatter(pos[:, 0], pos[:, 1], c=ia.index, cmap=cmap, s=250, zorder=2, edgecolors="black")
    for i, name in enumerate(ia.names):
        ax.annotate(name, pos[i], textcoords="offset points", xytext=(0, 10), ha="center", fontsize=9)

    plt.colorbar(sc, ax=ax, label="Identifiability Index")
    ax.set_xlim(-1.4, 1.4)
    ax.set_ylim(-1.4, 1.4)
    ax.set_aspect("equal")
    ax.axis("off")
    title = "Strong correlations" if ia.correlation_reliable else "Strong correlations (undefined -- too few points)"
    ax.set_title(title)
    return ax


def plot_identifiability_index(ia: "IdentifiabilityResult", ax: "Axes | None" = None) -> "Axes":
    """Horizontal bar chart of the per-parameter identifiability index, sorted ascending.

    Args:
        ia: Result from :func:`gsua_csb.identifiability_analysis`.
        ax: Axes to draw on. A new figure/axes is created if not given.

    Returns:
        The Axes drawn on.
    """
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots()

    order = np.argsort(ia.index)
    ax.barh(range(len(ia.index)), ia.index[order])
    ax.axvline(0.5, color="gray", linestyle="--", linewidth=1)
    ax.set_yticks(range(len(ia.index)))
    ax.set_yticklabels([ia.names[i] for i in order])
    ax.set_xlabel("Identifiability index (0 = well identified, 1 = poorly identified)")
    ax.set_title("Practical identifiability")
    return ax


def _ecdf_step(ax: "Axes", x: NDArray[np.float64], **kwargs) -> None:
    x = np.sort(np.asarray(x, dtype=np.float64))
    if x.size == 0:
        return
    y = np.arange(1, x.size + 1) / x.size
    ax.step(x, y, where="post", **kwargs)


def plot_mcf(
    mcf: "MCFResult", m_free: ArrayLike, axes: "NDArray[Axes] | None" = None
) -> "NDArray[Axes]":
    """Empirical CDF plots (prior vs. low vs. high) for a Monte Carlo filtering result.

    MATLAB equivalent: ``gsua_plot('MC_filtering', ...)``. One panel per free parameter.

    Args:
        mcf: Result from :func:`gsua_csb.monte_carlo_filter`.
        m_free: (N, Np_free) design matrix restricted to free columns, in the same column order as
            ``mcf.names`` (i.e. the same ``M[:, ~model.fixed]`` :func:`gsua_csb.monte_carlo_filter`
            was called with) -- the "prior" distribution each panel compares the low/high subsets
            against.
        axes: Flat array of Axes, one per free parameter. A new grid is created if not given.

    Returns:
        The flat array of Axes drawn on.
    """
    import matplotlib.pyplot as plt

    m_free = np.asarray(m_free, dtype=np.float64)
    Np = len(mcf.names)

    if axes is None:
        n_rows = max(1, int(np.floor(np.sqrt(Np))))
        n_cols = n_rows + int(np.ceil((Np - n_rows**2) / n_rows))
        _, axes = plt.subplots(n_rows, n_cols, squeeze=False)
        axes = axes.ravel()
    else:
        axes = np.atleast_1d(axes)

    for k, name in enumerate(mcf.names):
        ax = axes[k]
        _ecdf_step(ax, m_free[:, k], label="Prior")
        _ecdf_step(ax, mcf.low[:, k], label="Low")
        _ecdf_step(ax, mcf.high[:, k], label="High")
        ax.set_title(name)
        ax.set_xlabel("Value")
        ax.set_ylabel("CDF")
    axes[0].legend()
    return axes
