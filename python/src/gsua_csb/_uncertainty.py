"""Uncertainty analysis and Monte Carlo filtering: the Python replacement for ``gsua_ua``/``gsua_MCF``.

MATLAB's ``gsua_ua`` is a thin wrapper: it runs the model over every row of a design matrix (via
``gsua_pardeval``) and immediately hands the result to ``gsua_plot``/``gsua_MCF`` for plotting --
there is no separate "just give me the Monte Carlo output" data step. Here that data step is
:func:`uncertainty_analysis`, returning the raw ensemble instead of plotting it (plotting is a
separate, later module -- see the package README's status table).

MATLAB's ``gsua_MCF`` mixes two concerns: splitting the ensemble into behavioral ("low"/"high",
relative to a reference cost) subsets, and plotting their empirical CDFs against the prior. This
module ports only the data half as :func:`monte_carlo_filter` -- the classification a modeler
actually reasons about -- and drops the plotting.

One deliberate behavioral difference from MATLAB: ``gsua_MCF`` falls back to a single-element
``[min(M(:,i))]``/``[max(M(:,i))]`` array when a "low" or "high" subset is empty, purely so
``cdfplot`` (an early-return MATLAB plotting call) doesn't crash on an empty array. Since this
module doesn't plot, that fallback would silently fabricate a data point that was never actually
in the behavioral subset -- so an empty subset here stays a genuine empty array.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ._evalutils import eval_batch, nominal_output
from ._model import Model


@dataclass
class UncertaintyResult:
    """Result of :func:`uncertainty_analysis`.

    Attributes:
        Y: (N, Nd) Monte Carlo ensemble of model outputs, one row per row of the input design
            matrix ``M``.
        y_nom: (Nd,) reference output (``y_exp`` if given, otherwise the nominal run) that ``Y``
            can be compared/filtered against, e.g. via :func:`monte_carlo_filter`.
        xdata: (Nd,) domain points the ensemble was evaluated at.
    """

    Y: NDArray[np.float64]
    y_nom: NDArray[np.float64]
    xdata: NDArray[np.float64]


def uncertainty_analysis(
    model: Model,
    M: ArrayLike,
    xdata: ArrayLike | None = None,
    y_exp: ArrayLike | None = None,
    output_index: int = 0,
) -> UncertaintyResult:
    """Run a Monte Carlo ensemble of ``model`` over every row of a design matrix.

    Args:
        model: The model to evaluate.
        M: (N, Np) design/sample matrix, e.g. from :func:`gsua_csb.design_matrix`.
        xdata: Points to evaluate at. Defaults to ``model.domain``.
        y_exp: Reference output to attach to the result (``.y_nom``) for later comparison.
            Defaults to ``model.evaluate(model.nominal, xdata)``, matching MATLAB's default.
        output_index: For a multi-state model whose ``evaluate`` returns (n_states, Nd), which
            state to analyze. Ignored for single-output models.

    Returns:
        An :class:`UncertaintyResult` with the full (N, Nd) ensemble.
    """
    M = np.asarray(M, dtype=np.float64)
    d = model.domain if xdata is None else np.asarray(xdata, dtype=np.float64)

    if y_exp is None:
        y_nom = nominal_output(model, d, output_index)
    else:
        y_nom = np.atleast_1d(np.asarray(y_exp, dtype=np.float64))

    Y = eval_batch(model, M, d, output_index)
    return UncertaintyResult(Y=Y, y_nom=y_nom, xdata=np.asarray(d) if d is not None else None)


@dataclass
class MCFResult:
    """Result of :func:`monte_carlo_filter`.

    Attributes:
        names: Names of the free (non-fixed) parameters, in column order matching ``low``/``high``.
        low: (N_low, Np_free) parameter sets whose cost fell below the reference (``J < J_exp``).
        high: (N_high, Np_free) parameter sets whose cost fell above the reference (``J > J_exp``).
        J: (N,) scalar cost used for the split -- summed-output cost, or the output at a single
            time index if ``t_index`` was given.
        j_exp: The reference cost the split was made against.
    """

    names: list[str]
    low: NDArray[np.float64]
    high: NDArray[np.float64]
    J: NDArray[np.float64]
    j_exp: float


def monte_carlo_filter(
    model: Model,
    M: ArrayLike,
    Y: ArrayLike,
    y_exp: ArrayLike,
    t_index: int | None = None,
) -> MCFResult:
    """Split a Monte Carlo ensemble into "low"/"high" parameter subsets relative to a reference.

    For each run, a scalar cost is computed (summed output, or the output at ``t_index`` if
    given) and compared against the same summary of ``y_exp``. Runs below the reference go into
    ``low``, runs above go into ``high`` -- the parameter regions driving low vs. high model
    output, the question MATLAB's Monte Carlo filtering plot answers visually.

    Args:
        model: Source of parameter names and which are fixed (``model.fixed``) -- fixed parameters
            carry no information for filtering and are excluded from the output columns.
        M: (N, Np) design/sample matrix that produced ``Y`` (e.g. the same ``M`` passed to
            :func:`uncertainty_analysis`).
        Y: (N, Nd) Monte Carlo ensemble of model outputs (e.g. ``UncertaintyResult.Y``).
        y_exp: (Nd,) reference output (e.g. ``UncertaintyResult.y_nom``).
        t_index: If given, filter on the output at this single domain index instead of the
            summed output across the whole domain.

    Returns:
        An :class:`MCFResult`.
    """
    M = np.asarray(M, dtype=np.float64)
    Y = np.asarray(Y, dtype=np.float64)
    y_exp = np.atleast_1d(np.asarray(y_exp, dtype=np.float64))

    fixed = model.fixed
    M_free = M[:, ~fixed]
    names = [n for n, f in zip(model.names, fixed) if not f]

    if t_index is None:
        j_exp = float(np.sum(y_exp))
        J = np.sum(Y, axis=1) if Y.ndim > 1 else Y.copy()
    else:
        j_exp = float(y_exp[t_index])
        J = Y[:, t_index]

    low = M_free[J < j_exp]
    high = M_free[J > j_exp]

    return MCFResult(names=names, low=low, high=high, J=J, j_exp=j_exp)
