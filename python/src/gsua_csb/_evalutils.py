"""Shared batch-evaluation helpers for the Monte-Carlo-based analyses (SA, UA, MCF, ...).

Every one of these analyses runs a model over many parameter sets and needs the result collapsed
to a consistent (N, Nd) shape regardless of whether the underlying ``Model`` is scalar-output,
single time-series, or multi-state -- this is the one place that shape-collapsing logic lives.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from ._model import Model


def require_model(model: object, func_name: str) -> None:
    """Raise a clear error if ``model`` is not a :class:`Model`.

    Guards the argument-order trap under the MATLAB-parity aliases: MATLAB calls
    ``gsua_sa(M, T)`` / ``gsua_ua(M, T)`` with the design matrix first, but this package puts the
    model first (``sensitivity_analysis(model, M)``). Both arguments are positional, so a swap
    would otherwise run and return nonsense; catch it with a message that names the fix.
    """
    if not isinstance(model, Model):
        raise TypeError(
            f"{func_name}() expects the Model first, then the design matrix: "
            f"{func_name}(model, M). Got {type(model).__name__} as the first argument -- if you "
            "are porting a MATLAB call, note the order is reversed there (gsua_sa(M, T))."
        )


def eval_batch(
    model: Model,
    params: NDArray[np.float64],
    xdata: NDArray[np.float64] | None,
    output_index: int = 0,
    n_jobs: int = 1,
) -> NDArray[np.float64]:
    """Evaluate a batch of parameter sets, collapsing to a plain (N, Nd) array.

    A scalar-output model's ``evaluate_batch`` returns (N,), a single-output time-series model
    returns (N, Nd), and a multi-state model (e.g. ``SymbolicODEModel``) returns (N, n_states, Nd)
    -- ``output_index`` selects which state to use in the last case. ``n_jobs`` is forwarded to
    :meth:`Model.evaluate_batch` (parallel per-row evaluation; default serial).
    """
    Y = np.asarray(model.evaluate_batch(params, xdata, n_jobs=n_jobs), dtype=np.float64)
    if Y.ndim == 1:
        return Y[:, None]
    if Y.ndim == 3:
        return Y[:, output_index, :]
    return Y


def nominal_output(
    model: Model, xdata: NDArray[np.float64] | None, output_index: int = 0
) -> NDArray[np.float64]:
    """The model's output at ``model.nominal``, collapsed the same way as :func:`eval_batch`."""
    y0 = np.asarray(model.evaluate(model.nominal, xdata), dtype=np.float64)
    if y0.ndim == 0:
        return y0[None]
    if y0.ndim == 2:
        return y0[output_index]
    return y0
