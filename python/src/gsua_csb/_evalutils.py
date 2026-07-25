"""Shared batch-evaluation helpers for the Monte-Carlo-based analyses (SA, UA, MCF, ...).

Every one of these analyses runs a model over many parameter sets and needs the result collapsed
to a consistent (N, Nd) shape regardless of whether the underlying ``Model`` is scalar-output,
single time-series, or multi-state -- this is the one place that shape-collapsing logic lives.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from ._model import Model


def eval_batch(
    model: Model, params: NDArray[np.float64], xdata: NDArray[np.float64] | None, output_index: int = 0
) -> NDArray[np.float64]:
    """Evaluate a batch of parameter sets, collapsing to a plain (N, Nd) array.

    A scalar-output model's ``evaluate_batch`` returns (N,), a single-output time-series model
    returns (N, Nd), and a multi-state model (e.g. ``SymbolicODEModel``) returns (N, n_states, Nd)
    -- ``output_index`` selects which state to use in the last case.
    """
    Y = np.asarray(model.evaluate_batch(params, xdata), dtype=np.float64)
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
