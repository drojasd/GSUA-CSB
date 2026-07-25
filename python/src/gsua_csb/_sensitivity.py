"""Global sensitivity analysis: the Python replacement for ``gsua_sa``.

Ports the single-objective (single scalar/time-series output) branch of MATLAB's ``gsua_sa``: the
``Sobol``, ``Jansen``, ``Saltelli``, and ``Xiao`` variance-based estimators (each requiring
``N*(Np/2+1)`` or ``N*(Np+1)`` model evaluations via the classic Saltelli-style ``A``/``B``/``ABi``
sample matrices), plus the simpler ``OAT`` (one-at-a-time) method. Each estimator produces both a
scalar index (``Si``, ``STi``, based on the summed-squared-error cost against a reference output
``y_exp``) and a time-dependent index (``Si_vec``, ``STi_vec``, based on the raw output variance at
each domain point) -- mirroring the toolbox's two complementary views: "which parameter drives
overall fit error" vs. "which parameter drives the output at each point in time."

Not ported (documented gaps, not silent omissions):

- ``'brute-force'`` -- MATLAB's exhaustive method requiring ``N + 2*Np*N**2`` evaluations. Kept out
  because it is prohibitively expensive and superseded by the other four methods for any Np worth
  analyzing.
- Multi-objective SA (MATLAB's ``ParT.Properties.CustomProperties.output`` with more than one
  selected output, scored via ``gsua_costfMulti`` instead of raw SSE) -- callers with a multi-state
  model should select a single state via ``output_index`` and run analyses per state instead.
- MATLAB's ``'OAT'`` branch passes a single-column matrix to ``gsua_pardeval`` in a way that only
  works when the model has exactly one state variable -- almost certainly a latent bug in the
  MATLAB source, not a method to faithfully reproduce. The port here instead evaluates each free
  parameter at its full sampled range while holding every other parameter at ``model.nominal``,
  which is the standard definition of a one-at-a-time sensitivity sweep.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ._model import Model

SensMethod = Literal["sobol", "jansen", "saltelli", "xiao", "oat"]
_VALID_METHODS = {"sobol", "jansen", "saltelli", "xiao", "oat"}


@dataclass
class SensitivityResult:
    """Result of :func:`sensitivity_analysis`.

    Attributes:
        names: Parameter names, length Np (``model.names``).
        method: The method used, one of :data:`SensMethod`.
        Si: (Np,) first-order sensitivity index of the scalar SSE cost against ``y_exp``. Zero for
            fixed parameters (never estimated).
        STi: (Np,) total-order sensitivity index of the scalar SSE cost.
        Si_vec: (Np, Nd) time-dependent first-order index of the raw output at each domain point.
            ``None`` for ``method="oat"`` (MATLAB's OAT branch only ever produced scalar indices).
        STi_vec: (Np, Nd) time-dependent total-order index. ``None`` for ``method="oat"``.
        Y: Stacked model outputs actually evaluated (shape depends on method).
        J: Scalar SSE cost for each row of ``Y`` (vs. ``y_exp``).
    """

    names: list[str]
    method: str
    Si: NDArray[np.float64]
    STi: NDArray[np.float64]
    Si_vec: NDArray[np.float64] | None
    STi_vec: NDArray[np.float64] | None
    Y: NDArray[np.float64]
    J: NDArray[np.float64]


def _eval_batch(
    model: Model, params: NDArray[np.float64], xdata: NDArray[np.float64] | None, output_index: int
) -> NDArray[np.float64]:
    """Evaluate a batch of parameter sets, collapsing to a plain (N, Nd) array.

    Handles the shape variation across ``Model`` subclasses: a scalar-output model's
    ``evaluate_batch`` returns (N,), a single-output time-series model returns (N, Nd), and a
    multi-state model (e.g. ``SymbolicODEModel``) returns (N, n_states, Nd) -- ``output_index``
    selects which state to analyze in the last case.
    """
    Y = np.asarray(model.evaluate_batch(params, xdata), dtype=np.float64)
    if Y.ndim == 1:
        return Y[:, None]
    if Y.ndim == 3:
        return Y[:, output_index, :]
    return Y


def _nominal_output(
    model: Model, xdata: NDArray[np.float64] | None, output_index: int
) -> NDArray[np.float64]:
    y0 = np.asarray(model.evaluate(model.nominal, xdata), dtype=np.float64)
    if y0.ndim == 0:
        return y0[None]
    if y0.ndim == 2:
        return y0[output_index]
    return y0


def _oat(
    model: Model,
    M: NDArray[np.float64],
    xdata: NDArray[np.float64] | None,
    y_exp: NDArray[np.float64],
    output_index: int,
) -> SensitivityResult:
    N, Np = M.shape
    fixed = model.fixed
    Y = _eval_batch(model, M, xdata, output_index)
    J = np.zeros(Np)
    nominal_rows = np.tile(model.nominal, (N, 1))
    for k in range(Np):
        if fixed[k]:
            continue
        M_oat = nominal_rows.copy()
        M_oat[:, k] = M[:, k]
        Y_oat = _eval_batch(model, M_oat, xdata, output_index)
        J[k] = np.var(np.sum((Y_oat - y_exp) ** 2, axis=1))
    VT = np.sum(J)
    Si = J / VT if VT > 0 else np.zeros(Np)
    STi = J.copy()
    return SensitivityResult(
        names=list(model.names), method="oat", Si=Si, STi=STi, Si_vec=None, STi_vec=None, Y=Y, J=J
    )


def sensitivity_analysis(
    model: Model,
    M: ArrayLike,
    xdata: ArrayLike | None = None,
    y_exp: ArrayLike | None = None,
    method: SensMethod = "xiao",
    pod: float = 1.0,
    output_index: int = 0,
) -> SensitivityResult:
    """Estimate global sensitivity indices for ``model``'s free parameters.

    Args:
        model: The model to analyze (``model.range``/``model.fixed``/``model.nominal`` supply the
            factor metadata; free parameters are those with ``fixed[k] == False``).
        M: (N, Np) design matrix from :func:`gsua_csb.design_matrix`, ideally built with
            ``method="sobol"`` or ``"latin_hypercube"``. ``N`` is truncated to even if odd (the
            ``A``/``B`` split requires it), matching MATLAB's behavior.
        xdata: Points to evaluate at. Defaults to ``model.domain``.
        y_exp: Reference output the SSE cost (``Si``/``STi``) is measured against. Defaults to
            ``model.evaluate(model.nominal, xdata)`` (the nominal run), matching MATLAB's default.
        method: ``"xiao"`` (default, the toolbox's own robust-to-outliers estimator, requires
            ``N*(Np/2+1)`` evaluations), ``"sobol"`` (classic, ``N*(Np+1)``), ``"jansen"``,
            ``"saltelli"`` (both ``N*(Np/2+1)``), or ``"oat"`` (one-at-a-time, non-global, cheapest
            but ignores parameter interactions).
        pod: Exponent for the Xiao method's norm (``0 < pod <= 2``). Ignored for other methods.
        output_index: For a multi-state model whose ``evaluate`` returns (n_states, Nd), which
            state to analyze. Ignored for single-output models.

    Returns:
        A :class:`SensitivityResult`.

    Raises:
        ValueError: If ``method`` is not one of the five supported values, or ``pod`` is out of
            ``(0, 2]``.
    """
    method = method.lower()  # type: ignore[assignment]
    if method not in _VALID_METHODS:
        raise ValueError(f"Unknown sensitivity method: {method!r}, expected one of {sorted(_VALID_METHODS)}")
    if not (0 < pod <= 2):
        raise ValueError(f"pod must be in (0, 2], got {pod}")

    M = np.asarray(M, dtype=np.float64)
    N = M.shape[0]
    if N % 2 != 0:
        N -= 1
        M = M[:N]

    d = model.domain if xdata is None else np.asarray(xdata, dtype=np.float64)

    if y_exp is None:
        y_exp = _nominal_output(model, d, output_index)
    else:
        y_exp = np.atleast_1d(np.asarray(y_exp, dtype=np.float64))

    if method == "oat":
        return _oat(model, M, d, y_exp, output_index)

    Np = model.n_params
    fixed = model.fixed
    Nd = y_exp.shape[0]

    A = M[: N // 2]
    B = M[N // 2 :]
    YA = _eval_batch(model, A, d, output_index)
    YB = _eval_batch(model, B, d, output_index)
    Y = np.vstack([YA, YB])
    JA = np.sum((YA - y_exp) ** 2, axis=1)
    JB = np.sum((YB - y_exp) ** 2, axis=1)
    J = np.concatenate([JA, JB])

    V_vec = np.var(Y, axis=0)
    VJ = np.var(J)

    Si = np.zeros(Np)
    STi = np.zeros(Np)
    Si_vec = np.zeros((Np, Nd))
    STi_vec = np.zeros((Np, Nd))

    if method == "sobol":
        f02_vec = np.mean(Y, axis=0) ** 2
        f02 = np.mean(J) ** 2
    elif method == "saltelli":
        Y2 = (Y - y_exp) ** 2
        V_vec2 = np.var(Y2, axis=0)
        YA2 = (YA - y_exp) ** 2
        YB2 = (YB - y_exp) ** 2
    elif method == "xiao":
        den = np.mean(np.linalg.norm(YA - YB, axis=0) ** pod)
        denv = np.mean(np.abs(YA - YB) ** pod, axis=0)

    for k in range(Np):
        if fixed[k]:
            continue
        ABi_k = A.copy()
        ABi_k[:, k] = B[:, k]
        YABi_k = _eval_batch(model, ABi_k, d, output_index)
        JABi_k = np.sum((YABi_k - y_exp) ** 2, axis=1)

        if method == "sobol":
            BAi_k = B.copy()
            BAi_k[:, k] = A[:, k]
            YBAi_k = _eval_batch(model, BAi_k, d, output_index)
            JBAi_k = np.sum((YBAi_k - y_exp) ** 2, axis=1)
            Si_vec[k] = (np.mean(YA * YBAi_k, axis=0) - f02_vec) / V_vec
            STi_vec[k] = np.mean(YA * (YA - YABi_k), axis=0) / V_vec
            Si[k] = (np.mean(JA * JBAi_k) - f02) / VJ
            STi[k] = np.mean(JA * (JA - JABi_k)) / VJ
        elif method == "jansen":
            Si_vec[k] = 1 - np.mean((YB - YABi_k) ** 2, axis=0) / (2 * V_vec)
            STi_vec[k] = np.mean((YA - YABi_k) ** 2, axis=0) / (2 * V_vec)
            Si[k] = 1 - np.mean((JB - JABi_k) ** 2) / (2 * VJ)
            STi[k] = np.mean((JA - JABi_k) ** 2) / (2 * VJ)
        elif method == "saltelli":
            YABi2 = (YABi_k - y_exp) ** 2
            Si_vec[k] = np.mean(YB2 * (YABi2 - YA2), axis=0) / V_vec2
            STi_vec[k] = np.mean((YA2 - YABi2) ** 2, axis=0) / (2 * V_vec2)
            Si[k] = np.mean(JB * (JABi_k - JA)) / VJ
            STi[k] = np.mean((JA - JABi_k) ** 2) / (2 * VJ)
        elif method == "xiao":
            Si_vec[k] = (denv - np.mean(np.abs(YB - YABi_k) ** pod, axis=0)) / denv
            STi_vec[k] = np.mean(np.abs(YA - YABi_k) ** pod, axis=0) / denv
            Si[k] = (den - np.mean(np.linalg.norm(YB - YABi_k, axis=0) ** pod)) / den
            STi[k] = np.mean(np.linalg.norm(YA - YABi_k, axis=0) ** pod) / den

    return SensitivityResult(
        names=list(model.names), method=method, Si=Si, STi=STi, Si_vec=Si_vec, STi_vec=STi_vec, Y=Y, J=J
    )
