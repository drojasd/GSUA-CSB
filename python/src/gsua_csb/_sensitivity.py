"""Global sensitivity analysis: the Python replacement for ``gsua_sa``.

Ports the single-objective (single scalar/time-series output) branch of MATLAB's ``gsua_sa``: the
``Sobol``, ``Jansen``, ``Saltelli``, and ``Xiao`` variance-based estimators (each requiring
``N*(Np/2+1)`` or ``N*(Np+1)`` model evaluations via the classic Saltelli-style ``A``/``B``/``ABi``
sample matrices), plus the simpler ``OAT`` (one-at-a-time) method. Each estimator produces both a
scalar index (``Si``, ``STi``, based on the summed-squared-error cost against a reference output
``y_exp``) and a time-dependent index (``Si_vec``, ``STi_vec``, based on the raw output variance at
each domain point) -- mirroring the toolbox's two complementary views: "which parameter drives
overall fit error" vs. "which parameter drives the output at each point in time."

What each estimator measures (same as MATLAB): ``sobol`` and ``jansen`` estimate the variance-based
(Sobol') indices of the output; ``xiao`` estimates ``E|YA-YABi|^pod / E|YA-YB|^pod`` (distance
components, Xiao, Lu & Wang 2018) -- with ``pod=2`` it coincides with ``jansen``, with the default
``pod=1`` it is a moment-independent, distribution-based index that can rank factors differently;
``saltelli`` scores the squared-deviation field ``(Y - y_exp)**2`` rather than ``Y``. The scalar
``Si``/``STi`` are always indices of the cost ``J = sum((Y - y_exp)**2)``.

Indices outside [0, 1]: the indices lie in [0, 1] by definition, but their Monte Carlo *estimates*
do not -- with small ``N`` an ``STi`` above 1 or a negative ``Si`` is sampling noise, deliberately not
clipped (clipping would hide that ``N`` is too small). Most frequent for ``saltelli``, rare for
``xiao`` with ``pod=1``. Increase ``N`` until ``sum(Si)/sum(abs(Si))`` is stable near 1.

``method="shapley"`` estimates Shapley effects (Owen 2014; Song, Nelson & Staum 2016): the output
variance shared among the factors with every interaction split fairly between the factors involved.
They sum to 1, satisfy ``Si <= Sh_i <= STi`` for independent factors, and stay meaningful for
DEPENDENT factors. ``shapley_type="independent"`` (default) is the nested Monte Carlo estimator over
random permutations with factors drawn from their marginals (``N + m*(Np-1)*shapley_no*shapley_ni``
evaluations); ``shapley_type="given"`` is the nearest-neighbour estimator of Broto, Bachoc &
Depecker (2020) on a JOINT sample (the rows of ``M``, or ``pool``), keeping its dependence. Results go
to ``Sh``/``Sh_vec`` (``Si``/``STi`` are NaN and the ``*_vec`` Sobol fields ``None``).

``"oat"`` evaluates each free parameter at its sampled values while holding every other parameter at
``model.nominal`` (the standard one-at-a-time sweep; MATLAB's ``'OAT'`` uses the same definition since
v1.19, before which it passed a single column of ``M`` to the model and failed). ``Si`` is each
factor's share of the summed OAT variance of the cost, ``STi`` that variance itself; ``J`` is the cost
of the base sample (one value per row of ``Y``), as for every method.

Not ported (documented gaps, not silent omissions):

- ``'brute-force'`` -- MATLAB's exhaustive method requiring ``N + 2*Np*N**2`` evaluations. Kept out
  because it is prohibitively expensive and superseded by the other four methods for any Np worth
  analyzing.
- Multi-objective SA (MATLAB's ``ParT.Properties.CustomProperties.output`` with more than one
  selected output, scored via ``gsua_costfMulti`` instead of raw SSE) -- callers with a multi-state
  model should select a single state via ``output_index`` and run analyses per state instead.
"""

from __future__ import annotations

import itertools
import math
import warnings
from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.spatial import cKDTree

from ._evalutils import eval_batch, nominal_output, require_model
from ._model import Model

SensMethod = Literal["sobol", "jansen", "saltelli", "xiao", "oat", "shapley"]
_VALID_METHODS = {"sobol", "jansen", "saltelli", "xiao", "oat", "shapley"}
ShapleyType = Literal["independent", "given"]
_SHAPLEY_TYPES = {"independent", "given"}


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
        Sh: (Np,) Shapley effects of the scalar SSE cost (``method="shapley"`` only, else ``None``).
            Zero for fixed parameters; sums to 1.
        Sh_vec: (Np, Nd) time-dependent Shapley effects of the raw output (``method="shapley"``
            only, else ``None``); each column sums to 1.
    """

    names: list[str]
    method: str
    Si: NDArray[np.float64]
    STi: NDArray[np.float64]
    Si_vec: NDArray[np.float64] | None
    STi_vec: NDArray[np.float64] | None
    Y: NDArray[np.float64]
    J: NDArray[np.float64]
    Sh: NDArray[np.float64] | None = None
    Sh_vec: NDArray[np.float64] | None = None


# ---------------------------------------------------------------------------------------------------
# Shapley effects. Value of a coalition J of KNOWN factors: c(J) = E[Var(Y | X_J)], so c({}) = Var(Y)
# and c(all) = 0 for a deterministic model. Sh_i = (1/Var Y) * E_perm[c(P_i) - c(P_i U {i})], P_i the
# factors before i in a uniformly random permutation. The increments telescope from Var(Y) to 0 along
# every permutation, so the indices sum to 1 exactly. Same algorithms as MATLAB gsua_sa.
# ---------------------------------------------------------------------------------------------------


def _marginal_sample(model: Model, n: int, rng: np.random.Generator) -> NDArray[np.float64]:
    """Independent draws from each factor's marginal, matching :func:`gsua_csb.design_matrix`:
    uniform on ``model.range`` (log10-uniform for ``model.log_scale`` factors); fixed factors held."""
    lo, hi = model.range[:, 0], model.range[:, 1]
    fixed = model.fixed
    X = np.empty((n, model.n_params))
    for k in range(model.n_params):
        if fixed[k]:
            X[:, k] = lo[k]
        elif model.log_scale[k]:
            X[:, k] = 10.0 ** rng.uniform(np.log10(lo[k]), np.log10(hi[k]), n)
        else:
            X[:, k] = rng.uniform(lo[k], hi[k], n)
    return X


def _safe_ratio(num: NDArray[np.float64], den: NDArray[np.float64] | float) -> NDArray[np.float64]:
    with np.errstate(divide="ignore", invalid="ignore"):
        return num / den


def _shapley_independent(
    model: Model,
    M: NDArray[np.float64],
    d_x: NDArray[np.float64] | None,
    y_exp: NDArray[np.float64],
    output_index: int,
    n_jobs: int,
    m: int,
    n_outer: int,
    n_inner: int,
    rng: np.random.Generator,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Song, Nelson & Staum (2016) nested Monte Carlo over permutations, independent factors."""
    Np = model.n_params
    fixed = model.fixed
    free = np.flatnonzero(~fixed)
    d = free.size
    Y = eval_batch(model, M, d_x, output_index, n_jobs=n_jobs)
    J = np.sum((Y - y_exp) ** 2, axis=1)
    Nd = Y.shape[1]
    V = np.var(Y, axis=0, ddof=1)
    VJ = np.var(J, ddof=1)
    Sh = np.zeros(Np)
    Sh_vec = np.zeros((Np, Nd))
    if d == 1:
        Sh[free] = 1.0
        Sh_vec[free] = _safe_ratio(V, V)
    if d <= 1:
        return Sh, Sh_vec, Y, J

    # Permutations: every one equally often when they all fit in the budget (stratified), else random.
    n_perm = math.factorial(d)
    if n_perm <= m:
        P_all = np.array(list(itertools.permutations(range(d))))
        reps = m // n_perm
        extra = P_all[rng.choice(n_perm, m - reps * n_perm, replace=False)]
        P = np.vstack([np.tile(P_all, (reps, 1)), extra])
    else:
        P = np.array([rng.permutation(d) for _ in range(m)])

    # Group g = (i*(d-1) + j-1)*n_outer + o holds n_inner rows sharing the outer draw of the
    # coalition free[P[i, :j]] and differing in every other free factor.
    G = m * (d - 1) * n_outer
    known = np.zeros((G, Np), dtype=bool)
    g = 0
    for i in range(m):
        for j in range(1, d):
            known[g : g + n_outer, free[P[i, :j]]] = True
            g += n_outer
    outer = _marginal_sample(model, G, rng)
    X = _marginal_sample(model, G * n_inner, rng)
    grp = np.repeat(np.arange(G), n_inner)
    K = known[grp]
    O = outer[grp]
    X[K] = O[K]

    Yb = eval_batch(model, X, d_x, output_index, n_jobs=n_jobs)
    Jb = np.sum((Yb - y_exp) ** 2, axis=1)
    cv = np.var(Yb.reshape(G, n_inner, Nd), axis=1, ddof=1)
    cv = cv.reshape(m * (d - 1), n_outer, Nd).mean(axis=1).reshape(m, d - 1, Nd)
    cJ = np.var(Jb.reshape(G, n_inner), axis=1, ddof=1)
    cJ = cJ.reshape(m * (d - 1), n_outer).mean(axis=1).reshape(m, d - 1)

    # Increment of the j-th factor of permutation i: c(prefix_{j-1}) - c(prefix_j).
    C_full = np.concatenate([np.broadcast_to(V, (m, 1, Nd)), cv, np.zeros((m, 1, Nd))], axis=1)
    CJ_full = np.concatenate([np.full((m, 1), VJ), cJ, np.zeros((m, 1))], axis=1)
    who = free[P].ravel()
    acc = np.zeros((Np, Nd))
    accJ = np.zeros(Np)
    np.add.at(acc, who, (C_full[:, :-1] - C_full[:, 1:]).reshape(-1, Nd))
    np.add.at(accJ, who, (CJ_full[:, :-1] - CJ_full[:, 1:]).ravel())
    Sh = accJ / (m * VJ)
    Sh_vec = _safe_ratio(acc, m * V)
    Sh[fixed] = 0.0
    Sh_vec[fixed] = 0.0
    return Sh, Sh_vec, Y, J


def _shapley_given(
    model: Model,
    X: NDArray[np.float64],
    d_x: NDArray[np.float64] | None,
    y_exp: NDArray[np.float64],
    output_index: int,
    n_jobs: int,
    k: int,
    m: int,
    rng: np.random.Generator,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Broto, Bachoc & Depecker (2020) given-data estimator on one joint sample: c(J) is the mean
    variance of the output over each point's k nearest neighbours in the standardized X_J
    coordinates. Exact Shapley formula over all coalitions up to 12 players, random permutations
    above."""
    # A repeated vector adds no information and would zero the nearest-neighbour variance.
    _, first = np.unique(X, axis=0, return_index=True)
    X = X[np.sort(first)]
    n, Np = X.shape
    if n <= k:
        raise ValueError(f"The joint sample has {n} distinct vectors; shapley_type='given' needs more "
                         f"than shapley_k={k}.")
    if n < 100:
        warnings.warn(f"Only {n} distinct vectors in the joint sample: the nearest-neighbour Shapley "
                      "estimates will be noisy. Use a few hundred or more.", stacklevel=3)
    Y = eval_batch(model, X, d_x, output_index, n_jobs=n_jobs)
    J = np.sum((Y - y_exp) ** 2, axis=1)
    Nd = Y.shape[1]
    V = np.var(Y, axis=0, ddof=1)
    VJ = np.var(J, ddof=1)
    sd = X.std(axis=0, ddof=1)
    players = np.flatnonzero(~model.fixed & (sd > 0))   # a constant factor cannot carry variance
    d = players.size
    Sh = np.zeros(Np)
    Sh_vec = np.zeros((Np, Nd))
    if d == 1:
        Sh[players] = 1.0
        Sh_vec[players] = _safe_ratio(V, V)
    if d <= 1:
        return Sh, Sh_vec, Y, J
    Z = (X[:, players] - X[:, players].mean(axis=0)) / sd[players]

    def knn_var(mask: int) -> tuple[NDArray[np.float64], float]:
        cols = [b for b in range(d) if (mask >> b) & 1]
        Zj = Z[:, cols]
        _, idx = cKDTree(Zj).query(Zj, k=k)
        return np.var(Y[idx], axis=1, ddof=1).mean(axis=0), float(np.var(J[idx], axis=1, ddof=1).mean())

    if d <= 12:
        n_sub = 1 << d
        C = np.zeros((n_sub, Nd))
        CJ = np.zeros(n_sub)
        C[0], CJ[0] = V, VJ                     # empty coalition; the full one stays 0
        for mask in range(1, n_sub - 1):
            C[mask], CJ[mask] = knn_var(mask)
        size = np.array([bin(mask).count("1") for mask in range(n_sub)])
        for a in range(d):
            for mask in range(n_sub):
                if (mask >> a) & 1:
                    continue
                w = 1.0 / (d * math.comb(d - 1, int(size[mask])))
                with_a = mask | (1 << a)
                Sh_vec[players[a]] += w * (C[mask] - C[with_a])
                Sh[players[a]] += w * (CJ[mask] - CJ[with_a])
    else:
        memo: dict[int, tuple[NDArray[np.float64], float]] = {}
        zero = np.zeros(Nd)
        for _ in range(m):
            perm = rng.permutation(d)
            mask = 0
            prev, prevJ = V, VJ
            for j, a in enumerate(perm):
                mask |= 1 << int(a)
                if j == d - 1:
                    Cc, CJc = zero, 0.0
                else:
                    if mask not in memo:
                        memo[mask] = knn_var(mask)
                    Cc, CJc = memo[mask]
                Sh_vec[players[a]] += prev - Cc
                Sh[players[a]] += prevJ - CJc
                prev, prevJ = Cc, CJc
        Sh_vec /= m
        Sh /= m
    Sh = Sh / VJ
    Sh_vec = _safe_ratio(Sh_vec, V)
    non_players = np.ones(Np, dtype=bool)
    non_players[players] = False
    Sh[non_players] = 0.0
    Sh_vec[non_players] = 0.0
    return Sh, Sh_vec, Y, J


def _oat(
    model: Model,
    M: NDArray[np.float64],
    xdata: NDArray[np.float64] | None,
    y_exp: NDArray[np.float64],
    output_index: int,
    n_jobs: int = 1,
) -> SensitivityResult:
    N, Np = M.shape
    fixed = model.fixed
    Y = eval_batch(model, M, xdata, output_index, n_jobs=n_jobs)
    J = np.sum((Y - y_exp) ** 2, axis=1)       # cost of the base sample, as for every method
    V_oat = np.zeros(Np)
    nominal_rows = np.tile(model.nominal, (N, 1))
    for k in range(Np):
        if fixed[k]:
            continue
        M_oat = nominal_rows.copy()
        M_oat[:, k] = M[:, k]
        Y_oat = eval_batch(model, M_oat, xdata, output_index, n_jobs=n_jobs)
        V_oat[k] = np.var(np.sum((Y_oat - y_exp) ** 2, axis=1))
    VT = np.sum(V_oat)
    Si = V_oat / VT if VT > 0 else np.zeros(Np)
    STi = V_oat.copy()
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
    output_index: int | None = 0,
    n_jobs: int = 1,
    *,
    shapley_type: ShapleyType = "independent",
    shapley_perm: int | None = None,
    shapley_no: int = 1,
    shapley_ni: int = 3,
    shapley_k: int = 3,
    pool: ArrayLike | None = None,
    seed: int | np.random.Generator | None = None,
) -> "SensitivityResult | list[SensitivityResult]":
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
            ``"saltelli"`` (both ``N*(Np/2+1)``), ``"shapley"`` (Shapley effects, see the module
            docstring and the ``shapley_*`` options), or ``"oat"`` (one-at-a-time, non-global,
            cheapest but ignores parameter interactions).
        pod: Exponent for the Xiao method's norm (``0 < pod <= 2``). ``pod=2`` reproduces the
            variance-based Jansen index; the default ``1`` is moment-independent. Ignored for other
            methods.
        output_index: For a multi-state model whose ``evaluate`` returns (n_states, Nd), which
            state to analyze (default 0). Pass ``None`` to analyze EVERY active output series and get
            a list of per-series :class:`SensitivityResult` back (the per-series time-dependent
            indices that MATLAB's gsua_sa now returns stacked as ``Si_vec`` of shape
            ``(Np, Nd, nSeries)``). Ignored for single-output models.
        n_jobs: Parallel worker processes for the model evaluations, the counterpart of MATLAB
            ``gsua_sa``'s default-on ``'parallel'`` (default here is ``1``, serial; ``-1`` uses all
            cores). Worth setting only when a single ``model.evaluate`` is genuinely expensive.
        shapley_type: ``"independent"`` (default) -- factors drawn independently from their
            marginals (as :func:`gsua_csb.design_matrix` does); ``M`` is used only for ``Var(Y)``.
            ``"given"`` -- a JOINT sample whose dependence is kept: ``pool`` if given, else the rows
            of ``M`` (e.g. ``design_matrix(..., method="joint", joint_type="smooth_bootstrap")``;
            the default ``"bootstrap"`` only repeats pool rows, and repeats are dropped). With an
            accepted-estimate ensemble this measures sensitivity inside the identified region, not
            over the prior.
        shapley_perm: Number of permutations ``m`` (default ``N``). For ``"independent"``, when
            ``Np! <= m`` every permutation is used equally often; for ``"given"`` it is used only
            above 12 free factors (below that the exact formula over all coalitions is used).
        shapley_no: Outer draws per coalition (``"independent"``, default 1).
        shapley_ni: Inner draws per outer draw (``"independent"``, default 3, at least 2).
        shapley_k: Nearest neighbours (``"given"``, default 3, at least 2).
        pool: ``(n_pool, Np)`` joint sample for ``shapley_type="given"`` (rows are vectors, as
            ``PEResult.x``). Repeated vectors are evaluated once.
        seed: Seed or ``numpy.random.Generator`` for the Shapley draws.

    Returns:
        A :class:`SensitivityResult`, or -- when ``output_index is None`` -- a list of them, one per
        active output series.

    Raises:
        ValueError: If ``method`` is not one of the five supported values, or ``pod`` is out of
            ``(0, 2]``.
    """
    require_model(model, "sensitivity_analysis")
    method = method.lower()  # type: ignore[assignment]
    if method not in _VALID_METHODS:
        raise ValueError(f"Unknown sensitivity method: {method!r}, expected one of {sorted(_VALID_METHODS)}")
    if not (0 < pod <= 2):
        raise ValueError(f"pod must be in (0, 2], got {pod}")
    shapley_type = shapley_type.lower()  # type: ignore[assignment]
    if shapley_type not in _SHAPLEY_TYPES:
        raise ValueError(f"Unknown shapley_type: {shapley_type!r}, expected one of {sorted(_SHAPLEY_TYPES)}")
    if shapley_ni < 2 or shapley_k < 2 or shapley_no < 1:
        raise ValueError("shapley_ni and shapley_k must be >= 2 and shapley_no >= 1")

    # Per-series mode: analyze each active output and return a list (parity with MATLAB's native
    # series dimension). Each entry is an ordinary 2-D SensitivityResult, so plotting/consumers are
    # unchanged; y_exp is re-derived per series.
    if output_index is None:
        rng_shared = np.random.default_rng(seed)
        return [
            sensitivity_analysis(model, M, xdata=xdata, method=method, pod=pod,
                                 output_index=oi, n_jobs=n_jobs, shapley_type=shapley_type,
                                 shapley_perm=shapley_perm, shapley_no=shapley_no,
                                 shapley_ni=shapley_ni, shapley_k=shapley_k, pool=pool,
                                 seed=rng_shared)
            for oi in range(len(model.active_output_names))
        ]

    M = np.asarray(M, dtype=np.float64)
    N = M.shape[0]
    if N % 2 != 0:
        N -= 1
        M = M[:N]

    d = model.domain if xdata is None else np.asarray(xdata, dtype=np.float64)

    if y_exp is None:
        y_exp = nominal_output(model, d, output_index)
    else:
        y_exp = np.atleast_1d(np.asarray(y_exp, dtype=np.float64))

    if method == "oat":
        return _oat(model, M, d, y_exp, output_index, n_jobs=n_jobs)

    if method == "shapley":
        rng = np.random.default_rng(seed)
        m = N if shapley_perm is None else int(shapley_perm)
        if shapley_type == "independent":
            Sh, Sh_vec, Y, J = _shapley_independent(model, M, d, y_exp, output_index, n_jobs, m,
                                                    shapley_no, shapley_ni, rng)
        else:
            X = M if pool is None else np.asarray(pool, dtype=np.float64)
            if X.ndim != 2 or X.shape[1] != model.n_params:
                raise ValueError(f"pool must be (n_pool, Np) with Np = {model.n_params}; got {X.shape}")
            Sh, Sh_vec, Y, J = _shapley_given(model, X, d, y_exp, output_index, n_jobs, shapley_k,
                                              m, rng)
        nan = np.full(model.n_params, np.nan)
        return SensitivityResult(names=list(model.names), method="shapley", Si=nan, STi=nan.copy(),
                                 Si_vec=None, STi_vec=None, Y=Y, J=J, Sh=Sh, Sh_vec=Sh_vec)

    Np = model.n_params
    fixed = model.fixed
    Nd = y_exp.shape[0]

    A = M[: N // 2]
    B = M[N // 2 :]
    YA = eval_batch(model, A, d, output_index, n_jobs=n_jobs)
    YB = eval_batch(model, B, d, output_index, n_jobs=n_jobs)
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
        YABi_k = eval_batch(model, ABi_k, d, output_index, n_jobs=n_jobs)
        JABi_k = np.sum((YABi_k - y_exp) ** 2, axis=1)

        if method == "sobol":
            BAi_k = B.copy()
            BAi_k[:, k] = A[:, k]
            YBAi_k = eval_batch(model, BAi_k, d, output_index, n_jobs=n_jobs)
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
