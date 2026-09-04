"""Design-of-experiments sampling: the Python replacement for ``gsua_dmatrix``.

MATLAB's ``gsua_dmatrix`` builds an (N, Np) factor-space design matrix from a model's parameter
table, with a method switch between Latin hypercube (``lhsdesign``, default), plain uniform random
(``makedist('Uniform', ...)``), and scrambled Sobol (``sobolset``/``scramble``). This module reuses
``scipy.stats.qmc`` for the quasi-Monte-Carlo methods instead of MATLAB Statistics Toolbox
equivalents.

``method="joint"`` additionally provides correlation-preserving sampling, drawing whole parameter
vectors from an ensemble of accepted estimates instead of each factor independently. MATLAB
equivalent: ``gsua_dmatrix(T, N, 'Method', 'Joint')``.

Orientation note: MATLAB stores its ensemble as ``Np x nPool`` (``Tia.Est``, one row per parameter,
matching its table layout), while ``pool`` here is ``(n_pool, Np)`` -- one row per estimate, matching
this package's own ``PEResult.x`` and ``identifiability_analysis(estimates=...)`` convention. Each
language follows its own existing convention rather than importing the other's.

Not ported: ``gsua_dmatrix``'s ``rMethod == 'normal'`` branch, which samples each factor from an
independent normal distribution using ``T.Range`` columns as ``(mu, sigma)`` instead of
``(lower, upper)`` bounds -- a distinct convention from every other ``rMethod`` that :func:`gsua_csb.build_range`
does not preserve (it always returns bounds). A model wanting normal-distributed sampling can build
its own draws with ``numpy.random.Generator.normal`` directly.
"""

from __future__ import annotations

import warnings
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.stats import qmc

from ._model import Model

DesignMethod = Literal["latin_hypercube", "uniform", "sobol", "joint"]
JointType = Literal["bootstrap", "smooth_bootstrap", "gaussian"]

_JOINT_TYPES = ("bootstrap", "smooth_bootstrap", "gaussian")


def design_matrix(
    model: Model,
    n: int,
    method: DesignMethod = "latin_hypercube",
    seed: int | np.random.Generator | None = None,
    *,
    joint_type: JointType | None = None,
    pool: ArrayLike | None = None,
    min_pool_n: int | None = None,
    bandwidth: float | None = None,
    clip: ArrayLike | None = None,
) -> NDArray[np.float64]:
    """Sample an (n, Np) factor-space design matrix from ``model``'s parameter ranges.

    Fixed factors (``model.fixed``, i.e. ``range[:, 0] == range[:, 1]``) are held at that constant
    value in every row rather than passed through a sampler -- a degenerate (zero-width) dimension
    would either error or waste sampler budget for no informational gain.

    A free factor with ``model.log_scale[i]`` set is sampled in log10 space (the design is
    stratified/scrambled in log10(range), then exponentiated back), not linear space -- essential,
    not cosmetic, for a factor spanning many orders of magnitude: linear-uniform sampling over e.g.
    ``[1e-13, 1000]`` puts effectively all its mass above 1, so a true value like ``2e-12`` is
    essentially never sampled anywhere close to (a real, observed case -- see
    :func:`gsua_csb.load_petab`, which sets ``log_scale`` automatically from PEtab's
    ``parameterScale`` column).

    Args:
        model: Source of parameter ranges (``model.range``), which factors are fixed
            (``model.fixed``), and which free factors use log-scale sampling (``model.log_scale``).
        n: Number of samples (rows).
        method: ``"latin_hypercube"`` (default) -- stratified space-filling design via
            ``scipy.stats.qmc.LatinHypercube``. ``"uniform"`` -- independent uniform draws per
            factor. ``"sobol"`` -- scrambled Sobol low-discrepancy sequence via
            ``scipy.stats.qmc.Sobol``. ``"joint"`` -- correlation-preserving draws from ``pool``
            (see below); the three others are MARGINAL, sampling each factor independently.
        seed: Seed or ``numpy.random.Generator`` for reproducibility.
        joint_type: Only with ``method="joint"``. ``"bootstrap"`` (default) resamples whole pool
            rows with replacement: exact correlations, no distributional assumption, cannot leave
            the manifold, unaffected by a rank-deficient covariance, but limited to the distinct
            vectors already in the pool. ``"smooth_bootstrap"`` adds a variance-corrected Gaussian
            kernel so draws are new points that stay near the manifold; bandwidth defaults to
            Silverman's multivariate rule. ``"gaussian"`` uses the ensemble mean and covariance via
            eigendecomposition, sampling only non-null eigendirections so a rank-deficient
            covariance (``n_pool`` < number of free parameters) still works -- but it assumes the
            manifold is linear, and on a CURVED ridge it places mass beside it.
        pool: Only with ``method="joint"``. ``(n_pool, Np)`` ensemble of accepted estimates, e.g.
            ``IdentifiabilityResult.estimates_used`` or ``PEResult.x``.
        min_pool_n: Only with ``method="joint"``. Smallest pool accepted for joint sampling
            (default 5, matching ``identifiability_analysis``'s ``min_corr_n``). Below it the call
            warns and falls back to Latin hypercube sampling: correlation from very few points is
            degenerate rather than merely noisy -- two points always correlate at exactly +-1 -- so
            the result would be confidently wrong rather than approximately right.
        bandwidth: Only with ``joint_type="smooth_bootstrap"``. Kernel width; defaults to
            Silverman's multivariate rule.
        clip: Only with ``method="joint"``. Optional ``(Np, 2)`` bounds to clip draws to. Defaults
            to no clipping, deliberately: ``identifiability_analysis`` leaves ``range`` as a
            confidence interval of the MEDIAN of the pool, which narrows as the pool grows and is
            not the spread of values consistent with the data, so clipping joint draws to it would
            reimpose exactly the understated uncertainty joint sampling removes. Pass the original
            model's physical/feasible bounds if clipping is wanted.

    Returns:
        (n, Np) design matrix, one row per sample, columns ordered as ``model.names``.

    Raises:
        ValueError: If ``method`` or ``joint_type`` is unrecognized; if a joint-only option is
            given without ``method="joint"``; if ``method="joint"`` without a usable ``pool``; or
            if a free, log-scale factor has a lower bound that isn't strictly positive.

    Note:
        ``log_scale`` applies to the marginal methods only. Joint draws come from the pool's own
        parameter values, which are already in linear space, so no scale transform is involved.
    """
    joint_only = {
        "joint_type": joint_type,
        "pool": pool,
        "min_pool_n": min_pool_n,
        "bandwidth": bandwidth,
        "clip": clip,
    }
    if method != "joint":
        given = [k for k, v in joint_only.items() if v is not None]
        if given:
            raise ValueError(
                f"{', '.join(sorted(given))} only appl{'ies' if len(given) == 1 else 'y'} to "
                f"method='joint' (got method={method!r})"
            )

    lower = model.range[:, 0]
    upper = model.range[:, 1]
    free = ~model.fixed
    n_free = int(free.sum())

    M = np.tile(lower, (n, 1))
    if n_free == 0:
        return M

    if method == "joint":
        return _joint_matrix(
            model,
            n,
            M,
            free,
            np.random.default_rng(seed),
            "bootstrap" if joint_type is None else joint_type,
            pool,
            5 if min_pool_n is None else min_pool_n,
            bandwidth,
            clip,
        )

    log_scale = np.asarray(getattr(model, "log_scale", np.zeros(model.n_params, dtype=bool)), dtype=bool)
    free_log = log_scale[free]
    if np.any(free_log):
        bad_idx = np.where(free)[0][free_log & (lower[free] <= 0)]
        if bad_idx.size:
            raise ValueError(
                "log_scale factors must have a strictly positive lower bound; failed for "
                f"{[model.names[i] for i in bad_idx]}"
            )

    rng = np.random.default_rng(seed)
    if method == "latin_hypercube":
        u = qmc.LatinHypercube(d=n_free, seed=rng).random(n)
    elif method == "uniform":
        u = rng.random((n, n_free))
    elif method == "sobol":
        u = qmc.Sobol(d=n_free, seed=rng).random(n)
    else:
        raise ValueError(f"Unknown design method: {method!r}")

    with np.errstate(divide="ignore", invalid="ignore"):
        # np.where evaluates both branches eagerly, so log10 of a non-log-scale dimension's bound
        # (which may legitimately be <= 0) still runs -- its result is simply never selected.
        lo_free = np.where(free_log, np.log10(lower[free]), lower[free])
        hi_free = np.where(free_log, np.log10(upper[free]), upper[free])
    scaled = qmc.scale(u, lo_free, hi_free)
    scaled[:, free_log] = 10.0 ** scaled[:, free_log]
    M[:, free] = scaled
    return M


def _cov_factor(pool_free: NDArray[np.float64]) -> NDArray[np.float64]:
    """Square-root factor of the ensemble covariance, restricted to non-null eigendirections.

    Dropping the null directions is what makes a rank-deficient covariance (fewer pool rows than
    free parameters, which is routine for an expensive model) usable instead of an error.
    """
    C = np.atleast_2d(np.cov(pool_free, rowvar=False))
    C = (C + C.T) / 2.0
    d, V = np.linalg.eigh(C)
    tol = float(d.max()) * d.size * np.finfo(float).eps
    keep = d > max(tol, 0.0)
    if not keep.any():
        return np.zeros((C.shape[0], 1))
    return V[:, keep] * np.sqrt(d[keep])


def _joint_matrix(
    model: Model,
    n: int,
    M: NDArray[np.float64],
    free: NDArray[np.bool_],
    rng: np.random.Generator,
    joint_type: JointType,
    pool: ArrayLike | None,
    min_pool_n: int,
    bandwidth: float | None,
    clip: ArrayLike | None,
) -> NDArray[np.float64]:
    """Correlation-preserving sampling from an accepted-estimate ensemble."""
    if joint_type not in _JOINT_TYPES:
        raise ValueError(f"Unknown joint_type: {joint_type!r}; expected one of {_JOINT_TYPES}")
    if pool is None:
        raise ValueError(
            "method='joint' needs an ensemble of accepted estimates; pass pool=... "
            "(e.g. IdentifiabilityResult.estimates_used, or PEResult.x)"
        )

    P = np.atleast_2d(np.asarray(pool, dtype=np.float64))
    if P.ndim != 2 or P.shape[1] != model.n_params:
        raise ValueError(f"pool must be (n_pool, {model.n_params}); got {P.shape}")

    pool_free = P[:, free]
    keep = np.all(np.isfinite(pool_free), axis=1)
    if not keep.any():
        raise ValueError(
            "every pool row holds a non-finite value in at least one free parameter"
        )
    if not keep.all():
        warnings.warn(
            f"{int((~keep).sum())} of {keep.size} pool rows hold non-finite values "
            "and were dropped",
            stacklevel=3,
        )
        pool_free = pool_free[keep]

    n_pool, d_free = pool_free.shape
    if n_pool < min_pool_n:
        warnings.warn(
            f"pool has {n_pool} vectors (fewer than min_pool_n={min_pool_n}); correlation "
            "estimated from so few points is degenerate rather than merely noisy, so a joint "
            "sample would be confidently wrong -- falling back to marginal Latin hypercube "
            "sampling",
            stacklevel=3,
        )
        return design_matrix(model, n, "latin_hypercube", rng)

    if joint_type == "bootstrap":
        X = pool_free[rng.integers(0, n_pool, n)]
    else:
        L = _cov_factor(pool_free)
        mu = pool_free.mean(axis=0)
        Z = rng.standard_normal((n, L.shape[1])) @ L.T
        if joint_type == "smooth_bootstrap":
            h = (
                (4.0 / ((d_free + 2) * n_pool)) ** (1.0 / (d_free + 4))
                if bandwidth is None
                else float(bandwidth)
            )
            base = pool_free[rng.integers(0, n_pool, n)]
            # The 1/sqrt(1+h^2) is what keeps this a smoothing of the ensemble rather than an
            # inflation of it: without it the kernel would scale the covariance by (1+h^2).
            X = mu + (base - mu + h * Z) / np.sqrt(1.0 + h**2)
        else:
            X = mu + Z

    M[:, free] = X

    if clip is not None:
        C = np.asarray(clip, dtype=np.float64)
        if C.shape != (model.n_params, 2):
            raise ValueError(f"clip must be ({model.n_params}, 2); got {C.shape}")
        M = np.clip(M, C[:, 0], C[:, 1])
    return M
