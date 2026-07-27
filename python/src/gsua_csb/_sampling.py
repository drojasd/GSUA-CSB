"""Design-of-experiments sampling: the Python replacement for ``gsua_dmatrix``.

MATLAB's ``gsua_dmatrix`` builds an (N, Np) factor-space design matrix from a model's parameter
table, with a method switch between Latin hypercube (``lhsdesign``, default), plain uniform random
(``makedist('Uniform', ...)``), and scrambled Sobol (``sobolset``/``scramble``). This module reuses
``scipy.stats.qmc`` for the quasi-Monte-Carlo methods instead of MATLAB Statistics Toolbox
equivalents.

Not ported: ``gsua_dmatrix``'s ``rMethod == 'normal'`` branch, which samples each factor from an
independent normal distribution using ``T.Range`` columns as ``(mu, sigma)`` instead of
``(lower, upper)`` bounds -- a distinct convention from every other ``rMethod`` that :func:`gsua_csb.build_range`
does not preserve (it always returns bounds). A model wanting normal-distributed sampling can build
its own draws with ``numpy.random.Generator.normal`` directly.
"""

from __future__ import annotations

from typing import Literal

import numpy as np
from numpy.typing import NDArray
from scipy.stats import qmc

from ._model import Model

DesignMethod = Literal["latin_hypercube", "uniform", "sobol"]


def design_matrix(
    model: Model,
    n: int,
    method: DesignMethod = "latin_hypercube",
    seed: int | np.random.Generator | None = None,
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
            ``scipy.stats.qmc.Sobol``.
        seed: Seed or ``numpy.random.Generator`` for reproducibility.

    Returns:
        (n, Np) design matrix, one row per sample, columns ordered as ``model.names``.

    Raises:
        ValueError: If ``method`` is not one of the three supported values, or a free,
            log-scale factor has a lower bound that isn't strictly positive.
    """
    lower = model.range[:, 0]
    upper = model.range[:, 1]
    free = ~model.fixed
    n_free = int(free.sum())

    M = np.tile(lower, (n, 1))
    if n_free == 0:
        return M

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
