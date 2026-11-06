"""Parametric bootstrap of parameter uncertainty: the Python replacement for ``gsua_bootstrap``.

Refitting the SAME data from many starting points (multistart, :func:`gsua_csb.parameter_estimation`
with ``n > 1``) explores the cost function and defends against local minima, but it collapses to the
same optimum: it does not measure how much the DATA NOISE moves the estimate. On the discrete-SIR
test case, multistart drives S0 to a near-degenerate point (~2982.7) when the truth is 5000. A
parametric bootstrap measures that missing quantity -- it perturbs the DATA (not the start point)
with the model's own observation noise and re-estimates, so the spread of the resulting ensemble is
the parameter uncertainty.

Method:
    1. Take the best pool member (by cost) as theta*; simulate it to get the fitted curves f.
    2. Estimate a pooled observation-noise model from theta*'s residuals -- the SAME dispersion
       estimate :func:`gsua_csb.noise_floor` uses (shared ``_estimate_noise_model``).
    3. Draw ``b`` synthetic datasets from f under that noise model (shared ``_generate_bootstrap``):
       ``"poisson"``, ``"quasipoisson"`` (default) or ``"nbglobal"``, with cumulative-series support.
    4. Re-fit each synthetic dataset with :func:`gsua_csb.parameter_estimation` (same solver/margin/
       alpha, ``n`` starts per replica) and keep each replica's best estimate.
    5. Return the ensemble ``x`` (b, Np) and its ``cost``, ready to pass to
       :func:`gsua_csb.identifiability_analysis` as ``estimates``.

Multistart vs bootstrap: the multistart ensemble answers "did the optimizer find the global
optimum?"; the bootstrap ensemble answers "how much would the estimate change under a fresh draw of
the data's noise?". Pass the bootstrap ensemble to :func:`gsua_csb.identifiability_analysis` for
honest parameter-uncertainty intervals and correlations.
"""

from __future__ import annotations

import warnings
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ._estimation import parameter_estimation
from ._model import Model
from ._noisefloor import NoiseModel, _MODELS, _estimate_noise_model, _generate_bootstrap


@dataclass
class BootstrapResult:
    """Result of :func:`parametric_bootstrap`.

    Attributes:
        names: Parameter names, length Np (``model.names``).
        x: (b, Np) ensemble of best-per-replica estimates. Pass directly to
            :func:`gsua_csb.identifiability_analysis` as ``estimates``.
        cost: (b,) best cost per replica.
        theta_star: (Np,) best pool member used to seed the noise model.
        phi: Pooled Pearson dispersion estimate.
        k: Pooled method-of-moments global-NB dispersion (``inf`` if degenerate).
        noise: Which noise model the synthetic datasets were drawn from.
        quasi_fallback: True if quasi-Poisson degenerated to Poisson draws (``phi <= 1``).
        example_synthetic: (n_outputs, len(xdata)) first synthetic draw, for plotting a
            representative bootstrap dataset against the real data.
        solver: Solver used for the refits.
        margin: Margin used for the refits.
        alpha: Cost exponent used for the refits.
    """

    names: list[str]
    x: NDArray[np.float64]
    cost: NDArray[np.float64]
    theta_star: NDArray[np.float64]
    phi: float
    k: float
    noise: "NoiseModel | list[str]"
    quasi_fallback: bool
    example_synthetic: NDArray[np.float64]
    solver: str
    margin: float
    alpha: float


def parametric_bootstrap(
    model: Model,
    xdata: ArrayLike,
    ydata: ArrayLike,
    estimates: ArrayLike,
    cost: ArrayLike,
    *,
    b: int = 200,
    solver: str = "minimize",
    n: int = 20,
    margin: float = 0.0,
    alpha: float = 2.0,
    noise: "NoiseModel | Sequence[str]" = "quasipoisson",
    cumulative: ArrayLike | None = None,
    seed: int | np.random.Generator | None = None,
    n_jobs: int = 1,
    progress: bool = True,
    solver_kwargs: dict | None = None,
) -> BootstrapResult:
    """Parametric-bootstrap ensemble of parameter estimates for identifiability analysis.

    Args:
        model: The model the pool was fit for.
        xdata: Domain the pool was fit against.
        ydata: (n_outputs, len(xdata)) reference data (may contain NaN), same orientation as
            :func:`gsua_csb.parameter_estimation`.
        estimates: (N, Np) pool of estimated parameter vectors (e.g. ``PEResult.x``).
        cost: (N,) fit cost per pool member (e.g. ``PEResult.cost``). The best member seeds the
            noise model.
        b: Number of bootstrap datasets / refits. Default 200.
        solver: Solver for the refits (see :func:`gsua_csb.parameter_estimation`). Default
            ``"minimize"`` (the MATLAB ``fmincon`` analogue).
        n: Starts per replica (multistart within each refit). Default 20.
        margin: Margin for the refits (pass ``PEResult.margin`` to reuse the pool's objective).
        alpha: Cost exponent for the refits (pass ``PEResult.alpha``).
        noise: Noise model for the synthetic datasets. A single name applies to EVERY row -- the
            COUNT models ``"poisson"``, ``"quasipoisson"`` (default), ``"nbglobal"`` (they assume
            integer-count data). For a joint fit mixing count rows and fraction/continuous rows, pass
            a per-row sequence, adding ``"gaussian"`` (additive noise, std from that row's residuals)
            or ``"none"`` (no perturbation): e.g. ``["quasipoisson", "none"]``. This avoids drawing
            count noise (and zeros) on a 0-1 fraction row. See :func:`gsua_csb.noise_floor`.
        cumulative: (n_outputs,) bool, one flag per fitted output row (default: all False). A
            cumulative row is drawn on the incident series and re-accumulated.
        seed: Seed (or Generator) for the synthetic draws and the per-replica refit seeds.
        n_jobs: Parallel worker processes for the B refits (joblib). Default 1 (serial). ``-1`` uses
            all cores.
        progress: If True (default) and serial, print progress every ~10%.
        solver_kwargs: Extra keyword arguments passed through to each refit.

    Returns:
        A :class:`BootstrapResult`. ``result.x`` (b, Np) goes straight into
        :func:`gsua_csb.identifiability_analysis`.

    Raises:
        ValueError: If ``noise`` is unknown, ``cumulative`` has the wrong length, or ``cost``'s
            length doesn't match ``estimates``.

    Example::

        pe = parameter_estimation(model, xdata, ydata, n=200, solver="minimize", margin=0.3)
        bs = parametric_bootstrap(model, xdata, ydata, pe.x, pe.cost, b=200, n=20,
                                  margin=pe.margin, alpha=pe.alpha, seed=1)
        ia = identifiability_analysis(model, bs.x, cost=bs.cost, cluster=True)
    """
    xdata = np.asarray(xdata, dtype=np.float64)
    ydata = np.atleast_2d(np.asarray(ydata, dtype=np.float64))
    estimates = np.atleast_2d(np.asarray(estimates, dtype=np.float64))
    cost = np.atleast_1d(np.asarray(cost, dtype=np.float64))
    b = int(b)
    n = int(n)

    if cost.shape[0] != estimates.shape[0]:
        raise ValueError(
            f"cost must have {estimates.shape[0]} entries (one per row of estimates), got {cost.shape[0]}"
        )

    n_outputs = ydata.shape[0]
    if cumulative is None:
        cumulative_mask = np.zeros(n_outputs, dtype=bool)
    else:
        cumulative_mask = np.asarray(cumulative, dtype=bool)
        if cumulative_mask.shape[0] != n_outputs:
            raise ValueError(
                f"cumulative must have one entry per fitted output row ({n_outputs}), "
                f"got {cumulative_mask.shape[0]}"
            )

    # Resolve the noise model(s). A single name applies to every row (count models assume integer
    # counts); a sequence gives a PER-ROW model so a joint count+fraction fit can use e.g.
    # ("quasipoisson", "none"). 'gaussian' (additive, std from that row's residuals) and 'none' (no
    # perturbation) are only meaningful per row. See gsua_bootstrap in MATLAB for the rationale.
    bootstrap_models = tuple(_MODELS) + ("gaussian", "none")
    if isinstance(noise, str):
        if noise not in bootstrap_models:
            raise ValueError(f"Unknown noise model: {noise!r}, expected one of {list(bootstrap_models)}")
        noise_list = [noise] * n_outputs
    else:
        noise_list = list(noise)
        if len(noise_list) != n_outputs:
            raise ValueError(
                f"noise must be a single model or one per fitted output row ({n_outputs}), got {len(noise_list)}"
            )
        for m in noise_list:
            if m not in bootstrap_models:
                raise ValueError(f"Unknown noise model: {m!r}, expected one of {list(bootstrap_models)}")

    ibest = int(np.argmin(cost))
    theta_star = estimates[ibest]

    f, mu_clamped, phi, k_global, quasi_fallback, n_valid, p_free = _estimate_noise_model(
        model, xdata, ydata, theta_star, cumulative_mask
    )
    if n_valid - p_free <= 0:
        warnings.warn(
            f"Only {n_valid} valid (non-NaN, mu>0) pooled residuals for {p_free} free "
            "parameters; the noise model may be unreliable.",
            stacklevel=2,
        )
    if quasi_fallback and "quasipoisson" in noise_list:
        warnings.warn(
            f"Pooled Pearson dispersion phi={phi:.6g} <= 1; quasipoisson falls back to Poisson draws.",
            stacklevel=2,
        )
    if np.isinf(k_global) and "nbglobal" in noise_list:
        warnings.warn(
            "Method-of-moments global-NB k is degenerate; nbglobal falls back to Poisson draws.",
            stacklevel=2,
        )

    # Per-row additive-Gaussian std for any 'gaussian' row, from theta*'s residuals on that row.
    sigma_row = np.zeros(n_outputs)
    for i, m in enumerate(noise_list):
        if m == "gaussian":
            sigma_row[i] = float(np.nanstd(ydata[i] - f[i]))

    rng = np.random.default_rng(seed)
    noise_arg = noise if isinstance(noise, str) else noise_list
    yb = _generate_bootstrap(
        rng, noise_arg, mu_clamped, cumulative_mask, ydata, phi, k_global, quasi_fallback, b, sigma_row
    )
    # Distinct, reproducible per-replica seeds drawn from the same stream, so each refit explores
    # different starts yet the whole bootstrap is reproducible from `seed`.
    replica_seeds = rng.integers(0, 2**32 - 1, size=b)

    def _refit(j: int) -> tuple[NDArray[np.float64], float]:
        pe = parameter_estimation(
            model, xdata, yb[:, :, j], n=n, solver=solver, margin=margin, alpha=alpha,
            seed=int(replica_seeds[j]), solver_kwargs=solver_kwargs,
        )
        return pe.x[0], float(pe.cost[0])

    if progress:
        print(f"parametric_bootstrap: {b} refits ({solver}, n={n} per replica, noise={noise})")

    if n_jobs == 1:
        results = []
        next_tick = 0.0
        for j in range(b):
            results.append(_refit(j))
            if progress and (j + 1) / b * 100 >= next_tick:
                print(f"  {(j + 1) / b * 100:3.0f}% ({j + 1}/{b})")
                next_tick += 10
    else:
        from joblib import Parallel, delayed

        results = Parallel(n_jobs=n_jobs)(delayed(_refit)(j) for j in range(b))

    xs = np.array([r[0] for r in results], dtype=np.float64)
    costs = np.array([r[1] for r in results], dtype=np.float64)

    return BootstrapResult(
        names=list(model.names),
        x=xs,
        cost=costs,
        theta_star=theta_star,
        phi=phi,
        k=k_global,
        noise=noise,
        quasi_fallback=quasi_fallback,
        example_synthetic=yb[:, :, 0].copy(),
        solver=solver,
        margin=margin,
        alpha=alpha,
    )
