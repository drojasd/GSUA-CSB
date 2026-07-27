"""Ablation: does the cost-filter + cluster-before-outlier fix matter on real data?

Compares the current identifiability_analysis pipeline (cost filter -> cluster -> outlier, scoped
to the dominant cluster) against a faithful reconstruction of the pre-fix order (outlier removal on
the raw multistart estimates, THEN clustering on whatever survives, no cost filter) on real
experimental measurements from two independent, published PEtab benchmarks. "Faithful" means: the
old order is reconstructed by composing the current public API rather than resurrecting deleted
code -- _remove_multivariate_outliers on the raw estimates first, then identifiability_analysis
with cluster=True/outlier=False on whatever survives, which is exactly what the pre-fix code path
did structurally.

Two deliberately different regimes, both real data:

- Boehm_JProteomeRes2014 (systems biology, 6 free params): 60 multistart fits, most of which
  converge to a bad local optimum. This is the case the fix targets directly -- discrete multiple
  global minima, invisible to the old order because untreated bad fits swamp both the outlier
  detector (56/60 "bad" runs aren't outliers relative to each other) and the clustering step.

- Bertozzi_PNAS2020, California condition (epidemiology, 3 free params): 60 multistart fits that
  all converge to the *same* cost -- not a bug (see the direct least_squares check this was
  verified against: gamma_CA moves to the same value from every start while I0_CA/R0_CA don't move
  at all, a live demonstration of the well-known early-epidemic-curve identifiability limitation --
  only the growth rate is constrained by a short case-count window, not R0 and gamma individually).
  A genuine single-basin, continuous non-identifiability along 2 of 3 dimensions, not a discrete
  multi-modality problem. Included deliberately as a negative control: old and new pipelines should
  agree here (there is no real cluster structure for either to find), and they do.
"""

import time
import warnings
from pathlib import Path

import numpy as np

from gsua_csb import (
    UserFunctionModel,
    design_matrix,
    identifiability_analysis,
    load_petab,
    parameter_estimation,
    plot_identifiability_graph,
)
from gsua_csb._identifiability import _remove_multivariate_outliers

warnings.filterwarnings("ignore")

SEED = 0
FIGDIR = Path(__file__).parent / "figures"
FIGDIR.mkdir(exist_ok=True)


def finite_initial_points(model, xdata, n, seed, oversample=4):
    """design_matrix rows that evaluate to finite output, for a model whose PEtab-derived bounds
    span wide-enough ranges that some random draws hit a numerically pathological ODE regime
    (LSODA failing to integrate) -- a real robustness gap in the wide-bounds regime, not something
    to paper over inside parameter_estimation itself (see its module docstring on why it doesn't
    swallow exceptions), but fine to filter for here since this script's purpose is comparing the
    identifiability pipeline, not stress-testing the solver.
    """
    free_idx = np.where(~model.fixed)[0]
    candidates = design_matrix(model, n * oversample, seed=seed)
    keep = []
    for row in candidates:
        full = model.nominal.copy()
        full[free_idx] = row[free_idx]
        if np.all(np.isfinite(model.evaluate(full, xdata))):
            keep.append(row)
        if len(keep) == n:
            break
    if len(keep) < n:
        raise RuntimeError(f"Only found {len(keep)}/{n} finite starting points")
    return np.array(keep)


def run_ablation(label, model, xdata, ydata, n, safe_starts=False):
    print(f"\n{'#'*70}\n# {label}: {len(model.names)} params, {int((~model.fixed).sum())} free, N={n}\n{'#'*70}")
    t0 = time.time()
    initial_points = finite_initial_points(model, xdata, n, SEED) if safe_starts else None
    pe = parameter_estimation(
        model, xdata, ydata, n=n, solver="least_squares", seed=SEED, initial_points=initial_points
    )
    print(f"multistart PE done in {time.time()-t0:.1f}s")
    print(f"cost range: min={pe.cost.min():.4g} max={pe.cost.max():.4g} median={np.median(pe.cost):.4g}")

    X = pe.x

    new_result = identifiability_analysis(
        model, X, cost=pe.cost, cost_method="rtol", cost_rtol=0.1,
        cluster=True, max_k=5, sil_threshold=0.5,
        outlier=True, outlier_alpha=0.025,
        seed=SEED,
    )

    _, keep_mask = _remove_multivariate_outliers(X, alpha=0.025)
    old_result = identifiability_analysis(
        model, X[keep_mask], cost=None,
        cluster=True, max_k=5, sil_threshold=0.5,
        outlier=False,
        seed=SEED,
    )

    def summarize(sub_label, result, n_removed_outlier):
        width0 = model.range[:, 1] - model.range[:, 0]
        with np.errstate(invalid="ignore", divide="ignore"):
            norm_width = (result.range[:, 1] - result.range[:, 0]) / width0
        print(f"\n=== {sub_label} ===")
        print(f"  runs in: {n}, bad-fit removed: {result.n_bad_fit_removed}, outliers removed: {n_removed_outlier}")
        print(f"  clusters found: {result.cluster.num_clusters}, sizes: {result.cluster.cluster_sizes.tolist()}, "
              f"silhouette: {result.cluster.silhouette}")
        print(f"  identifiability index (mean, free params): {np.nanmean(result.index[~model.fixed]):.4f}")
        print(f"  CI width (mean, normalized, free params): {np.nanmean(norm_width[~model.fixed]):.4f}")

    summarize("NEW pipeline (cost filter -> cluster -> outlier)", new_result, new_result.n_outliers_removed)
    summarize("OLD pipeline (outlier -> cluster, no cost filter)", old_result, int((~keep_mask).sum()))

    print("\n  per-parameter (NEW pipeline), free params only:")
    for i in np.where(~model.fixed)[0]:
        print(f"    {model.names[i]:<20} CI [{new_result.range[i,0]:.4g}, {new_result.range[i,1]:.4g}]"
              f"  (original range {model.range[i].tolist()})  index={new_result.index[i]:.4f}")
    print(f"\n  correlation_reliable (NEW pipeline, min_corr_n default): {new_result.correlation_reliable}")
    return new_result, old_result, X, pe.cost


def save_graph_before_after(label, model, X, cost):
    """The dominant-cluster-size gap, in one figure: correlation computed naively from however
    many points landed in the dominant cluster (min_corr_n=1 -- the pre-fix behavior, since the
    old code had no floor at all) vs. with the min_corr_n floor active (default 5). A 2-point
    dominant cluster (a real, observed outcome here) makes every pairwise correlation exactly +-1
    regardless of any real relationship -- naive computation renders a spurious complete graph;
    the floor correctly reports it as undefined instead (no edges).
    """
    import matplotlib.pyplot as plt

    naive = identifiability_analysis(
        model, X, cost=cost, cost_method="rtol", cost_rtol=0.1, cluster=True, max_k=5,
        sil_threshold=0.5, outlier=True, outlier_alpha=0.025, min_corr_n=1, seed=SEED,
    )
    fixed = identifiability_analysis(
        model, X, cost=cost, cost_method="rtol", cost_rtol=0.1, cluster=True, max_k=5,
        sil_threshold=0.5, outlier=True, outlier_alpha=0.025, min_corr_n=5, seed=SEED,
    )
    print(f"\n  dominant cluster size: {int(np.max(fixed.cluster.cluster_sizes))}")
    print(f"  naive (min_corr_n=1) correlation_reliable: {naive.correlation_reliable}")
    print(f"  floored (min_corr_n=5) correlation_reliable: {fixed.correlation_reliable}")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
    plot_identifiability_graph(naive, ax=axes[0])
    axes[0].set_title("Before fix (min_corr_n=1): naive correlation")
    plot_identifiability_graph(fixed, ax=axes[1])
    axes[1].set_title("After fix (min_corr_n=5): flagged unreliable")
    fig.suptitle(f"{label}: dominant-cluster correlation degeneracy")
    fig.tight_layout()
    out = FIGDIR / f"{label}_correlation_graph_before_after.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  saved {out}")


# --- Boehm_JProteomeRes2014 (systems biology, 6 free params) ---
boehm = load_petab("tests/data/petab/Boehm_JProteomeRes2014/Boehm_JProteomeRes2014.yaml")
boehm_new, boehm_old, boehm_X, boehm_cost = run_ablation(
    "Boehm_JProteomeRes2014", boehm.model, boehm.xdata, boehm.ydata, n=60
)
if boehm_new.cluster.num_clusters > 1:
    save_graph_before_after("Boehm_JProteomeRes2014", boehm.model, boehm_X, boehm_cost)

# --- Bertozzi_PNAS2020, California condition (epidemiology, 3 free params) ---
# Only y_I_CA is measured in this condition (y_I_NY row is all-NaN); wrap the model down to the
# one observable that actually has data so least_squares gets finite residuals everywhere.
bertozzi = load_petab("tests/data/petab/Bertozzi_PNAS2020/problem.yaml")["u_CA"]
obs_idx = int(np.where(~np.all(np.isnan(bertozzi.ydata), axis=1))[0][0])
single_obs_model = UserFunctionModel(
    func=lambda params, xd: bertozzi.model.evaluate(params, xd)[[obs_idx]],
    names=bertozzi.model.names, range=bertozzi.model.range, nominal=bertozzi.model.nominal,
    domain=bertozzi.xdata, output_names=[bertozzi.observable_names[obs_idx]],
)
run_ablation(
    "Bertozzi_PNAS2020 (u_CA)", single_obs_model, bertozzi.xdata, bertozzi.ydata[[obs_idx]], n=60,
    safe_starts=True,
)
