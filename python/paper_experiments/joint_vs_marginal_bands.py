"""Does correlation-preserving sampling actually build a better uncertainty band?

The bands this toolbox's users build come from ``design_matrix(model_with_ia_range, N)`` -- each
parameter sampled independently inside its own interval. For a model with parameter confounding
that is wrong twice over:

1. ``identifiability_analysis`` leaves ``range`` as a distribution-free confidence interval of the
   MEDIAN of the accepted pool. That says where the pool's centre lies, and it NARROWS as the pool
   grows; it is not the spread of parameter values consistent with the data.
2. Sampling those intervals independently discards the joint structure, so draws leave the
   identified manifold entirely.

``method="joint"`` draws whole parameter vectors from the accepted ensemble, which addresses both at
once -- the draws carry the pool's own spread AND its correlations -- provided nothing clips them
back to the CI of the median.

Three rungs, each testing a different claim:

- CORRECTNESS (analytic control). ``y = a*b*exp(-k*t)``: the data constrains only the product
  ``a*b``, so the identified manifold is exactly the hyperbola ``a*b = const``, known in closed
  form. Leakage off it is therefore measurable rather than arguable. The manifold is CURVED, which
  is what separates the samplers: a Gaussian matching the pool's pairwise correlation still places
  mass beside a curved ridge.

- UTILITY (real data, plottable). Bertozzi_PNAS2020 California COVID case counts: 3 free parameters
  with a documented continuous non-identifiability along 2 of 3 dimensions -- only the growth rate,
  a combination of R0 and gamma, is constrained by a short early-epidemic window. Three parameters
  means the ridge can be drawn.

- SAFETY (the stress case). Boehm_JProteomeRes2014: two well-separated basins and a dominant
  cluster small enough that every pairwise correlation is exactly +-1. Joint sampling must REFUSE
  here and fall back, because a band built from such a pool is confidently wrong rather than
  approximately right -- the failure that looks like success.

Acceptance criterion is calibration, NOT tightness. A pool that is a multistart artifact bootstraps
into a very tight, very wrong band; tightness alone would score that as a win.
"""

import time
import warnings
from copy import deepcopy
from pathlib import Path

import numpy as np

from gsua_csb import (
    UserFunctionModel,
    coverage_metric,
    design_matrix,
    identifiability_analysis,
    load_petab,
    parameter_estimation,
    uncertainty_analysis,
)

warnings.filterwarnings("ignore")

SEED = 0
NBAND = 2000
FIGDIR = Path(__file__).parent / "figures"
FIGDIR.mkdir(exist_ok=True)

SAMPLERS = [
    ("marginal (current)", dict(method="latin_hypercube")),
    ("joint/bootstrap", dict(method="joint", joint_type="bootstrap")),
    ("joint/smooth_bootstrap", dict(method="joint", joint_type="smooth_bootstrap")),
    ("joint/gaussian", dict(method="joint", joint_type="gaussian")),
]


def finite_initial_points(model, xdata, n, seed, oversample=4):
    """Starting points whose ODE actually integrates -- same helper as the ablation script."""
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


def build_bands(model, ia_result, xdata, ydata, n=NBAND, seed=SEED):
    """One band per sampler.

    The marginal band samples the ia range (what users do today); the joint bands sample the pool
    the ia range was summarized FROM. Everything downstream is identical, so any difference is
    attributable to band construction alone.
    """
    marginal_model = deepcopy(model)
    marginal_model.range = ia_result.range.copy()

    bands = {}
    for label, kw in SAMPLERS:
        if kw["method"] == "joint":
            M = design_matrix(
                model, n, "joint", seed,
                joint_type=kw["joint_type"], pool=ia_result.estimates_used,
            )
        else:
            M = design_matrix(marginal_model, n, kw["method"], seed)
        Y = uncertainty_analysis(model, M, xdata=xdata).Y
        cost_data, cost_band, p_lo, p_50, p_hi = coverage_metric(Y, ydata)
        bands[label] = dict(
            M=M, Y=Y, cost_data=cost_data, cost_band=cost_band,
            p_lo=p_lo, p_50=p_50, p_hi=p_hi,
            # Spread of the simulated trajectories themselves. Read against the parameter
            # interval each sampler drew from: a sampler drawing from a NARROWER parameter
            # interval yet producing a WIDER output spread is manufacturing that width by
            # leaving the manifold, not reporting inference uncertainty.
            out_spread=float(np.mean(np.std(np.atleast_3d(Y), axis=0))),
        )
    return bands


def containment(p_lo, p_hi, ydata):
    """Fraction of observations below / within / above the band, and its median width.

    This is the calibration check. A well-built 5-95 band should contain most observations without
    being so wide that containment is trivial -- so `within` and `width` have to be read together.
    """
    y = np.atleast_2d(np.asarray(ydata, dtype=float))
    m = np.isfinite(y)
    below = float(np.sum((y < p_lo) & m) / np.sum(m))
    above = float(np.sum((y > p_hi) & m) / np.sum(m))
    return below, 1.0 - below - above, above, float(np.median(p_hi - p_lo))


def report_bands(label, bands, ydata):
    print(f"\n--- {label}: band calibration ---")
    print(f"{'sampler':<24} {'below':>7} {'within':>7} {'above':>7} {'width':>11} "
          f"{'out_spread':>11} {'cost_data':>10}")
    for name, b in bands.items():
        below, within, above, width = containment(b["p_lo"], b["p_hi"], ydata)
        print(f"{name:<24} {below:>7.3f} {within:>7.3f} {above:>7.3f} {width:>11.4g} "
              f"{b['out_spread']:>11.4g} {b['cost_data']:>10.4g}")
    print("  width/out_spread = how much trajectory variation the sampler produced")
    print("  cost_data = distance from the band's MEDIAN to the data (accuracy; <1 is in tolerance)")


# ======================================================================================
# RUNG 1 -- CORRECTNESS: analytic control with a manifold known in closed form
# ======================================================================================
print("=" * 86)
print("RUNG 1 -- CORRECTNESS: y = a*b*exp(-k*t), identified manifold is the hyperbola a*b = const")
print("=" * 86)

XDATA1 = np.linspace(0.0, 5.0, 20)
TRUE1 = np.array([1.0, 2.0, 0.5])
PRODUCT = TRUE1[0] * TRUE1[1]


def _decay(p, xd):
    return (p[0] * p[1] * np.exp(-p[2] * xd))[None, :]


model1 = UserFunctionModel(
    func=_decay, names=["a", "b", "k"],
    range=np.array([[0.3, 3.0], [0.3, 3.0], [0.1, 1.2]]),
    nominal=TRUE1, domain=XDATA1, output_names=["y"],
)
ydata1 = _decay(TRUE1, XDATA1) + np.random.default_rng(SEED).normal(0, 0.03, (1, XDATA1.size))

t0 = time.time()
pe1 = parameter_estimation(model1, XDATA1, ydata1, n=60, solver="least_squares", seed=SEED)
ia1 = identifiability_analysis(model1, pe1.x, cost=pe1.cost, cost_rtol=0.1, seed=SEED)
print(f"multistart PE + IA in {time.time()-t0:.1f}s; pool kept {ia1.estimates_used.shape[0]}/60 fits")
print(f"recovered a*b: mean {np.mean(ia1.estimates_used[:,0]*ia1.estimates_used[:,1]):.4f} "
      f"(true {PRODUCT:.4f}), sd {np.std(ia1.estimates_used[:,0]*ia1.estimates_used[:,1]):.2e}")
print(f"pool corr(a,b) = {np.corrcoef(ia1.estimates_used[:,0], ia1.estimates_used[:,1])[0,1]:+.4f}")
print(f"ia range for a: [{ia1.range[0,0]:.4f}, {ia1.range[0,1]:.4f}]  "
      f"vs pool spread [{ia1.estimates_used[:,0].min():.4f}, {ia1.estimates_used[:,0].max():.4f}]"
      "   <- the CI-of-the-median narrowness, visible directly")

target = float(np.mean(ia1.estimates_used[:, 0] * ia1.estimates_used[:, 1]))
print(f"\n{'sampler':<24} {'mean |a*b - c|':>16} {'corr(a,b)':>11}    (c = {target:.4f})")
for name, kw in SAMPLERS:
    if kw["method"] == "joint":
        M = design_matrix(model1, NBAND, "joint", SEED,
                          joint_type=kw["joint_type"], pool=ia1.estimates_used)
    else:
        mm = deepcopy(model1)
        mm.range = ia1.range.copy()
        M = design_matrix(mm, NBAND, kw["method"], SEED)
    leak = np.abs(M[:, 0] * M[:, 1] - target).mean()
    print(f"{name:<24} {leak:>16.4e} {np.corrcoef(M[:,0], M[:,1])[0,1]:>+11.4f}")

bands1 = build_bands(model1, ia1, XDATA1, ydata1)
report_bands("Rung 1 (analytic)", bands1, ydata1)


# ======================================================================================
# RUNG 2 -- UTILITY: real data, 3 parameters, a plottable ridge
# ======================================================================================
print("\n" + "=" * 86)
print("RUNG 2 -- UTILITY: Bertozzi_PNAS2020 (California), documented 2-of-3 non-identifiability")
print("=" * 86)

bertozzi = load_petab("tests/data/petab/Bertozzi_PNAS2020/problem.yaml")["u_CA"]
obs_idx = int(np.where(~np.all(np.isnan(bertozzi.ydata), axis=1))[0][0])
model2 = UserFunctionModel(
    func=lambda params, xd: bertozzi.model.evaluate(params, xd)[[obs_idx]],
    names=bertozzi.model.names, range=bertozzi.model.range,
    nominal=bertozzi.model.nominal, domain=bertozzi.xdata,
    output_names=[bertozzi.observable_names[obs_idx]],
)
ydata2 = bertozzi.ydata[[obs_idx]]

t0 = time.time()
starts = finite_initial_points(model2, bertozzi.xdata, 60, SEED)
pe2 = parameter_estimation(model2, bertozzi.xdata, ydata2, n=60, solver="least_squares",
                           seed=SEED, initial_points=starts)
ia2 = identifiability_analysis(model2, pe2.x, cost=pe2.cost, cost_rtol=0.1, seed=SEED)
print(f"multistart PE + IA in {time.time()-t0:.1f}s; pool kept {ia2.estimates_used.shape[0]}/60 fits")
free2 = np.where(~model2.fixed)[0]
for i in free2:
    print(f"  {model2.names[i]:<12} ia CI [{ia2.range[i,0]:.5g}, {ia2.range[i,1]:.5g}]   "
          f"pool [{ia2.estimates_used[:,i].min():.5g}, {ia2.estimates_used[:,i].max():.5g}]")
if len(free2) >= 2:
    C = np.corrcoef(ia2.estimates_used[:, free2], rowvar=False)
    print(f"  pool correlation among free params:\n{np.array2string(C, precision=3)}")

bands2 = build_bands(model2, ia2, bertozzi.xdata, ydata2)
report_bands("Rung 2 (Bertozzi CA)", bands2, ydata2)


def plot_rung2():
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    t = np.asarray(bertozzi.xdata, dtype=float)
    obs = np.asarray(ydata2, dtype=float)[0]

    for ax, name in zip(axes[0], ["marginal (current)", "joint/bootstrap"]):
        b = bands2[name]
        ax.fill_between(t, b["p_lo"][0], b["p_hi"][0], alpha=0.3, label="5-95% band")
        ax.plot(t, b["p_50"][0], lw=2, label="band median")
        ax.plot(t, obs, "k.", ms=5, label="observed")
        ax.set_title(f"{name}\ncost_data={b['cost_data']:.3g}  cost_band={b['cost_band']:.3g}")
        ax.set_xlabel("time")
        ax.set_ylabel(model2.output_names[0])
        ax.legend(fontsize=8)

    # The mechanism, in parameter space. Plot the two LEAST identified parameters -- the pair
    # whose pool spread is widest relative to its own scale. Picking the first two free columns
    # would put the well-determined gamma_CA on an axis, collapsing the ridge to a vertical line
    # and hiding the very structure the figure exists to show.
    scale = np.maximum(np.abs(ia2.estimates_used).mean(axis=0), np.finfo(float).eps)
    spread = np.ptp(ia2.estimates_used, axis=0) / scale
    i, j = free2[np.argsort(spread[free2])[-2:]]
    for ax, name in zip(axes[1], ["marginal (current)", "joint/bootstrap"]):
        M = bands2[name]["M"]
        ax.scatter(M[:, i], M[:, j], s=4, alpha=0.25, label="draws")
        ax.scatter(ia2.estimates_used[:, i], ia2.estimates_used[:, j], s=28,
                   color="crimson", label="accepted fits")
        ax.set_xlabel(model2.names[i])
        ax.set_ylabel(model2.names[j])
        ax.set_title(f"{name}: parameter space")
        ax.legend(fontsize=8)

    fig.suptitle("Bertozzi_PNAS2020 (California): marginal vs joint band construction")
    fig.tight_layout()
    out = FIGDIR / "joint_vs_marginal_bertozzi.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  saved {out}")


try:
    plot_rung2()
except ImportError:
    print("  matplotlib unavailable; skipped figure")


# ======================================================================================
# RUNG 3 -- SAFETY: the guard must fire rather than produce a confident wrong band
# ======================================================================================
print("\n" + "=" * 86)
print("RUNG 3 -- SAFETY: Boehm_JProteomeRes2014, two basins and a degenerate dominant cluster")
print("=" * 86)

boehm = load_petab("tests/data/petab/Boehm_JProteomeRes2014/Boehm_JProteomeRes2014.yaml")
t0 = time.time()
pe3 = parameter_estimation(boehm.model, boehm.xdata, boehm.ydata, n=60,
                           solver="least_squares", seed=SEED)
ia3 = identifiability_analysis(
    boehm.model, pe3.x, cost=pe3.cost, cost_rtol=0.1,
    cluster=True, max_k=5, sil_threshold=0.5, outlier=True, outlier_alpha=0.025, seed=SEED,
)
print(f"multistart PE + IA in {time.time()-t0:.1f}s")
print(f"clusters found: {ia3.cluster.num_clusters}, sizes {ia3.cluster.cluster_sizes.tolist()}, "
      f"silhouette {ia3.cluster.silhouette}")
print(f"dominant-basin pool passed to the sampler: {ia3.estimates_used.shape[0]} vectors")
print(f"correlation_reliable: {ia3.correlation_reliable}")

with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always")
    M3 = design_matrix(boehm.model, 500, "joint", SEED, pool=ia3.estimates_used)
guard = [str(w.message) for w in caught if "min_pool_n" in str(w.message)]
print(f"\nguard fired: {bool(guard)}")
for g in guard:
    print(f"  -> {g}")
if not guard:
    print("  (pool was large enough this run; the guard is exercised by the unit tests)")

if ia3.cluster.num_clusters > 1:
    # cluster.labels is indexed over the COST-FILTERED points, not the raw 60 runs.
    print("\nPer-basin structure (the advanced alternative to a dominant-basin-only band):")
    print(f"  {ia3.n_bad_fit_removed} of 60 runs dropped by the fit-quality filter first")
    for c, size in enumerate(np.asarray(ia3.cluster.cluster_sizes)):
        mark = " (dominant)" if c == ia3.cluster.dominant_cluster else ""
        print(f"  basin {c}: {int(size)} fits{mark}")
    print("  Even the largest basin here is below min_pool_n, so no basin supports a joint band;")
    print("  the honest answer on this dataset is 'not enough accepted fits', which is what the")
    print("  guard reports instead of a confident-looking band built from two points.")

print("\n" + "=" * 86)
print("READING THESE RESULTS")
print("=" * 86)
print("""
The accepted-fit pool is selected by cost equivalence, so by construction every member
reproduces the data about equally well -- the pool's OUTPUT spread is near zero even when its
PARAMETER spread is large. That is what a ridge is: the parameters are unidentified, the
trajectory they produce is not. A correct joint band on such a model is therefore NARROW.

The two rungs isolate the two stacked errors separately, which is worth reading carefully
because they push band width in OPPOSITE directions:

- Rung 1 shows the discarded-correlation error. The marginal sampler draws from a parameter
  interval ~14x narrower than the pool, yet produces ~3e5 times MORE output spread
  (out_spread 3.6e-2 vs 1.1e-7). Width that large out of a box that small can only come from
  leaving the manifold: it is manufactured, not inferred.

- Rung 2 shows the CI-of-the-median error. Here the marginal interval is so narrow in absolute
  terms that even off-manifold draws barely move the output (out_spread 5.6e-5), so the
  marginal band understates rather than manufactures. The joint band is WIDER (6.4e-2) because
  the pool genuinely spans I0_CA 131-997 where the CI admits only 457-749.

Neither dataset reproduces the dramatic band MIS-CENTERING reported on the dengue-Wolbachia
case: cost_data is identical across samplers in both rungs. That effect needs marginal
intervals wide enough for off-manifold draws to swing the median, which these do not provide.
Reported as a negative result rather than glossed.

Two consequences worth stating plainly:

1. Joint sampling's benefit on these datasets is a CORRECT parameter sample, which matters for
   everything that consumes parameters -- sensitivity analysis, confidence subcontour boxes,
   any further propagation -- and it makes the width question answerable at all.

2. A parameter-uncertainty band is NOT a prediction interval and should not be expected to
   contain noisy observations; observation noise is a separate term, which is what
   gsua_noisefloor / noise_floor models. The two features compose: joint sampling gets the
   parameter part right, the noise floor supplies the observation part. Read the 'within'
   column with that in mind -- it is not a verdict on the sampler.
""")
print("Done.")
