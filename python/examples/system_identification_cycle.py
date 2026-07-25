"""Worked example: the GSUA-CSB semi-automation system-identification cycle.

This script demonstrates the toolbox's core workflow end to end, using only the Python API built
so far: given a model and data, decide whether the model is well identified, and if it isn't, take
a principled corrective action instead of just reporting failure.

    Phase 1 -- Reachability check
        Before trusting any single estimate, sample the *current* parameter box and check whether
        the real data even falls inside the range of plausible model behavior. If it doesn't, no
        amount of optimization will fix that -- the model structure or ranges need to change first.

    Phase 2 -- Parameter estimation (multistart)
        Fit the model repeatedly from independent starting points. One estimate alone can't tell
        you whether the optimizer found *the* answer or *an* answer.

    Phase 3 -- Practical identifiability analysis
        Turn the repeated estimates into per-parameter confidence intervals and an identifiability
        index. A parameter whose repeated estimates scatter across most of its allowed range, or
        that is strongly correlated with another parameter, is not something the data actually
        pins down -- no matter how good the fit looks.

    Phase 4 -- Confidence-interval thinness check
        Two independent questions, both worth asking before declaring victory: does the identified
        region *fit the data* (coverage_metric's cost_data), and is the identified region *tight*
        (cost_band)? A model can fit beautifully while remaining practically unidentifiable, and it
        can be tightly converged while still not fitting the data well -- they fail for different
        reasons and call for different fixes.

    Phase 5 -- Corrective action and repeat
        When a parameter is flagged poorly identified (usually because it's structurally entangled
        with another one), the standard fix is not "run the optimizer harder" -- it's to fix that
        parameter at its best estimate and re-run the cycle on the remaining free parameters. This
        script does exactly that and shows the previously-poorly-identified parameters become well
        identified once the entangled one is fixed.

The running example is a deliberately non-identifiable model: y(t) = a * b * exp(-k * t). The data
only constrains the *product* a*b, not a and b individually -- a classic structural
non-identifiability that shows up constantly in real kinetic/pharmacological models (e.g. a
combined rate constant that's actually a product of two physical quantities). Watch the
correlation between a and b in Phase 3, and the fix in Phase 5.

Run with: python examples/system_identification_cycle.py
"""

from __future__ import annotations

import numpy as np

from gsua_csb import (
    UserFunctionModel,
    confidence_subcontour_box,
    coverage_metric,
    design_matrix,
    identifiability_analysis,
    parameter_estimation,
    uncertainty_analysis,
)

RNG = np.random.default_rng(0)


def _model(a: float, b: float, k: float) -> np.ndarray:
    def func(p, d):
        return p[0] * p[1] * np.exp(-p[2] * d)

    return UserFunctionModel(
        func=func,
        names=["a", "b", "k"],
        range=np.array([[0.3, 3.0], [0.3, 3.0], [0.1, 1.2]]),
        nominal=np.array([a, b, k]),
        domain=np.linspace(0, 5, 20),
    )


def section(title: str) -> None:
    print(f"\n{'=' * 70}\n{title}\n{'=' * 70}")


def main() -> None:
    xdata = np.linspace(0, 5, 20)

    # "True" system: amplitude a*b = 2.0, decay k = 0.5. Note a and b are individually made up --
    # only their product is a real, data-grounded quantity.
    true_model = _model(a=1.0, b=2.0, k=0.5)
    ydata = true_model.evaluate(true_model.nominal, xdata) + RNG.normal(0, 0.03, size=xdata.shape)

    # Start the cycle from a wide, deliberately uninformative guess.
    model = _model(a=1.0, b=1.0, k=0.6)

    # ---- Phase 1: Reachability check --------------------------------------------------------
    section("Phase 1 -- Reachability check")
    # coverage_metric's cost_data/cost_band are scaled to a tight (default 10%) tolerance band --
    # meaningful as a *post-convergence* thinness check (Phase 4), but not a useful pass/fail gate
    # against the wide prior range still in play here. For a coarse "is the data even in the right
    # ballpark" check before spending any optimizer budget, a simple containment fraction is more
    # robust: what fraction of the real data points fall inside the Monte Carlo ensemble's P5-P95
    # band at all?
    M0 = design_matrix(model, 500, seed=RNG)
    ua0 = uncertainty_analysis(model, M0, y_exp=ydata)
    _, _, p_low0, _, p_high0 = coverage_metric(ua0.Y, ydata)
    contained = np.mean((ydata >= p_low0[0]) & (ydata <= p_high0[0]))
    print(f"{contained:.0%} of the real data points fall inside the current P5-P95 band.")
    if contained < 0.5:
        raise SystemExit("Most of the data falls outside the plausible band -- widen the ranges or "
                          "revisit the model structure before spending any optimizer budget.")

    # ---- Phase 2: Parameter estimation (multistart) ------------------------------------------
    section("Phase 2 -- Parameter estimation (multistart, n=8)")
    pe = parameter_estimation(model, xdata, ydata, n=8, solver="least_squares", seed=RNG.integers(1 << 30))
    print("Best fit:", dict(zip(pe.names, np.round(pe.x[0], 3))), f"cost={pe.cost[0]:.4g}")
    print(f"{pe.x.shape[0]} independent estimates collected for identifiability analysis.")

    # ---- Phase 3: Practical identifiability analysis ------------------------------------------
    section("Phase 3 -- Practical identifiability analysis")
    ia = identifiability_analysis(model, pe.x, cost=pe.cost, correction=True)
    if ia.n_bad_fit_removed:
        print(f"Dropped {ia.n_bad_fit_removed} multistart run(s) that converged to a much worse "
              f"fit than the best one -- their parameters would have contaminated the statistics below.")
    for name, rng, idx in zip(ia.names, ia.range, ia.index):
        print(f"  {name:>2}: range=[{rng[0]:.3f}, {rng[1]:.3f}]  index={idx:.3f}")
    print("Correlation matrix:")
    print(np.round(ia.correlation, 2))
    # The composite `index` blends interval width with correlation, so a pair can be strongly
    # entangled without necessarily crossing an arbitrary index cutoff (their *joint* uncertainty
    # is what's poorly constrained, not each one's marginal interval). A strong pairwise
    # correlation is the more direct signal for *this specific* failure mode -- two parameters
    # trading off against each other -- so that's what drives the corrective action below.
    entangled_pairs = [
        (ia.names[i], ia.names[j], ia.correlation[i, j])
        for i in range(len(ia.names))
        for j in range(i + 1, len(ia.names))
        if abs(ia.correlation[i, j]) > 0.9
    ]
    for n1, n2, r in entangled_pairs:
        print(f"\ncorr({n1}, {n2}) = {r:.3f} -- strongly entangled: the data constrains their "
              f"combination, not either one individually.")

    # ---- Phase 4: Confidence-interval thinness check ------------------------------------------
    section("Phase 4 -- Confidence-interval thinness check")
    model.range = ia.range
    model.nominal = ia.nominal
    M1 = design_matrix(model, 500, seed=RNG)
    ua1 = uncertainty_analysis(model, M1, y_exp=ydata)
    cost_data1, cost_band1, *_ = coverage_metric(ua1.Y, ydata)
    print(f"cost_data = {cost_data1:.3f} (accuracy), cost_band = {cost_band1:.3f} (precision)")
    print("Both are on gsua_costf's tight-tolerance scale (below 1 = within a 10% shape tolerance) --")
    print("still above 1 here because Phase 3 alone narrowed the range, not a dedicated tightening")
    print("pass; that's exactly what Phase 5's correction and the CSB step below are for.")

    # ---- Phase 5: Corrective action and repeat -------------------------------------------------
    section("Phase 5 -- Corrective action: fix one of each entangled pair")
    if entangled_pairs:
        # 'a' and 'b' are symmetric in the model (only a*b matters) -- fix whichever has the wider
        # normalized interval, the standard "give up on estimating this one" move.
        widths = (ia.range[:, 1] - ia.range[:, 0]) / (model.range[:, 1] - model.range[:, 0])
        n1, n2, _ = entangled_pairs[0]
        fix_name = n1 if widths[ia.names.index(n1)] >= widths[ia.names.index(n2)] else n2
        fix_idx = ia.names.index(fix_name)
        print(f"Fixing '{fix_name}' at its estimated value and re-running the cycle on the rest.")
        model.fix(fix_idx)

        pe2 = parameter_estimation(model, xdata, ydata, n=8, solver="least_squares", seed=RNG.integers(1 << 30))
        ia2 = identifiability_analysis(model, pe2.x, cost=pe2.cost, correction=True)
        print("\nAfter fixing:")
        for name, rng, idx in zip(ia2.names, ia2.range, ia2.index):
            print(f"  {name:>2}: range=[{rng[0]:.3f}, {rng[1]:.3f}]  index={idx:.3f}")
        max_remaining_corr = np.nanmax(np.abs(ia2.correlation - np.eye(len(ia2.names))))
        print(f"Largest remaining pairwise |correlation|: {max_remaining_corr:.3f} (was "
              f"{abs(entangled_pairs[0][2]):.3f} before fixing '{fix_name}')")

        recovered_product = pe2.x[0][0] * pe2.x[0][1]
        print(f"\nRecovered a*b = {recovered_product:.3f} (true value: 2.0)")
        print(f"Recovered k   = {pe2.x[0][2]:.3f} (true value: 0.5)")
    else:
        print("Nothing to fix -- no strongly entangled parameter pairs found.")

    # ---- Optional: tighten the surviving free parameters' ranges via the CSB algorithm --------
    section("Bonus -- Confidence sub-contour box on the reduced model")
    csb = confidence_subcontour_box(model, n=300, y_exp=ydata, reps=25, lim=0.3, stop=0.6, seed=0)
    print(f"CSB converged: {csb.converged}")
    for name, rng in zip(csb.names, csb.range):
        print(f"  {name:>2}: [{rng[0]:.3f}, {rng[1]:.3f}]")


if __name__ == "__main__":
    main()
