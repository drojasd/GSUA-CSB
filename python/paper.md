---
title: 'gsua-csb: A Python Toolbox for Global Sensitivity, Uncertainty, and Practical Identifiability Analysis of Dynamical Models'
tags:
  - Python
  - sensitivity analysis
  - uncertainty quantification
  - parameter estimation
  - identifiability analysis
  - system identification
  - ordinary differential equations
authors:
  - name: Daniel Rojas-Diaz
    orcid: 0000-0000-0000-0000  # TODO: replace with real ORCID iD before submission
    corresponding: true
    affiliation: 1
  - name: Carlos Mario Vélez-Sánchez
    orcid: 0000-0000-0000-0000  # TODO: replace with real ORCID iD before submission
    affiliation: 1
affiliations:
  - name: Universidad EAFIT, Medellín, Colombia
    index: 1
date: 25 July 2026
bibliography: paper.bib
---

<!--
TODO before submission:
- Fill in real ORCID iDs for both authors (or remove the field for whichever author does not have
  one -- it is optional in JOSS, not something to leave as a placeholder).
- Confirm the Zenodo DOI cited for the original MATLAB toolbox (@RojasDiaz2019) -- the form
  currently on record in the repository README, `10.5755281/zenodo.3383316`, does not match
  Zenodo's usual `10.5281/zenodo.<record>` pattern and should be double-checked against the
  actual Zenodo record before this paper is submitted anywhere.
- Double check every citation below resolves to the intended paper (author, venue, year) -- they
  were assembled from established, well-known references but were not individually verified
  against the publisher record at the time of writing.
-->

# Summary

`gsua-csb` is a Python package for validating mathematical models of dynamical systems against
data. Given a model (an ordinary differential equation system or an arbitrary user-defined
function) and a set of uncertain parameters, it answers four questions researchers routinely need
answered together: which parameters actually drive the model's output (global sensitivity
analysis); how parameter uncertainty propagates to output uncertainty (Monte Carlo uncertainty
analysis and filtering); what parameter values best explain observed data (multistart parameter
estimation); and, critically, whether those fitted parameters are actually determined by the data
at all (practical identifiability analysis) [@Raue2009]. A parameter can fit beautifully and still
be practically unidentifiable -- structurally entangled with another parameter so that only their
combination, not either one individually, is constrained by the data. `gsua-csb` treats this as a
first-class outcome rather than an edge case: its identifiability analysis reports a per-parameter
index and a full correlation matrix, and includes spectral-clustering-based detection of multiple
global minima across repeated estimation runs, surfacing when a "good fit" landed in one of several
distinct, equally-valid basins of attraction rather than a single optimum. A worked example
distributed with the package (`examples/system_identification_cycle.py`) demonstrates the complete
cycle end to end -- reachability check, multistart estimation, identifiability analysis,
confidence-interval thinness check, and corrective action -- on a deliberately non-identifiable
model, recovering the true, data-constrained quantity even though the individual entangled
parameters are not recoverable.

# Statement of need

Global sensitivity analysis, uncertainty quantification, parameter estimation, and identifiability
analysis are frequently needed together on the same dynamical-systems model, in epidemiology,
pharmacology, systems biology, and engineering, yet the Python scientific ecosystem does not
currently offer one coherent package spanning all four. `SALib` [@Herman2017] is the standard
Python library for variance-based sensitivity analysis but is intentionally scoped to sensitivity
indices alone, with no built-in model-evaluation layer, parameter estimation, or identifiability
tooling. Bayesian frameworks such as `PyMC` address uncertainty and parameter inference jointly but
require committing to a full probabilistic model and MCMC sampling, a heavier workflow than the
point-estimate-plus-repeated-multistart approach common in applied model calibration. Practical
identifiability analysis via profile likelihood or repeated-estimation clustering has no dedicated,
general-purpose Python package at all as far as the authors are aware; researchers currently
re-implement it ad hoc per project.

`gsua-csb` is a Python port of the MATLAB toolbox of the same name [@RojasDiaz2019], developed at
Universidad EAFIT and, at the time of writing, distributed on MATLAB File Exchange. Porting it to
Python removes the MATLAB license requirement and puts the toolbox's methodology -- most notably
Xiao et al.'s distance-based multivariate sensitivity decomposition [@Xiao2018] alongside the
classical Saltelli/Jansen/Sobol estimators [@Saltelli2010], and the toolbox's namesake Confidence
Sub-contour Box algorithm for iteratively narrowing a parameter box to the region consistent with
observed data -- within reach of the SciPy ecosystem [@Virtanen2020; @Harris2020] and, by extension,
of automated and agentic research tooling that increasingly targets Python over MATLAB. The port
also carries forward, and documents explicitly, several fixes and clarifications found by
re-deriving the MATLAB implementation line by line: a missing correlation-undefined guard in the
core cost function, a unit inconsistency between the acceptance threshold and the inner objective
in the profile-likelihood routine, and a NaN-cascade in the identifiability computation once a
parameter is fixed mid-workflow -- exactly the situation the toolbox's own corrective-action step
produces. The package replaces MATLAB's single struct-with-custom-properties-and-integer-dispatch
data model with a small `Model` class hierarchy (real polymorphism), exposes every analysis result
as a typed dataclass rather than positional outputs, and provides composable Matplotlib plotting
functions [@Hunter2007] instead of one monolithic, string-dispatched plotting function. Symbolic
ODE models are supported via SymPy and SciPy's `solve_ivp` [@Meurer2017; @Virtanen2020]; clustering
uses scikit-learn [@Pedregosa2011]. Simulink-backed models have no Python equivalent and are
intentionally out of scope; users with a Simulink model are expected to call MATLAB directly via
its Python-MATLAB Engine API rather than have Simulink's execution semantics reimplemented.

# Acknowledgements

`gsua-csb` builds directly on the GSUA-CSB MATLAB toolbox and the global sensitivity analysis
methodology it implements, developed at Universidad EAFIT.

# References
