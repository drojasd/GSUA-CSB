"""GSUA-CSB: sensitivity analysis, uncertainty analysis, parameter estimation, and practical
identifiability analysis for symbolic-ODE and user-defined models.

Python port of the MATLAB GSUA-CSB toolbox (https://github.com/drojasd/GSUA-CSB). This package
exposes a clean, idiomatic API (:class:`Model`, :func:`costf`, :func:`median_ci`, ...) as the real
implementation, plus ``gsua_*``-named aliases below for readers coming from the MATLAB toolbox or
its citations, who want literal name parity rather than a lookup table.

Not ported: Simulink-backed models (no Python equivalent -- use the MATLAB Engine API to call
MATLAB directly instead of reimplementing Simulink). Every other MATLAB capability has a Python
port here -- see the project README's status table for the full list and the deliberate deviations
from the MATLAB originals documented in each module.
"""

from __future__ import annotations

from ._costs import costf, costf_multi, coverage_metric, likecost, rcostf
from ._csb import CSBResult, RangeRefinementResult, confidence_subcontour_box, range_refinement
from ._estimation import PEResult, parameter_estimation
from ._identifiability import ClusterInfo, IdentifiabilityResult, identifiability_analysis
from ._likelihood import ProfileLikelihoodResult, profile_likelihood
from ._model import Model, UserFunctionModel, build_range
from ._plotting import (
    plot_identifiability_correlation,
    plot_identifiability_index,
    plot_mcf,
    plot_sensitivity_area,
    plot_sensitivity_bar,
    plot_uncertainty,
)
from ._sampling import design_matrix
from ._sensitivity import SensitivityResult, sensitivity_analysis
from ._stats import band_depth, median_ci
from ._uncertainty import MCFResult, UncertaintyResult, monte_carlo_filter, uncertainty_analysis

__version__ = "0.2.0"

try:
    from ._symbolic import SymbolicODEModel
except ImportError:  # pragma: no cover - exercised only when the `symbolic` extra is absent
    SymbolicODEModel = None  # type: ignore[assignment,misc]

try:
    from ._petab import PEtabProblem, load_petab
    from ._sbml import SBMLODESystem, parse_sbml
except ImportError:  # pragma: no cover - exercised only when the `petab` extra is absent
    PEtabProblem = None  # type: ignore[assignment,misc]
    load_petab = None  # type: ignore[assignment,misc]
    SBMLODESystem = None  # type: ignore[assignment,misc]
    parse_sbml = None  # type: ignore[assignment,misc]

# MATLAB-name-parity aliases -- same functions, for cross-reference with the MATLAB toolbox/citations.
gsua_costf = costf
gsua_rcostf = rcostf
gsua_costfmulti = costf_multi
gsua_likecost = likecost
gsua_covmetric = coverage_metric
gsua_medianci = median_ci
gsua_depth = band_depth
gsua_dmatrix = design_matrix
gsua_sa = sensitivity_analysis
gsua_ua = uncertainty_analysis
gsua_mcf = monte_carlo_filter
gsua_pe = parameter_estimation
gsua_ia = identifiability_analysis
gsua_dia = identifiability_analysis
gsua_oatr = range_refinement
gsua_oatr2 = range_refinement
gsua_csb = confidence_subcontour_box
gsua_likelihood = profile_likelihood
# No gsua_plot alias: MATLAB's single dispatcher is replaced by the plot_* functions above, each
# taking a result dataclass directly instead of a plot-type string + positional-count-dependent args.

__all__ = [
    "Model",
    "UserFunctionModel",
    "SymbolicODEModel",
    "build_range",
    "costf",
    "rcostf",
    "costf_multi",
    "likecost",
    "coverage_metric",
    "median_ci",
    "band_depth",
    "design_matrix",
    "SensitivityResult",
    "sensitivity_analysis",
    "UncertaintyResult",
    "MCFResult",
    "uncertainty_analysis",
    "monte_carlo_filter",
    "PEResult",
    "parameter_estimation",
    "ClusterInfo",
    "IdentifiabilityResult",
    "identifiability_analysis",
    "RangeRefinementResult",
    "range_refinement",
    "CSBResult",
    "confidence_subcontour_box",
    "ProfileLikelihoodResult",
    "profile_likelihood",
    "plot_uncertainty",
    "plot_sensitivity_bar",
    "plot_sensitivity_area",
    "plot_identifiability_correlation",
    "plot_identifiability_index",
    "plot_mcf",
    "gsua_costf",
    "gsua_rcostf",
    "gsua_costfmulti",
    "gsua_likecost",
    "gsua_covmetric",
    "gsua_medianci",
    "gsua_depth",
    "gsua_dmatrix",
    "gsua_sa",
    "gsua_ua",
    "gsua_mcf",
    "gsua_pe",
    "gsua_ia",
    "gsua_dia",
    "gsua_oatr",
    "gsua_oatr2",
    "gsua_csb",
    "gsua_likelihood",
    "SBMLODESystem",
    "parse_sbml",
    "PEtabProblem",
    "load_petab",
]
