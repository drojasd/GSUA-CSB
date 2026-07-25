"""GSUA-CSB: sensitivity analysis, uncertainty analysis, parameter estimation, and practical
identifiability analysis for symbolic-ODE and user-defined models.

Python port of the MATLAB GSUA-CSB toolbox (https://github.com/drojasd/GSUA-CSB). This package
exposes a clean, idiomatic API (:class:`Model`, :func:`costf`, :func:`median_ci`, ...) as the real
implementation, plus ``gsua_*``-named aliases below for readers coming from the MATLAB toolbox or
its citations, who want literal name parity rather than a lookup table.

Not ported: Simulink-backed models (no Python equivalent -- use the MATLAB Engine API to call
MATLAB directly instead of reimplementing Simulink). Everything else in the MATLAB toolbox is
in scope; this is an early, partial release -- see the project README for what's implemented so far.
"""

from __future__ import annotations

from ._costs import costf, costf_multi, coverage_metric, likecost, rcostf
from ._estimation import PEResult, parameter_estimation
from ._model import Model, UserFunctionModel, build_range
from ._sampling import design_matrix
from ._sensitivity import SensitivityResult, sensitivity_analysis
from ._stats import band_depth, median_ci
from ._uncertainty import MCFResult, UncertaintyResult, monte_carlo_filter, uncertainty_analysis

__version__ = "0.1.0"

try:
    from ._symbolic import SymbolicODEModel
except ImportError:  # pragma: no cover - exercised only when the `symbolic` extra is absent
    SymbolicODEModel = None  # type: ignore[assignment,misc]

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
]
