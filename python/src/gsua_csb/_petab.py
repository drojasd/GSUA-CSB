"""PEtab problem loading: turn a PEtab parameter-estimation problem into a :class:`gsua_csb.Model`.

No MATLAB equivalent -- this is new capability, not a port. `PEtab <https://petab.readthedocs.io/>`_
is a community standard (SBML + a handful of TSV tables) for specifying "fit this dynamical model
to this data" problems in systems biology, used by tools including AMICI/pyPESTO, COPASI, and
Data2Dynamics, with a curated benchmark collection (`Benchmark-Models-PEtab
<https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab>`_) of real published models
with real data. Loading one directly gives this package interoperability with that ecosystem and,
more immediately, a source of externally-validated test problems instead of only hand-built ones.

Scope -- this is deliberately not a general PEtab importer (there is no full PEtab implementation
in Python outside the reference `petab <https://github.com/PEtab-dev/PEtab>`_ library itself, which
this module does not depend on or wrap). Supported, and confirmed against the two benchmark models
this was built and tested against (``Perelson_Science1996``, ``Boehm_JProteomeRes2014``):

- Exactly one problem, one SBML model, one condition/measurement/observable file set.
- Exactly one simulation condition, with no condition-table parameter overrides (i.e. the SBML
  model's own initial values are used as-is) and no preequilibration.
- ``observableFormula`` as an algebraic expression of species and constant SBML parameters (a bare
  species name, as in Perelson, or a derived multi-species expression, as in Boehm) -- no
  ``observableParameters`` placeholder substitution, since neither confirmed benchmark uses it.
- Parameters: every PEtab parameter table row whose ID matches an SBML constant parameter becomes
  part of the returned :class:`gsua_csb.Model` (free if ``estimate=1``, using ``lowerBound``/
  ``upperBound``/``nominalValue`` directly -- these are always given on the *linear* scale in the
  PEtab format regardless of ``parameterScale``, which only affects how an optimizer's *own* search
  space is transformed; fixed at ``nominalValue`` if ``estimate=0``). Rows with no matching SBML
  parameter (typically noise-model parameters like a measurement's standard deviation) are skipped
  -- this importer does not estimate noise-model parameters, only the dynamical ones.
- Species initial conditions are always treated as fixed (not estimated), from :func:`gsua_csb.
  parse_sbml`'s ``initial_conditions`` (literal or evaluated from an ``<initialAssignment>`` at
  nominal parameter values -- see that function's docstring for the scope this inherits).

``observableTransformation`` (e.g. ``log10``, common for viral-load-scale data) is reported as
metadata but not applied -- the returned model simulates in linear (untransformed) space, and
``ydata`` is the raw linear measurement values from the file. Apply any transform yourself before
fitting if the benchmark's own analysis does.

The internal :class:`~gsua_csb.SymbolicODEModel` is integrated with ``method="LSODA"`` rather than
the package default ``"RK45"``. Systems biology reaction networks (like Boehm, whose rate constants
span roughly 1e-5 to 1e5) are routinely stiff, and an explicit non-stiff solver can become
prohibitively slow or effectively stall on them; LSODA auto-switches between non-stiff (Adams) and
stiff (BDF) methods as needed, which is a safe general-purpose default for PEtab-imported models.
Tolerances are tightened to ``rtol=1e-8``/``atol=1e-10`` (looser than that was enough to blur the
cross-validation against Boehm's published reference simulation out past 1e-4 relative error).
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import sympy as sp
from numpy.typing import NDArray

from ._model import UserFunctionModel
from ._sbml import parse_sbml
from ._symbolic import SymbolicODEModel


@dataclass
class PEtabProblem:
    """Result of :func:`load_petab`.

    Attributes:
        model: A :class:`gsua_csb.UserFunctionModel` wrapping ODE simulation + observable
            computation. ``model.names``/``range``/``nominal`` are ``[state initial conditions
            (all fixed), SBML constant parameters (free or fixed per the PEtab table)]``.
            ``model.evaluate(params, xdata)`` returns ``(n_observables, len(xdata))``.
        xdata: Sorted unique measurement time points across every observable.
        ydata: ``(n_observables, len(xdata))`` measurement values aligned to ``xdata``; ``NaN``
            where a given observable has no measurement at a given time point.
        observable_names: Observable IDs, in the same order as ``ydata``'s first axis.
        observable_transformation: This PEtab problem's ``observableTransformation`` per
            observable (e.g. ``"lin"``, ``"log10"``) -- informational only, see the module
            docstring.
    """

    model: UserFunctionModel
    xdata: NDArray[np.float64]
    ydata: NDArray[np.float64]
    observable_names: list[str]
    observable_transformation: list[str]


def load_petab(yaml_path: str) -> PEtabProblem:
    """Load a PEtab problem (YAML + SBML + TSV tables) into a :class:`PEtabProblem`.

    Args:
        yaml_path: Path to the PEtab problem's YAML configuration file. Every file it references
            is resolved relative to its directory.

    Returns:
        A :class:`PEtabProblem`.

    Raises:
        ImportError: If the ``petab`` extra (``python-libsbml``, ``pandas``, ``PyYAML``) is not
            installed.
        ValueError: If the problem uses a structure outside this loader's scope (multiple
            problems/conditions, condition-table parameter overrides, etc. -- see the module
            docstring) or references a parameter not defined anywhere.
    """
    try:
        import pandas as pd
        import yaml
    except ImportError as exc:  # pragma: no cover - exercised only without the petab extra
        raise ImportError("load_petab requires the 'petab' extra: pip install gsua-csb[petab]") from exc

    yaml_path = Path(yaml_path)
    base = yaml_path.parent
    with open(yaml_path) as f:
        config = yaml.safe_load(f)

    if len(config["problems"]) != 1:
        raise ValueError("load_petab only supports a single-problem PEtab file (see module docstring)")
    problem = config["problems"][0]
    for key in ("sbml_files", "condition_files", "measurement_files", "observable_files"):
        if len(problem[key]) != 1:
            raise ValueError(f"load_petab only supports exactly one {key[:-1]} (got {len(problem[key])})")

    ode_sys = parse_sbml(base / problem["sbml_files"][0])

    conditions = pd.read_csv(base / problem["condition_files"][0], sep="\t")
    if len(conditions) != 1:
        raise ValueError("load_petab only supports a single simulation condition (see module docstring)")
    extra_cols = set(conditions.columns) - {"conditionId", "conditionName"}
    if extra_cols:
        raise ValueError(
            f"Condition-table parameter overrides are not supported (found columns {sorted(extra_cols)})"
        )
    condition_id = conditions["conditionId"].iloc[0]

    param_table = pd.read_csv(base / (config["parameter_file"] if isinstance(config["parameter_file"], str)
                                       else config["parameter_file"][0]), sep="\t")
    param_table = param_table.set_index("parameterId")

    free_state_names = ode_sys.state_names
    n_states = len(free_state_names)
    ic_values = np.array([ode_sys.initial_conditions[name] for name in free_state_names])

    kept_param_names: list[str] = []
    param_lb: list[float] = []
    param_ub: list[float] = []
    param_nom: list[float] = []
    for name in ode_sys.param_names:
        if name not in param_table.index:
            warnings.warn(
                f"SBML parameter {name!r} has no entry in the PEtab parameter table; using its "
                f"SBML nominal value ({ode_sys.param_values[name]}) as a fixed constant.",
                stacklevel=2,
            )
            kept_param_names.append(name)
            param_lb.append(ode_sys.param_values[name])
            param_ub.append(ode_sys.param_values[name])
            param_nom.append(ode_sys.param_values[name])
            continue
        row = param_table.loc[name]
        kept_param_names.append(name)
        nominal = float(row["nominalValue"])
        param_nom.append(nominal)
        if int(row["estimate"]) == 1:
            param_lb.append(float(row["lowerBound"]))
            param_ub.append(float(row["upperBound"]))
        else:
            param_lb.append(nominal)
            param_ub.append(nominal)

    names = free_state_names + kept_param_names
    lb = np.concatenate([ic_values, np.array(param_lb)])
    ub = np.concatenate([ic_values, np.array(param_ub)])
    nominal = np.concatenate([ic_values, np.array(param_nom)])
    model_range = np.column_stack([lb, ub])

    observables = pd.read_csv(base / problem["observable_files"][0], sep="\t").set_index("observableId")
    obs_symtab = {name: sym for name, sym in zip(ode_sys.state_names, ode_sys.state_vars)}
    obs_symtab.update({name: sym for name, sym in zip(ode_sys.param_names, ode_sys.params)})
    observable_names = list(observables.index)
    observable_transformation = [
        str(observables.loc[oid].get("observableTransformation", "lin")) for oid in observable_names
    ]
    obs_funcs = []
    for oid in observable_names:
        formula = str(observables.loc[oid, "observableFormula"])
        expr = sp.sympify(formula, locals=obs_symtab)
        obs_funcs.append(sp.lambdify(list(ode_sys.state_vars) + list(ode_sys.params), expr, modules="numpy"))

    measurements = pd.read_csv(base / problem["measurement_files"][0], sep="\t")
    measurements = measurements[measurements["simulationConditionId"] == condition_id]
    xdata = np.array(sorted(measurements["time"].unique()), dtype=np.float64)
    time_index = {t: i for i, t in enumerate(xdata)}
    ydata = np.full((len(observable_names), len(xdata)), np.nan)
    obs_index = {oid: i for i, oid in enumerate(observable_names)}
    for _, row in measurements.iterrows():
        ydata[obs_index[row["observableId"]], time_index[row["time"]]] = row["measurement"]

    ode_model = SymbolicODEModel(
        odes=ode_sys.odes,
        state_vars=ode_sys.state_vars,
        t=ode_sys.t,
        params=ode_sys.params,
        domain=xdata,
        names=names,
        range=model_range,
        method="LSODA",
        solver_kwargs={"rtol": 1e-8, "atol": 1e-10},
    )

    def func(full_params: NDArray[np.float64], xd: NDArray[np.float64]) -> NDArray[np.float64]:
        states = ode_model.evaluate(full_params, xd)  # (n_states, len(xd))
        true_params = full_params[n_states:]
        out = np.zeros((len(obs_funcs), len(xd)))
        for i, obs_func in enumerate(obs_funcs):
            out[i] = obs_func(*states, *true_params)
        return out

    model = UserFunctionModel(
        func=func,
        names=names,
        range=model_range,
        nominal=nominal,
        domain=xdata,
        output_names=observable_names,
    )

    return PEtabProblem(
        model=model,
        xdata=xdata,
        ydata=ydata,
        observable_names=observable_names,
        observable_transformation=observable_transformation,
    )
