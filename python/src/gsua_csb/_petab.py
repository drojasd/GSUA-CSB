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
this module does not depend on or wrap). Supported, and confirmed against the benchmark models this
was built and tested against (``Perelson_Science1996``, ``Boehm_JProteomeRes2014``,
``Giordano_Nature2020``, ``Okuonghae_ChaosSolitonsFractals2020``, ``Bertozzi_PNAS2020``):

- Exactly one problem, one SBML model, one condition/measurement/observable file set.
- One or more simulation conditions, no preequilibration (a non-empty
  ``preequilibrationConditionId`` column raises). A single-condition problem returns one
  :class:`PEtabProblem`; a multi-condition problem (e.g. one model fit jointly across several
  regions or intervention scenarios, common in epidemiology PEtab files) returns a
  ``{conditionId: PEtabProblem}`` dict, one independent model per condition.
- Condition-table overrides (a column per overridden species or parameter id, common for
  region/scenario-specific dynamics): a numeric cell fixes that target at that literal value for
  the condition; a cell containing another parameter's id (string) makes that referenced
  parameter -- with its own ``parameters.tsv`` bounds, free or fixed -- the thing being
  estimated/fixed for that target, in that condition (e.g. Bertozzi's ``I0_`` column holding
  ``"I0_CA"``/``"I0_NY"``: the same SBML parameter position, a different estimated quantity per
  region). If the override target feeds an SBML ``<initialAssignment>`` and is itself a reference
  to a *free* parameter, the initial condition is computed once from that parameter's PEtab
  ``nominalValue`` (see :func:`gsua_csb.parse_sbml`'s docstring -- the same pre-existing limitation
  as any other parameter feeding an initial assignment, just also reachable through a condition
  override now).
- ``observableFormula`` as an algebraic expression of species and constant SBML parameters (a bare
  species name, as in Perelson, or a derived multi-species expression, as in Boehm) -- no
  ``observableParameters`` placeholder substitution, since no confirmed benchmark uses it.
- Parameters: every PEtab parameter table row whose ID matches an SBML constant parameter (either
  directly, or indirectly through a condition-table override, see above) becomes part of the
  returned :class:`gsua_csb.Model` (free if ``estimate=1``, using ``lowerBound``/``upperBound``/
  ``nominalValue`` directly -- these are always given on the *linear* scale in the PEtab format
  regardless of ``parameterScale``; fixed at ``nominalValue`` if ``estimate=0``). A free
  parameter's ``parameterScale`` (``"log10"`` or otherwise) is read and becomes that position's
  :attr:`Model.log_scale`, so :func:`gsua_csb.design_matrix`/:func:`gsua_csb.parameter_estimation`
  sample and search it in log space -- essential, not optional, for the epidemiology/systems-
  biology rate constants this importer routinely encounters, which can span many orders of
  magnitude (e.g. Okuonghae's ``theta``: bound ``[1e-13, 1000]``, true value ``2e-12`` -- linear-
  uniform sampling essentially never lands anywhere near it). Rows with no matching SBML parameter
  and no condition-table reference (typically noise-model parameters like a measurement's standard
  deviation) are skipped -- this importer does not estimate noise-model parameters, only the
  dynamical ones.
- Species initial conditions are always treated as fixed (not estimated), from :func:`gsua_csb.
  parse_sbml`'s ``initial_conditions`` (literal, evaluated from an ``<initialAssignment>`` at
  nominal/overridden parameter values, or a direct condition-table override -- see that function's
  docstring for the scope this inherits).

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
    """Result of :func:`load_petab` for a single simulation condition (:func:`load_petab` returns
    a ``{conditionId: PEtabProblem}`` dict instead when the PEtab problem has more than one).

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


def _resolve_condition_overrides(
    cond_row, override_cols: list[str], param_table, condition_id: str
) -> tuple[dict[str, float], dict[str, str]]:
    """Split one condition-table row's override cells into literal values and free-parameter
    aliases (see the module docstring for the two forms an override cell can take)."""
    import pandas as pd

    literal: dict[str, float] = {}
    free_aliases: dict[str, str] = {}
    for col in override_cols:
        cell = cond_row[col]
        if pd.isna(cell):
            continue
        try:
            literal[col] = float(cell)
            continue
        except (TypeError, ValueError):
            pass
        ref_id = str(cell)
        if ref_id not in param_table.index:
            raise ValueError(
                f"Condition {condition_id!r} column {col!r} references parameter {ref_id!r}, "
                "which has no entry in the PEtab parameter table"
            )
        row = param_table.loc[ref_id]
        if int(row["estimate"]) == 1:
            free_aliases[col] = ref_id
        else:
            literal[col] = float(row["nominalValue"])
    return literal, free_aliases


def _load_condition_problem(
    sbml_path,
    condition_id: str,
    literal_overrides: dict[str, float],
    free_aliases: dict[str, str],
    param_table,
    observables,
    measurements,
) -> PEtabProblem:
    """Build one :class:`PEtabProblem` for a single condition, given its resolved overrides."""
    # Free-parameter aliases still need a numeric stand-in for initial-assignment recomputation
    # (see parse_sbml's docstring on this pre-existing limitation for any parameter feeding one).
    sbml_overrides = dict(literal_overrides)
    for col, ref_id in free_aliases.items():
        sbml_overrides[col] = float(param_table.loc[ref_id, "nominalValue"])
    ode_sys = parse_sbml(sbml_path, overrides=sbml_overrides)

    unmatched = (set(literal_overrides) | set(free_aliases)) - set(ode_sys.state_names) - set(
        ode_sys.param_names
    )
    if unmatched:
        raise ValueError(
            f"Condition {condition_id!r} overrides columns {sorted(unmatched)}, which match no "
            "species or parameter in the SBML model"
        )

    free_state_names = ode_sys.state_names
    n_states = len(free_state_names)
    ic_values = np.array([ode_sys.initial_conditions[name] for name in free_state_names])

    kept_param_names: list[str] = []
    param_lb: list[float] = []
    param_ub: list[float] = []
    param_nom: list[float] = []
    param_log_scale: list[bool] = []
    for name in ode_sys.param_names:
        if name in free_aliases:
            ref_id = free_aliases[name]
            row = param_table.loc[ref_id]
            kept_param_names.append(ref_id)
            param_nom.append(float(row["nominalValue"]))
            param_lb.append(float(row["lowerBound"]))
            param_ub.append(float(row["upperBound"]))
            param_log_scale.append(str(row.get("parameterScale", "lin")) == "log10")
        elif name in literal_overrides:
            kept_param_names.append(name)
            param_nom.append(literal_overrides[name])
            param_lb.append(literal_overrides[name])
            param_ub.append(literal_overrides[name])
            param_log_scale.append(False)  # fixed for this condition -- scale is moot
        elif name in param_table.index:
            row = param_table.loc[name]
            kept_param_names.append(name)
            nominal = float(row["nominalValue"])
            param_nom.append(nominal)
            if int(row["estimate"]) == 1:
                param_lb.append(float(row["lowerBound"]))
                param_ub.append(float(row["upperBound"]))
                param_log_scale.append(str(row.get("parameterScale", "lin")) == "log10")
            else:
                param_lb.append(nominal)
                param_ub.append(nominal)
                param_log_scale.append(False)  # fixed -- scale is moot
        else:
            warnings.warn(
                f"SBML parameter {name!r} has no entry in the PEtab parameter table; using its "
                f"SBML nominal value ({ode_sys.param_values[name]}) as a fixed constant.",
                stacklevel=3,
            )
            kept_param_names.append(name)
            param_lb.append(ode_sys.param_values[name])
            param_ub.append(ode_sys.param_values[name])
            param_nom.append(ode_sys.param_values[name])
            param_log_scale.append(False)  # fixed -- scale is moot

    names = free_state_names + kept_param_names
    lb = np.concatenate([ic_values, np.array(param_lb)])
    ub = np.concatenate([ic_values, np.array(param_ub)])
    nominal = np.concatenate([ic_values, np.array(param_nom)])
    model_range = np.column_stack([lb, ub])
    log_scale = np.concatenate([np.zeros(n_states, dtype=bool), np.array(param_log_scale, dtype=bool)])

    obs_symtab = {name: sym for name, sym in zip(ode_sys.state_names, ode_sys.state_vars)}
    obs_symtab.update({name: sym for name, sym in zip(ode_sys.param_names, ode_sys.params)})
    observable_names = list(observables.index)
    observable_transformation = [
        str(observables.loc[oid].get("observableTransformation", "lin")) for oid in observable_names
    ]
    derived_subs = {obs_symtab.get(name, sp.Symbol(name)): expr for name, expr in ode_sys.derived.items()}
    obs_funcs = []
    for oid in observable_names:
        formula = str(observables.loc[oid, "observableFormula"])
        expr = sp.sympify(formula, locals=obs_symtab)
        if derived_subs:
            # An observableFormula can name an assignment-rule-defined "reporter" species directly
            # (e.g. Giordano's observable_CurrentCases = CurrentDiagnosedInfected, a derived
            # aggregate with no ODE of its own) rather than a real state -- resolve it the same way
            # parse_sbml already resolves such references internally.
            expr = expr.subs(derived_subs)
        obs_funcs.append(sp.lambdify(list(ode_sys.state_vars) + list(ode_sys.params), expr, modules="numpy"))

    cond_measurements = measurements[measurements["simulationConditionId"] == condition_id]
    xdata = np.array(sorted(cond_measurements["time"].unique()), dtype=np.float64)
    time_index = {t: i for i, t in enumerate(xdata)}
    ydata = np.full((len(observable_names), len(xdata)), np.nan)
    obs_index = {oid: i for i, oid in enumerate(observable_names)}
    for _, row in cond_measurements.iterrows():
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
        log_scale=log_scale,
    )

    return PEtabProblem(
        model=model,
        xdata=xdata,
        ydata=ydata,
        observable_names=observable_names,
        observable_transformation=observable_transformation,
    )


def load_petab(yaml_path: str) -> PEtabProblem | dict[str, PEtabProblem]:
    """Load a PEtab problem (YAML + SBML + TSV tables) into one or more :class:`PEtabProblem`.

    Args:
        yaml_path: Path to the PEtab problem's YAML configuration file. Every file it references
            is resolved relative to its directory.

    Returns:
        A single :class:`PEtabProblem` if the problem has exactly one simulation condition, or a
        ``{conditionId: PEtabProblem}`` dict if it has more than one (see the module docstring).

    Raises:
        ImportError: If the ``petab`` extra (``python-libsbml``, ``pandas``, ``PyYAML``) is not
            installed.
        ValueError: If the problem uses a structure outside this loader's scope (multiple
            problems, preequilibration, an unresolvable condition-table override -- see the module
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

    sbml_path = base / problem["sbml_files"][0]

    conditions = pd.read_csv(base / problem["condition_files"][0], sep="\t")
    override_cols = [c for c in conditions.columns if c not in ("conditionId", "conditionName")]

    param_table = pd.read_csv(base / (config["parameter_file"] if isinstance(config["parameter_file"], str)
                                       else config["parameter_file"][0]), sep="\t")
    param_table = param_table.set_index("parameterId")

    observables = pd.read_csv(base / problem["observable_files"][0], sep="\t").set_index("observableId")

    measurements = pd.read_csv(base / problem["measurement_files"][0], sep="\t")
    if "preequilibrationConditionId" in measurements.columns:
        preeq = measurements["preequilibrationConditionId"].dropna().astype(str).str.strip()
        if (preeq != "").any():
            raise ValueError("Preequilibration is not supported (see module docstring)")

    condition_ids_with_data = set(measurements["simulationConditionId"])

    results: dict[str, PEtabProblem] = {}
    for _, cond_row in conditions.iterrows():
        condition_id = str(cond_row["conditionId"])
        if condition_id not in condition_ids_with_data:
            # A condition can legitimately have zero measurements -- e.g. a forward-projection
            # "what if" scenario in an epidemiology PEtab file, alongside the one condition that
            # was actually fit to data. Nothing to estimate or evaluate against, so skip it rather
            # than raise; the PEtab problem itself is otherwise perfectly well-formed.
            warnings.warn(
                f"Condition {condition_id!r} has no measurements; skipping (see module docstring)",
                stacklevel=2,
            )
            continue
        literal_overrides, free_aliases = _resolve_condition_overrides(
            cond_row, override_cols, param_table, condition_id
        )
        results[condition_id] = _load_condition_problem(
            sbml_path, condition_id, literal_overrides, free_aliases, param_table, observables, measurements
        )

    if not results:
        raise ValueError("No simulation condition in this PEtab problem has any measurements")
    if len(results) == 1:
        return next(iter(results.values()))
    return results
