"""A scoped SBML reaction-network parser producing explicit SymPy ODEs.

There is no MATLAB equivalent to port here -- this exists purely to make :func:`gsua_csb.
load_petab` possible, since PEtab problems define their dynamics as SBML reaction networks
(species + reactions + stoichiometry + kinetic laws), not the explicit ``dx/dt = ...`` SymPy
expressions :class:`gsua_csb.SymbolicODEModel` expects.

This is **not** a general SBML importer. Full SBML support (function definitions, events, delays,
algebraic rules, spatial/hierarchical compartments, unit conversion, `hasOnlySubstanceUnits`
species, multiple simultaneous conditions, ...) is a large, separate undertaking better served by
dedicated tools (libSBML itself, AMICI, COPASI, tellurium). This module supports exactly the
subset needed to import the reaction-network SBML models used by the PEtab benchmark collection's
simpler problems (confirmed against ``Perelson_Science1996`` and ``Boehm_JProteomeRes2014``):

- Reactions with constant stoichiometry and an algebraic kinetic law (species, compartment IDs,
  and constant parameters as operands; no ``delay``, no ``piecewise`` -- though ``piecewise``
  parses fine through SymPy's ``Piecewise`` if it does appear, it just isn't tested here).
- Multiple compartments with fixed, non-unity sizes -- each species' ODE is the stoichiometric sum
  of the reactions it participates in, divided by *that species'* own compartment size (the
  standard SBML semantics for a concentration-valued species; this is what makes reactions moving
  species between differently-sized compartments, e.g. cytoplasm to nucleus, come out correct).
- Assignment rules (``<assignmentRule>``) defining a non-constant parameter as a formula, which may
  reference the SBML time symbol -- substituted symbolically into every kinetic law that uses that
  parameter before building the ODEs, since it is not a state variable and has no ODE of its own.
- Initial assignments (``<initialAssignment>``) defining a species' initial value as a formula
  (typically depending on other parameters) instead of a literal ``initialConcentration``. Evaluated
  once, substituting each parameter's *nominal* SBML value -- if a parameter used in an initial
  assignment is later estimated by PEtab, this importer's initial condition will not update as that
  parameter varies during estimation. This did not affect either confirmed benchmark (in both, the
  parameters feeding initial assignments are fixed, not estimated) but is a real scope limitation
  for other PEtab problems.

Not supported at all (raises or silently produces a wrong model -- always spot-check a newly
imported model, e.g. against ``simulatedData_*.tsv`` if the benchmark source provides one):
species with ``hasOnlySubstanceUnits="true"`` (amount, not concentration, semantics), rate rules,
algebraic rules, function definitions, events, and multi-level substitution chains where one
assignment rule's target appears inside another assignment rule's formula (only a single
substitution pass is applied).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import sympy as sp
from sympy.parsing.sympy_parser import convert_xor, parse_expr, standard_transformations

_TRANSFORMATIONS = standard_transformations + (convert_xor,)


@dataclass
class SBMLODESystem:
    """Explicit ODE system extracted from an SBML reaction network.

    Attributes:
        t: Time symbol used in ``odes`` (and, if the model has a time-dependent assignment rule,
            already substituted in).
        state_vars: Species symbols, in a stable order matching ``state_names``/``odes``.
        state_names: Species IDs, same order as ``state_vars``.
        odes: ``d(state)/dt`` SymPy expression for each entry of ``state_vars``, same order.
        params: Symbols for every ``constant="true"`` SBML parameter (i.e. every parameter that
            is not itself defined by an assignment rule -- those are substituted directly into
            ``odes`` and never appear as free symbols).
        param_names: Parameter IDs, same order as ``params``.
        param_values: SBML nominal value for every entry in ``param_names``.
        initial_conditions: Numeric initial value for every entry in ``state_names`` -- from the
            species' literal ``initialConcentration``/``initialAmount``, or (if present) its
            ``<initialAssignment>`` formula evaluated at nominal parameter values.
    """

    t: sp.Symbol
    state_vars: list[sp.Symbol]
    state_names: list[str]
    odes: list[sp.Expr]
    params: list[sp.Symbol]
    param_names: list[str]
    param_values: dict[str, float] = field(default_factory=dict)
    initial_conditions: dict[str, float] = field(default_factory=dict)


def parse_sbml(path: str) -> SBMLODESystem:
    """Parse an SBML reaction-network model into an explicit :class:`SBMLODESystem`.

    Args:
        path: Path to an SBML file.

    Returns:
        The extracted ODE system.

    Raises:
        ImportError: If the ``petab`` extra (``python-libsbml``) is not installed.
        ValueError: If the file fails to parse, or uses a feature outside this module's scope
            (rate rules, algebraic rules -- see the module docstring).
    """
    try:
        import libsbml
    except ImportError as exc:  # pragma: no cover - exercised only without the petab extra
        raise ImportError(
            "parse_sbml requires the 'petab' extra: pip install gsua-csb[petab]"
        ) from exc

    doc = libsbml.readSBMLFromFile(str(path))
    if doc.getNumErrors(libsbml.LIBSBML_SEV_ERROR) > 0:
        msgs = [doc.getError(i).getMessage() for i in range(doc.getNumErrors())]
        raise ValueError(f"Failed to parse SBML file {path!r}: {'; '.join(msgs)}")
    model = doc.getModel()
    if model is None:
        raise ValueError(f"No <model> element found in {path!r}")

    t = sp.Symbol("t")
    compartment_sizes = {c.getId(): float(c.getSize()) for c in model.getListOfCompartments()}

    species_compartment: dict[str, str] = {}
    species_initial: dict[str, float] = {}
    for s in model.getListOfSpecies():
        species_compartment[s.getId()] = s.getCompartment()
        if s.isSetInitialConcentration():
            species_initial[s.getId()] = float(s.getInitialConcentration())
        elif s.isSetInitialAmount():
            species_initial[s.getId()] = float(s.getInitialAmount())
        else:
            species_initial[s.getId()] = 0.0

    param_values: dict[str, float] = {}
    const_param_ids: list[str] = []
    for p in model.getListOfParameters():
        param_values[p.getId()] = float(p.getValue())
        if p.getConstant():
            const_param_ids.append(p.getId())

    for r in model.getListOfRules():
        if not r.isAssignment():
            raise ValueError(
                f"Unsupported SBML rule type for variable {r.getVariable()!r} in {path!r} "
                "(only assignment rules are supported -- see the parse_sbml docstring)"
            )

    # Symbol table used to parse every formula (kinetic laws, assignment rules, initial
    # assignments): species, all parameters (including non-constant ones, substituted away
    # below), compartment IDs (substituted to numeric sizes below), and SBML's time symbol.
    all_ids = list(species_compartment) + list(param_values) + list(compartment_sizes)
    symtab = {name: sp.Symbol(name) for name in all_ids}
    symtab["time"] = t

    def parse_formula(math_ast) -> sp.Expr:
        formula = libsbml.formulaToL3String(math_ast)
        return parse_expr(formula, local_dict=symtab, transformations=_TRANSFORMATIONS)

    assignment_rules = {r.getVariable(): parse_formula(r.getMath()) for r in model.getListOfRules()}
    initial_assignments = {
        ia.getSymbol(): parse_formula(ia.getMath()) for ia in model.getListOfInitialAssignments()
    }

    def substitute_known(expr: sp.Expr) -> sp.Expr:
        # Assignment-rule-defined parameters are not free symbols -- substitute their formula in
        # (one pass; chained assignment rules referencing each other are out of scope).
        expr = expr.subs({symtab[var]: rule for var, rule in assignment_rules.items()})
        # Compartment IDs are always just numeric constants once resolved.
        expr = expr.subs({symtab[cid]: size for cid, size in compartment_sizes.items()})
        return expr

    raw_rhs = {name: sp.Integer(0) for name in species_compartment}
    for rxn in model.getListOfReactions():
        kl = rxn.getKineticLaw()
        if kl is None:
            continue
        rate = substitute_known(parse_formula(kl.getMath()))
        for sr in rxn.getListOfReactants():
            sid = sr.getSpecies()
            if sid in raw_rhs:
                raw_rhs[sid] -= sr.getStoichiometry() * rate
        for sr in rxn.getListOfProducts():
            sid = sr.getSpecies()
            if sid in raw_rhs:
                raw_rhs[sid] += sr.getStoichiometry() * rate

    state_names = list(species_compartment)
    state_vars = [symtab[name] for name in state_names]
    # sp.expand (not sp.simplify -- much cheaper, and simplify() was observed not to fully cancel
    # exp(a)*exp(-a) pairs from assignment-rule substitution anyway; numeric evaluation via
    # lambdify is correct either way, this is purely about keeping parsing fast) tidies up the
    # sum-of-reaction-rates form without the combinatorial cost of full simplification, which was
    # slow enough on Boehm-sized models (8 states, exp() terms) to noticeably slow down every call.
    odes = [sp.expand(raw_rhs[name] / compartment_sizes[species_compartment[name]]) for name in state_names]

    initial_conditions: dict[str, float] = {}
    for name in state_names:
        if name in initial_assignments:
            value = substitute_known(initial_assignments[name]).subs(param_values)
            initial_conditions[name] = float(value)
        else:
            initial_conditions[name] = species_initial[name]

    param_names = [pid for pid in const_param_ids if pid not in assignment_rules]
    params = [symtab[name] for name in param_names]

    return SBMLODESystem(
        t=t,
        state_vars=state_vars,
        state_names=state_names,
        odes=odes,
        params=params,
        param_names=param_names,
        param_values={name: param_values[name] for name in param_names},
        initial_conditions=initial_conditions,
    )
