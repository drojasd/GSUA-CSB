"""A scoped SBML reaction-network parser producing explicit SymPy ODEs.

There is no MATLAB equivalent to port here -- this exists purely to make :func:`gsua_csb.
load_petab` possible, since PEtab problems define their dynamics as SBML reaction networks
(species + reactions + stoichiometry + kinetic laws), not the explicit ``dx/dt = ...`` SymPy
expressions :class:`gsua_csb.SymbolicODEModel` expects.

This is **not** a general SBML importer. Full SBML support (events, delays, algebraic rules,
spatial/hierarchical compartments, unit conversion, ``hasOnlySubstanceUnits`` species, ...) is a
large, separate undertaking better served by dedicated tools (libSBML itself, AMICI, COPASI,
tellurium). This module supports exactly the subset needed to import the reaction-network SBML
models used by the PEtab benchmark collection's simpler problems (confirmed against
``Perelson_Science1996``, ``Boehm_JProteomeRes2014``, ``Giordano_Nature2020``,
``Okuonghae_ChaosSolitonsFractals2020``, ``Bertozzi_PNAS2020``):

- Reactions with constant stoichiometry and an algebraic kinetic law (species, compartment IDs,
  and constant parameters as operands; no ``delay``). ``piecewise(...)`` (including its use in a
  time-dependent assignment rule, e.g. a rate parameter that steps to a new value at known
  intervention dates -- common in epidemiology models) and ``&&``/``||`` are supported.
- Multiple compartments with fixed, non-unity sizes -- each species' ODE is the stoichiometric sum
  of the reactions it participates in, divided by *that species'* own compartment size (the
  standard SBML semantics for a concentration-valued species; this is what makes reactions moving
  species between differently-sized compartments, e.g. cytoplasm to nucleus, come out correct).
- ``<functionDefinition>`` (a named, reusable rate-law formula called from one or more kinetic
  laws) is resolved by inlining: its body is substituted at every call site, with its bound
  variable names (which, for epidemiology models, are routinely single letters like ``S``/``I``/
  ``D`` that collide with SymPy globals such as the imaginary unit -- handled explicitly) scoped to
  that call.
- Assignment rules (``<assignmentRule>``) defining a non-constant parameter *or species* (a
  "reporter" species with no ODE of its own, e.g. a derived aggregate like "current active cases"
  used only by an observable) as a formula, which may reference the SBML time symbol --
  substituted symbolically into every formula that uses it (including another assignment rule's
  formula, resolved iteratively) before building the ODEs. A species that is itself an
  assignment-rule target is excluded from the returned ``state_vars``/``odes`` entirely (it isn't
  an integrated state), and every assignment-rule target -- species or parameter -- is exposed,
  fully resolved, via ``derived`` for a caller (e.g. :func:`gsua_csb.load_petab`, for a PEtab
  ``observableFormula`` naming one directly) to substitute into its own formulas.
- A constant parameter's value can likewise be given by an ``<initialAssignment>`` formula instead
  of a literal ``<parameter value="...">`` (e.g. ``beta = R0*gamma/N``, a "derived constant" built
  from directly-meaningful rate parameters) -- resolved once (iteratively) per :func:`parse_sbml`
  call using ``overrides``, then substituted away the same way; it never appears in ``param_names``.
- Initial assignments (``<initialAssignment>``) defining a species' initial value as a formula
  (typically depending on other parameters) instead of a literal ``initialConcentration``. Evaluated
  once per :func:`parse_sbml` call, substituting each parameter's *effective* (nominal, or
  ``overrides``-supplied) value -- if a parameter used in an initial assignment is estimated by
  PEtab (not fixed, not condition-overridden), the initial condition uses its PEtab *nominal* value
  and will not update as that parameter varies during estimation; a real scope limitation, but a
  narrow one (only reachable when an estimated parameter also happens to feed a species' initial
  value), not something either confirmed benchmark exercises.

Not supported at all (raises or silently produces a wrong model -- always spot-check a newly
imported model, e.g. against ``simulatedData_*.tsv`` if the benchmark source provides one):
species with ``hasOnlySubstanceUnits="true"`` (amount, not concentration, semantics), rate rules,
algebraic rules, events, and a function definition or assignment rule whose formula contains
anything beyond arithmetic/``piecewise``/comparison/logical operators and other resolvable calls
(e.g. a MathML function libSBML's L3 printer renders as a bare, unhandled name).
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
        state_vars: Species symbols, in a stable order matching ``state_names``/``odes``. Excludes
            any species that is itself an assignment-rule target -- see ``derived``.
        state_names: Species IDs, same order as ``state_vars``.
        odes: ``d(state)/dt`` SymPy expression for each entry of ``state_vars``, same order.
        params: Symbols for every ``constant="true"`` SBML parameter that is not itself defined by
            an assignment rule or an initialAssignment (those are substituted directly into
            ``odes`` and never appear as free symbols -- see ``derived``).
        param_names: Parameter IDs, same order as ``params``.
        param_values: SBML nominal (or ``overrides``-supplied) value for every entry in
            ``param_names``.
        initial_conditions: Numeric initial value for every entry in ``state_names`` -- from the
            species' literal ``initialConcentration``/``initialAmount``, its ``<initialAssignment>``
            formula evaluated at effective parameter values, or a species-id entry in ``overrides``.
        derived: Every assignment-rule target (species or parameter id), mapped to its formula
            fully resolved (rule chains included) down to ``state_vars``/``params``/``t`` -- for a
            caller whose own formula (e.g. a PEtab ``observableFormula``) references one of these
            "reporter" quantities directly, without duplicating this module's substitution logic.
    """

    t: sp.Symbol
    state_vars: list[sp.Symbol]
    state_names: list[str]
    odes: list[sp.Expr]
    params: list[sp.Symbol]
    param_names: list[str]
    param_values: dict[str, float] = field(default_factory=dict)
    initial_conditions: dict[str, float] = field(default_factory=dict)
    derived: dict[str, sp.Expr] = field(default_factory=dict)


def parse_sbml(path: str, overrides: dict[str, float] | None = None) -> SBMLODESystem:
    """Parse an SBML reaction-network model into an explicit :class:`SBMLODESystem`.

    Args:
        path: Path to an SBML file.
        overrides: Optional ``{species or parameter id: value}`` map, used by
            :func:`gsua_csb.load_petab` to apply a PEtab condition table's per-condition
            overrides. A species id in ``overrides`` fixes that species' initial condition
            directly, overriding any ``<initialAssignment>`` or literal SBML value. A parameter
            id in ``overrides`` is substituted for that parameter's SBML nominal value when
            evaluating any ``<initialAssignment>`` formula that depends on it (this is the only
            place a constant parameter's value is baked in at parse time -- ``odes`` itself stays
            symbolic in every parameter regardless of ``overrides``), and is reflected in the
            returned ``param_values``.

    Returns:
        The extracted ODE system.

    Raises:
        ImportError: If the ``petab`` extra (``python-libsbml``) is not installed.
        ValueError: If the file fails to parse, or uses a feature outside this module's scope
            (rate rules, algebraic rules -- see the module docstring).
    """
    overrides = overrides or {}
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
    # libsbml.formulaToL3String renders a degenerate (single- or zero-operand) MathML <times>/
    # <plus> apply as literal L3 function-call syntax (e.g. "times(x)" for a one-operand product)
    # rather than infix -- these two names are reserved in the L3 grammar for exactly that case, so
    # provide them as ordinary Python callables sympy's parser can invoke directly.
    symtab.setdefault("times", lambda *args: sp.Mul(*args) if args else sp.Integer(1))
    symtab.setdefault("plus", lambda *args: sp.Add(*args) if args else sp.Integer(0))

    def _piecewise(*args):
        # SBML L3 piecewise(v1, c1, v2, c2, ..., [default]) -- an even arg count is (value,
        # condition) pairs with no fallback; an odd count has a trailing unconditional default.
        pairs = [(args[i], args[i + 1]) for i in range(0, len(args) - 1, 2)]
        if len(args) % 2:
            pairs.append((args[-1], True))
        return sp.Piecewise(*pairs)

    symtab.setdefault("piecewise", _piecewise)

    def parse_formula(math_ast, extra_locals: dict | None = None) -> sp.Expr:
        # && / || are libsbml's L3 spelling of logical and/or -- not valid Python syntax, so
        # rewrite to the bitwise operators sympy overloads for symbolic boolean combination
        # (Python's own `and`/`or` would instead try immediate truthiness on a symbolic relational
        # and raise). Every operand libsbml emits is already fully parenthesized, so this is a safe
        # textual substitution regardless of Python's tighter-than-comparison `&`/`|` precedence.
        formula = libsbml.formulaToL3String(math_ast).replace("&&", "&").replace("||", "|")
        local_dict = {**symtab, **extra_locals} if extra_locals else symtab
        return parse_expr(formula, local_dict=local_dict, transformations=_TRANSFORMATIONS)

    # <functionDefinition> (a named lambda, e.g. a shared rate-law helper called from several
    # kinetic laws) is resolved by inlining: parse its body once against the bound variable names
    # as ordinary symbols, then register a Python callable under its id that substitutes the actual
    # call-site arguments for those symbols. Must happen before any formula that might call it.
    # extra_locals is essential here, not an optimization: SIR-style epidemiology models routinely
    # name bound variables S/I/D/E/N/O/Q/C, which are also sympy globals (S the singleton registry,
    # I the imaginary unit, ...) -- without shadowing them for this parse, "S * (...)" would try to
    # multiply sympy.S itself rather than treat S as this function's first argument.
    for fd in model.getListOfFunctionDefinitions():
        arg_names = [fd.getArgument(i).getName() for i in range(fd.getNumArguments())]
        arg_syms = [sp.Symbol(a) for a in arg_names]
        body_expr = parse_formula(fd.getBody(), extra_locals=dict(zip(arg_names, arg_syms)))

        def _make_call(arg_syms=arg_syms, body_expr=body_expr):
            def _call(*call_args):
                return body_expr.subs(dict(zip(arg_syms, call_args)), simultaneous=True)
            return _call

        symtab[fd.getId()] = _make_call()

    assignment_rules = {r.getVariable(): parse_formula(r.getMath()) for r in model.getListOfRules()}
    initial_assignments = {
        ia.getSymbol(): parse_formula(ia.getMath()) for ia in model.getListOfInitialAssignments()
    }

    # Resolved once, outside substitute_known's definition so the closure below picks up the
    # final values (populated further down, before it is first called).
    resolved_derived: dict[str, float] = {}

    def substitute_known(expr: sp.Expr) -> sp.Expr:
        # Assignment-rule-defined species/parameters are not free symbols -- substitute their
        # formula in, repeating until nothing changes (bounded by the rule count) so a rule whose
        # own formula references another rule's target (e.g. Giordano's CumulativeInfected, defined
        # in terms of CurrentTotalInfected, itself an assignment rule) fully resolves too.
        rule_subs = {symtab[var]: rule for var, rule in assignment_rules.items()}
        for _ in range(len(rule_subs) + 1):
            new_expr = expr.subs(rule_subs)
            if new_expr == expr:
                break
            expr = new_expr
        # Compartment IDs are always just numeric constants once resolved.
        expr = expr.subs({symtab[cid]: size for cid, size in compartment_sizes.items()})
        # Constant parameters whose value is itself an <initialAssignment> formula (see below).
        if resolved_derived:
            expr = expr.subs(resolved_derived)
        return expr

    # A constant SBML parameter's value can be given by an <initialAssignment> formula instead of
    # (or in addition to, in which case the formula wins) a literal <parameter value="...">
    # attribute -- a "derived constant" (e.g. beta = R0*gamma/N), common in epidemiology SBML
    # models built from a handful of directly-meaningful rate parameters. Resolve these to numbers
    # (iteratively, since one derived constant's formula may reference another) using this
    # condition's overrides, then fold them into eff_param_values and substitute them away in odes/
    # initial assignments exactly like an assignment-rule-defined parameter -- they are not a free
    # degree of freedom, so they never appear in the returned param_names/params.
    eff_param_values = dict(param_values)
    eff_param_values.update({k: v for k, v in overrides.items() if k in param_values})
    remaining = {
        name: expr for name, expr in initial_assignments.items()
        if name in param_values and name not in species_compartment
    }
    for _ in range(len(remaining) + 1):
        if not remaining:
            break
        progressed = False
        for name, expr in list(remaining.items()):
            value = substitute_known(expr).subs(eff_param_values)
            if not value.free_symbols:
                fval = float(value)
                resolved_derived[name] = fval
                eff_param_values[name] = fval
                del remaining[name]
                progressed = True
        if not progressed:
            break
    if remaining:
        raise ValueError(
            f"Could not resolve initialAssignment-derived parameter(s) {sorted(remaining)} in "
            f"{path!r} (unsupported dependency structure -- see the parse_sbml docstring)"
        )

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

    # A species that is itself an assignment-rule target (a "reporter" species -- e.g. Giordano's
    # CurrentDiagnosedInfected = Diagnosed + Recognized + Threatened, used only to name a derived
    # aggregate for an observable) is not an ODE state: its value is the rule, not an integral of
    # reaction rates, even if it happens to also carry an <initialConcentration>/<initialAmount>.
    # Treating it as one would integrate a constant zero derivative and silently freeze it at its
    # initial value instead of tracking the real states it's defined from.
    state_names = [name for name in species_compartment if name not in assignment_rules]
    state_vars = [symtab[name] for name in state_names]
    # sp.expand (not sp.simplify -- much cheaper, and simplify() was observed not to fully cancel
    # exp(a)*exp(-a) pairs from assignment-rule substitution anyway; numeric evaluation via
    # lambdify is correct either way, this is purely about keeping parsing fast) tidies up the
    # sum-of-reaction-rates form without the combinatorial cost of full simplification, which was
    # slow enough on Boehm-sized models (8 states, exp() terms) to noticeably slow down every call.
    odes = [sp.expand(raw_rhs[name] / compartment_sizes[species_compartment[name]]) for name in state_names]

    initial_conditions: dict[str, float] = {}
    for name in state_names:
        if name in overrides:
            initial_conditions[name] = overrides[name]
        elif name in initial_assignments:
            value = substitute_known(initial_assignments[name]).subs(eff_param_values)
            initial_conditions[name] = float(value)
        else:
            initial_conditions[name] = species_initial[name]

    param_names = [pid for pid in const_param_ids if pid not in assignment_rules and pid not in resolved_derived]
    params = [symtab[name] for name in param_names]

    # Every assignment-rule target (species or parameter), fully resolved (chains included) down
    # to real state variables and constant parameters -- so callers parsing an unrelated formula
    # that references one directly (e.g. a PEtab observableFormula naming a "reporter" species like
    # CurrentDiagnosedInfected above) can substitute it in without duplicating substitute_known.
    derived = {name: substitute_known(symtab[name]) for name in assignment_rules}

    return SBMLODESystem(
        t=t,
        state_vars=state_vars,
        state_names=state_names,
        odes=odes,
        params=params,
        param_names=param_names,
        param_values={name: eff_param_values[name] for name in param_names},
        initial_conditions=initial_conditions,
        derived=derived,
    )
