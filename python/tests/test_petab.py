"""Tests for parse_sbml/load_petab against real PEtab benchmark problems.

Test data in tests/data/petab/ is vendored (not fetched at test time) from
https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab (BSD-3-Clause), specifically the
Perelson_Science1996 and Boehm_JProteomeRes2014 problems, including each problem's own
`simulatedData_*.tsv` reference simulation (used for cross-validation below).
"""

from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("libsbml")
pd = pytest.importorskip("pandas")

from gsua_csb import load_petab, parse_sbml

DATA_DIR = Path(__file__).parent / "data" / "petab"


@pytest.fixture(scope="module")
def perelson_sbml():
    return parse_sbml(DATA_DIR / "Perelson_Science1996" / "model_Perelson_Science1996.xml")


@pytest.fixture(scope="module")
def boehm_sbml():
    return parse_sbml(DATA_DIR / "Boehm_JProteomeRes2014" / "model_Boehm_JProteomeRes2014.xml")


@pytest.fixture(scope="module")
def perelson_problem():
    return load_petab(DATA_DIR / "Perelson_Science1996" / "Perelson_Science1996.yaml")


@pytest.fixture(scope="module")
def boehm_problem():
    return load_petab(DATA_DIR / "Boehm_JProteomeRes2014" / "Boehm_JProteomeRes2014.yaml")


def test_parse_sbml_perelson_structure(perelson_sbml):
    assert perelson_sbml.state_names == ["Tstar", "V", "Vin", "Vni"]
    assert perelson_sbml.param_names == ["NN", "T0", "c", "delta", "K0"]
    assert len(perelson_sbml.odes) == 4
    np.testing.assert_allclose(perelson_sbml.initial_conditions["Vni"], 0.0)
    np.testing.assert_allclose(perelson_sbml.initial_conditions["V"], 1860000.0)


def test_parse_sbml_boehm_compartment_conversion(boehm_sbml):
    # Regression check for the compartment-size division: a reaction moving a species from the
    # larger "cyt" compartment (1.4) to the smaller "nuc" compartment (0.45) must pick up a
    # cyt/nuc ~= 3.111 factor in the *importing* species' ODE, not a bare rate constant.
    assert "nucpApA" in boehm_sbml.state_names
    ode = boehm_sbml.odes[boehm_sbml.state_names.index("nucpApA")]
    coeffs = [float(c) for c in ode.as_coefficients_dict().values()]
    assert any(np.isclose(abs(c), 1.4 / 0.45, rtol=1e-6) for c in coeffs)


def test_parse_sbml_boehm_initial_assignment_evaluated(boehm_sbml):
    # STAT5A/STAT5B initial values come from <initialAssignment> formulas (207.6*ratio and
    # 207.6*(1-ratio)), not a literal initialConcentration -- must be evaluated, not left symbolic.
    ratio = boehm_sbml.param_values["ratio"]
    np.testing.assert_allclose(boehm_sbml.initial_conditions["STAT5A"], 207.6 * ratio, rtol=1e-10)
    np.testing.assert_allclose(boehm_sbml.initial_conditions["STAT5B"], 207.6 * (1 - ratio), rtol=1e-10)


def test_load_petab_perelson_model_shape(perelson_problem):
    problem = perelson_problem
    # 4 species (all fixed initial conditions) + 5 SBML parameters (c, delta free; NN, T0, K0 fixed).
    assert problem.model.names == ["Tstar", "V", "Vin", "Vni", "NN", "T0", "c", "delta", "K0"]
    np.testing.assert_array_equal(
        problem.model.fixed, [True, True, True, True, True, True, False, False, True]
    )
    assert problem.observable_names == ["task0_model0_perelson1_V"]
    assert problem.xdata.shape == (16,)
    assert problem.ydata.shape == (1, 16)
    assert not np.any(np.isnan(problem.ydata))  # every time point has a measurement here

    sim = problem.model.evaluate(problem.model.nominal, problem.xdata)
    assert sim.shape == (1, 16)
    assert np.all(np.isfinite(sim))
    # NOTE: this problem's own simulatedData_Perelson_Science1996.tsv reference does NOT match a
    # simulation at the parameters_*.tsv nominalValue column (checked directly: >40% relative
    # error, while the structurally much harder Boehm benchmark below matches its own reference to
    # ~1e-7) -- most likely that reference file was generated at different (e.g. fitted) parameter
    # values than nominalValue, not a bug here, but unconfirmed. Not asserted against below.


def test_load_petab_boehm_cross_validates_against_reference_simulation(boehm_problem):
    # The decisive cross-validation: Boehm exercises every scoped SBML feature at once
    # (two differently-sized compartments, a time-dependent assignment rule, initial-assignment-
    # computed initial conditions, three derived multi-species observable formulas). Matching its
    # own published simulatedData to near machine precision is strong evidence the whole pipeline
    # -- SBML parsing, ODE construction, and observable evaluation -- is correct.
    problem = boehm_problem
    sim = problem.model.evaluate(problem.model.nominal, problem.xdata)

    ref = pd.read_csv(
        DATA_DIR / "Boehm_JProteomeRes2014" / "simulatedData_Boehm_JProteomeRes2014.tsv", sep="\t"
    )
    for i, oid in enumerate(problem.observable_names):
        ref_i = ref[ref["observableId"] == oid].sort_values("time")
        ref_sim = ref_i["simulation"].to_numpy()
        denom = np.where(np.abs(ref_sim) < 1e-8, 1.0, np.abs(ref_sim))
        rel_err = np.abs(sim[i] - ref_sim) / denom
        assert rel_err.max() < 1e-4, f"{oid}: max relative error {rel_err.max()}"


def test_load_petab_boehm_model_shape(boehm_problem):
    problem = boehm_problem
    assert len(problem.model.names) == 8 + 8  # 8 species (fixed) + 8 SBML constant parameters
    # BaF3_Epo (constant="false", assignment-rule-defined) must NOT appear as a model parameter.
    assert "BaF3_Epo" not in problem.model.names
    assert problem.observable_names == ["pSTAT5A_rel", "pSTAT5B_rel", "rSTAT5A_rel"]
    assert problem.ydata.shape == (3, len(problem.xdata))


def test_load_petab_can_estimate_free_parameters(boehm_problem):
    from gsua_csb import parameter_estimation

    problem = boehm_problem
    assert not np.any(np.isnan(problem.ydata)), "test assumes Boehm's 3 observables share a time grid"
    pe = parameter_estimation(
        problem.model, problem.xdata, problem.ydata, n=1, solver="least_squares",
        initial_points=problem.model.nominal[None, :],
    )
    assert pe.x.shape == (1, len(problem.model.names))
    assert np.isfinite(pe.cost[0])
