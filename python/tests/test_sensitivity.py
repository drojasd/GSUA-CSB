import numpy as np
import pytest

from gsua_csb import UserFunctionModel, design_matrix, sensitivity_analysis

# Linear, no-interaction model: y = 2*p0 + 3*p1, p0,p1 ~ U[0,1].
# Analytic first-order Sobol indices (Var(a*p)=a^2/12 for U[0,1]):
#   Var(2 p0) = 4/12,  Var(3 p1) = 9/12,  Total = 13/12
#   Si0 = (4/12)/(13/12) = 4/13 ~= 0.3077
#   Si1 = (9/12)/(13/12) = 9/13 ~= 0.6923
# No interaction terms, so total-order indices should match first-order closely.
TRUE_SI = np.array([4 / 13, 9 / 13])


def _linear_model(extra_fixed=False):
    names = ["p0", "p1"]
    rng = np.array([[0.0, 1.0], [0.0, 1.0]])
    if extra_fixed:
        names = ["p0", "p1", "p2"]
        rng = np.array([[0.0, 1.0], [0.0, 1.0], [5.0, 5.0]])

    def func(p):
        y = 2.0 * p[0] + 3.0 * p[1]
        return np.array([y])

    return UserFunctionModel(func=func, names=names, range=rng)


@pytest.mark.parametrize("method", ["sobol", "jansen", "saltelli", "xiao"])
def test_variance_based_methods_recover_known_linear_indices(method):
    # Si/STi are estimated on the scalar SSE cost J=(y-y_nom)^2, a quadratic transform of the
    # underlying linear model -- it does not reproduce the classic Sobol fractions of Var(y)
    # exactly, but larger-coefficient parameters must still dominate (checked below for every
    # method). Si_vec/STi_vec for 'sobol'/'jansen' are estimated directly on the raw output Y
    # (the same quantity the classic formulas decompose), so those are checked precisely.
    model = _linear_model()
    M = design_matrix(model, 4000, method="latin_hypercube", seed=0)
    result = sensitivity_analysis(model, M, method=method)

    assert result.Si.shape == (2,)
    assert result.STi.shape == (2,)
    # p1 has the larger coefficient -> it should be identified as more influential everywhere.
    assert result.Si[1] > result.Si[0]
    assert result.STi[1] > result.STi[0]

    if method in ("sobol", "jansen"):
        np.testing.assert_allclose(result.Si_vec[:, 0], TRUE_SI, atol=0.12)
        np.testing.assert_allclose(result.STi_vec[:, 0], TRUE_SI, atol=0.12)


def test_fixed_parameter_gets_zero_indices():
    model = _linear_model(extra_fixed=True)
    M = design_matrix(model, 2000, method="latin_hypercube", seed=1)
    result = sensitivity_analysis(model, M, method="xiao")
    assert result.Si[2] == 0.0
    assert result.STi[2] == 0.0
    assert np.all(result.Si_vec[2] == 0.0)


def test_oat_orders_parameters_by_influence():
    model = _linear_model()
    M = design_matrix(model, 500, method="latin_hypercube", seed=2)
    result = sensitivity_analysis(model, M, method="oat")
    assert result.Si_vec is None
    assert result.STi_vec is None
    assert result.Si[1] > result.Si[0]
    np.testing.assert_allclose(np.sum(result.Si), 1.0)


def test_time_dependent_indices_shape_matches_domain():
    names = ["k"]
    rng = np.array([[0.1, 1.0]])
    domain = np.linspace(0, 5, 10)

    def func(p, d):
        return np.exp(-p[0] * d)

    model = UserFunctionModel(func=func, names=names, range=rng, domain=domain)
    M = design_matrix(model, 200, method="latin_hypercube", seed=3)
    result = sensitivity_analysis(model, M, method="jansen")
    assert result.Si_vec.shape == (1, 10)
    assert result.STi_vec.shape == (1, 10)


def test_odd_n_is_truncated_to_even():
    model = _linear_model()
    M = design_matrix(model, 201, method="latin_hypercube", seed=4)
    result = sensitivity_analysis(model, M, method="xiao")
    assert result.Y.shape[0] == 200


def test_unknown_method_raises():
    model = _linear_model()
    M = design_matrix(model, 100, seed=0)
    with pytest.raises(ValueError, match="Unknown sensitivity method"):
        sensitivity_analysis(model, M, method="bogus")


def test_pod_out_of_range_raises():
    model = _linear_model()
    M = design_matrix(model, 100, seed=0)
    with pytest.raises(ValueError, match="pod must be"):
        sensitivity_analysis(model, M, method="xiao", pod=3.0)


def test_swapped_model_and_matrix_raises_clear_error():
    # gsua_sa(M, T) in MATLAB but sensitivity_analysis(model, M) here; a swap must be caught.
    import numpy as np
    import pytest
    from gsua_csb import (UserFunctionModel, design_matrix, sensitivity_analysis,
                          uncertainty_analysis)
    m = UserFunctionModel(func=lambda p, t: (p[0] * np.exp(-p[1] * t))[None, :],
                          names=["a", "r"], range=np.array([[1.0, 3.0], [0.1, 1.0]]),
                          domain=np.linspace(0, 5, 20), output_names=["o"])
    M = design_matrix(m, 20, seed=0)
    with pytest.raises(TypeError, match="reversed"):
        sensitivity_analysis(M, m)
    with pytest.raises(TypeError, match="reversed"):
        uncertainty_analysis(M, m)


def test_n_jobs_matches_serial():
    import numpy as np
    from gsua_csb import UserFunctionModel, design_matrix, uncertainty_analysis, sensitivity_analysis
    m = UserFunctionModel(func=lambda p, t: (p[0] * np.exp(-p[1] * t))[None, :],
                          names=["a", "r"], range=np.array([[1.0, 3.0], [0.1, 1.0]]),
                          domain=np.linspace(0, 5, 20), output_names=["o"])
    M = design_matrix(m, 24, seed=0)
    np.testing.assert_allclose(uncertainty_analysis(m, M, n_jobs=1).Y,
                               uncertainty_analysis(m, M, n_jobs=2).Y)
    s1 = sensitivity_analysis(m, M, method="sobol", n_jobs=1)
    s2 = sensitivity_analysis(m, M, method="sobol", n_jobs=2)
    np.testing.assert_allclose(s1.STi, s2.STi)
