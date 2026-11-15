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
    # J is the cost of the base sample, one value per row of Y (as for every method).
    assert result.J.shape == (result.Y.shape[0],)


def test_oat_linear_model_matches_analytic_shares():
    # y = x1 + 2*x2 on [-4,4]^2, nominal 0: OAT cost of factor k is c_k^2 x_k^2, so
    # Si = c^4 / sum(c^4) = [1/17, 16/17] -- same value MATLAB's gsua_sa 'OAT' gives (tests/gsuaSaOatTest.m).
    model = _linear12_model()
    M = design_matrix(model, 4000, seed=21)
    result = sensitivity_analysis(model, M, method="oat")
    np.testing.assert_allclose(result.Si, [1 / 17, 16 / 17], atol=0.01)


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


# ---------------------------------------------------------------------------------------------------
# Shapley effects (method="shapley"), mirroring MATLAB tests/gsuaShapleyTest.m.
# Ishigami (a=7, b=0.1) on [-pi, pi]^3: S = [0.3139 0.4424 0], ST = [0.5576 0.4424 0.2437],
# Shapley Sh = [0.4358 0.4424 0.1218] (the x1-x3 interaction split evenly).
# ---------------------------------------------------------------------------------------------------
SH_ISHIGAMI = np.array([0.4358, 0.4424, 0.1218])
S_ISHIGAMI = np.array([0.3139, 0.4424, 0.0])
ST_ISHIGAMI = np.array([0.5576, 0.4424, 0.2437])


def _ishigami_model(x3_range=(-np.pi, np.pi)):
    def func(p):
        return np.array([np.sin(p[0]) + 7.0 * np.sin(p[1]) ** 2 + 0.1 * p[2] ** 4 * np.sin(p[0])])

    rng = np.array([[-np.pi, np.pi], [-np.pi, np.pi], list(x3_range)])
    return UserFunctionModel(func=func, names=["x1", "x2", "x3"], range=rng)


def _linear12_model():
    # y = x1 + 2*x2. With corr(x1, x2) = 0.5 (unit normals): Sh = [0.3393, 0.6607];
    # with independent identical uniforms: Sh = [0.2, 0.8].
    return UserFunctionModel(func=lambda p: np.array([p[0] + 2.0 * p[1]]), names=["x1", "x2"],
                             range=np.array([[-4.0, 4.0], [-4.0, 4.0]]))


def test_shapley_independent_matches_ishigami():
    model = _ishigami_model()
    M = design_matrix(model, 8000, seed=11)
    r = sensitivity_analysis(model, M, method="shapley", seed=11)
    np.testing.assert_allclose(r.Sh_vec[:, 0], SH_ISHIGAMI, atol=0.04)
    np.testing.assert_allclose(r.Sh_vec.sum(axis=0), 1.0, atol=1e-10)
    np.testing.assert_allclose(r.Sh.sum(), 1.0, atol=1e-10)
    assert np.all(np.isnan(r.Si)) and r.Si_vec is None   # Sobol fields are not computed


def test_shapley_given_matches_ishigami_and_lies_between_first_and_total():
    model = _ishigami_model()
    M = design_matrix(model, 4000, seed=12)
    r = sensitivity_analysis(model, M, method="shapley", shapley_type="given")
    sh = r.Sh_vec[:, 0]
    np.testing.assert_allclose(sh, SH_ISHIGAMI, atol=0.03)
    assert np.all(sh + 0.03 >= S_ISHIGAMI) and np.all(sh - 0.03 <= ST_ISHIGAMI)
    np.testing.assert_allclose(r.Sh.sum(), 1.0, atol=1e-10)


def test_shapley_given_uses_the_dependence_of_a_correlated_pool():
    model = _linear12_model()
    rng = np.random.default_rng(14)
    pool = rng.multivariate_normal([0.0, 0.0], [[1.0, 0.5], [0.5, 1.0]], size=4000)   # (n_pool, Np)
    M = design_matrix(model, 10, seed=0)
    r = sensitivity_analysis(model, M, method="shapley", shapley_type="given", pool=pool)
    np.testing.assert_allclose(r.Sh_vec[:, 0], [0.3393, 0.6607], atol=0.03)
    # Same model, independent factors: the dependence, not only the model, decides the attribution.
    ri = sensitivity_analysis(model, design_matrix(model, 4000, seed=15), method="shapley", seed=15)
    np.testing.assert_allclose(ri.Sh_vec[:, 0], [0.2, 0.8], atol=0.04)


@pytest.mark.parametrize("shapley_type", ["independent", "given"])
def test_shapley_fixed_factor_is_zero(shapley_type):
    model = _ishigami_model(x3_range=(1.0, 1.0))      # Sh = [0.0899, 0.9101, 0]
    M = design_matrix(model, 4000, seed=16)
    r = sensitivity_analysis(model, M, method="shapley", shapley_type=shapley_type, seed=16)
    assert r.Sh[2] == 0.0 and np.all(r.Sh_vec[2] == 0.0)
    np.testing.assert_allclose(r.Sh_vec[:2, 0], [0.0899, 0.9101], atol=0.05)


def test_shapley_given_evaluates_repeated_vectors_once():
    model = _ishigami_model()
    base = design_matrix(model, 300, seed=17)
    pool = np.vstack([base, base[:100]])
    r = sensitivity_analysis(model, base, method="shapley", shapley_type="given", pool=pool)
    assert r.Y.shape[0] == 300


def test_shapley_is_reproducible_with_seed():
    model = _ishigami_model()
    M = design_matrix(model, 300, seed=18)
    a = sensitivity_analysis(model, M, method="shapley", seed=5)
    b = sensitivity_analysis(model, M, method="shapley", seed=5)
    np.testing.assert_array_equal(a.Sh_vec, b.Sh_vec)


def test_shapley_bad_pool_or_type_raises():
    model = _ishigami_model()
    M = design_matrix(model, 50, seed=0)
    with pytest.raises(ValueError, match="pool must be"):
        sensitivity_analysis(model, M, method="shapley", shapley_type="given", pool=np.ones((20, 2)))
    with pytest.raises(ValueError, match="Unknown shapley_type"):
        sensitivity_analysis(model, M, method="shapley", shapley_type="bogus")
