import numpy as np
import pytest

from gsua_csb import UserFunctionModel, design_matrix


def _model(fixed_last=False):
    rng = np.array([[0.0, 10.0], [-5.0, 5.0], [100.0, 200.0]])
    if fixed_last:
        rng[2, :] = 150.0
    return UserFunctionModel(
        func=lambda p: p,
        names=["a", "b", "c"],
        range=rng,
    )


@pytest.mark.parametrize("method", ["latin_hypercube", "uniform", "sobol"])
def test_shape_and_bounds(method):
    model = _model()
    M = design_matrix(model, 16, method=method, seed=0)
    assert M.shape == (16, 3)
    assert np.all(M >= model.range[:, 0])
    assert np.all(M <= model.range[:, 1])


@pytest.mark.parametrize("method", ["latin_hypercube", "uniform", "sobol"])
def test_fixed_factor_held_constant(method):
    model = _model(fixed_last=True)
    M = design_matrix(model, 12, method=method, seed=0)
    assert np.all(M[:, 2] == 150.0)
    # Free factors still vary across rows.
    assert np.ptp(M[:, 0]) > 0
    assert np.ptp(M[:, 1]) > 0


def test_all_fixed_returns_constant_rows():
    model = _model()
    model.range = np.tile([[3.0, 3.0]], (3, 1))
    M = design_matrix(model, 5)
    np.testing.assert_array_equal(M, np.tile([3.0, 3.0, 3.0], (5, 1)))


def test_seed_reproducibility():
    model = _model()
    M1 = design_matrix(model, 20, method="latin_hypercube", seed=42)
    M2 = design_matrix(model, 20, method="latin_hypercube", seed=42)
    np.testing.assert_array_equal(M1, M2)


def test_unknown_method_raises():
    model = _model()
    with pytest.raises(ValueError, match="Unknown design method"):
        design_matrix(model, 5, method="bogus")


def test_latin_hypercube_stratified_per_dimension():
    # LHS guarantees exactly one sample per equal-width stratum, per dimension.
    model = _model()
    n = 20
    M = design_matrix(model, n, method="latin_hypercube", seed=1)
    for col, (lo, hi) in zip(M.T, model.range):
        strata = np.floor((col - lo) / (hi - lo) * n).astype(int)
        strata = np.clip(strata, 0, n - 1)
        assert sorted(strata) == list(range(n))
