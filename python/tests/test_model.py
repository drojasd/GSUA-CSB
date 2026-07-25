import numpy as np
import pytest

from gsua_csb import UserFunctionModel, build_range


class TestBuildRange:
    def test_range_method_passthrough(self):
        bounds = np.array([[0.0, 1.0], [2.0, 4.0]])
        r = build_range(bounds, None, "range")
        np.testing.assert_array_equal(r, bounds)

    def test_percent_method(self):
        r = build_range([10.0, 20.0], [10.0, 50.0], "percent")
        np.testing.assert_allclose(r, [[9.0, 11.0], [10.0, 30.0]])

    def test_std_method(self):
        r = build_range([5.0], [1.0], "std")
        np.testing.assert_allclose(r, [[4.0, 6.0]])

    def test_range_method_requires_2d_shape(self):
        with pytest.raises(ValueError, match="range"):
            build_range([1.0, 2.0], None, "range")

    def test_missing_spread_raises(self):
        with pytest.raises(ValueError, match="spread"):
            build_range([1.0], None, "percent")


class TestUserFunctionModel:
    def test_evaluate_domain_less(self):
        model = UserFunctionModel(
            func=lambda p: np.array([p[0] + p[1]]),
            names=["a", "b"],
            range=np.array([[0.0, 1.0], [0.0, 1.0]]),
        )
        out = model.evaluate(np.array([0.3, 0.4]))
        np.testing.assert_allclose(out, [0.7])

    def test_evaluate_with_domain(self):
        t = np.linspace(0, 1, 5)
        model = UserFunctionModel(
            func=lambda p, d: p[0] * d,
            names=["k"],
            range=np.array([[0.0, 2.0]]),
            domain=t,
        )
        out = model.evaluate(np.array([2.0]))
        np.testing.assert_allclose(out, 2.0 * t)

    def test_evaluate_batch_default_loops(self):
        t = np.linspace(0, 1, 4)
        model = UserFunctionModel(
            func=lambda p, d: p[0] * d,
            names=["k"],
            range=np.array([[0.0, 2.0]]),
            domain=t,
        )
        params = np.array([[1.0], [2.0], [3.0]])
        out = model.evaluate_batch(params)
        assert out.shape == (3, 4)
        np.testing.assert_allclose(out[1], 2.0 * t)

    def test_evaluate_batch_vectorized_calls_func_once(self):
        calls = []

        def batched_func(p, d):
            calls.append(p.shape)
            return p[:, [0]] * d[None, :]

        t = np.linspace(0, 1, 4)
        model = UserFunctionModel(
            func=batched_func,
            names=["k"],
            range=np.array([[0.0, 2.0]]),
            domain=t,
            vectorized=True,
        )
        params = np.array([[1.0], [2.0]])
        out = model.evaluate_batch(params)
        assert len(calls) == 1  # one batched call, not a Python loop
        assert out.shape == (2, 4)

    def test_nominal_defaults_to_midpoint(self):
        model = UserFunctionModel(
            func=lambda p: p,
            names=["a"],
            range=np.array([[0.0, 10.0]]),
        )
        np.testing.assert_allclose(model.nominal, [5.0])

    def test_fixed_mask_from_collapsed_range(self):
        model = UserFunctionModel(
            func=lambda p: p,
            names=["a", "b"],
            range=np.array([[1.0, 1.0], [0.0, 2.0]]),
        )
        np.testing.assert_array_equal(model.fixed, [True, False])

    def test_fix_collapses_range_and_updates_nominal(self):
        model = UserFunctionModel(
            func=lambda p: p,
            names=["a", "b"],
            range=np.array([[0.0, 10.0], [0.0, 10.0]]),
        )
        model.fix("a")
        assert model.fixed[0]
        assert model.range[0, 0] == model.range[0, 1] == model.nominal[0]

    def test_fix_by_index_with_explicit_value(self):
        model = UserFunctionModel(
            func=lambda p: p,
            names=["a"],
            range=np.array([[0.0, 10.0]]),
        )
        model.fix(0, value=7.0)
        np.testing.assert_allclose(model.range[0], [7.0, 7.0])
        assert model.nominal[0] == 7.0

    def test_from_bounds_generates_default_names(self):
        model = UserFunctionModel.from_bounds(
            func=lambda p: p,
            values=np.array([[0.0, 1.0], [0.0, 1.0]]),
        )
        assert model.names == ["0", "1"]
