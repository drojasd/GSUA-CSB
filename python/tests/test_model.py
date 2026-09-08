import numpy as np
import pytest

from gsua_csb import UserFunctionModel, build_range
from gsua_csb._model import Model


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


class TestOutputSelection:
    """`Model.output` is the single place output selection lives.

    MATLAB sets it once on the summary table (`CustomProperties.output`) and every downstream
    function follows. Before this existed the Python port had no model-level equivalent:
    `output_index` was a per-call argument on sensitivity/uncertainty analysis only, and
    `parameter_estimation`/`profile_likelihood` could not select an output at all, so a
    multi-output model had to be wrapped in a second model just to be fitted.
    """

    @staticmethod
    def _three_output_model(**kw):
        return UserFunctionModel(
            func=lambda p, t: np.vstack([p[0] * t, p[1] * t**2, p[0] + p[1] + 0 * t]),
            names=["a", "b"],
            range=np.array([[1.0, 2.0], [1.0, 2.0]]),
            domain=np.linspace(0, 1, 5),
            output_names=["lin", "quad", "const"],
            **kw,
        )

    def test_defaults_to_every_output(self):
        m = self._three_output_model()
        assert m.output is None
        assert m.evaluate(np.array([1.0, 2.0])).shape == (3, 5)
        assert m.active_output_names == ["lin", "quad", "const"]

    def test_select_by_index_name_and_list(self):
        m = self._three_output_model()
        m.set_output(1)
        assert m.evaluate(np.array([1.0, 2.0])).shape == (1, 5)
        m.set_output("const")
        assert m.active_output_names == ["const"]
        m.set_output([0, 2])
        assert m.evaluate(np.array([1.0, 2.0])).shape == (2, 5)
        assert m.active_output_names == ["lin", "const"]

    def test_selection_is_reversible(self):
        m = self._three_output_model()
        full = m.evaluate(np.array([1.0, 2.0]))
        m.set_output("quad")
        m.set_output(None)
        np.testing.assert_array_equal(m.evaluate(np.array([1.0, 2.0])), full)

    def test_constructor_argument(self):
        m = self._three_output_model(output="quad")
        assert m.active_output_names == ["quad"]
        assert m.evaluate(np.array([1.0, 2.0])).shape == (1, 5)

    def test_selection_matches_manual_indexing(self):
        m = self._three_output_model()
        full = m.evaluate(np.array([1.0, 2.0]))
        m.set_output(2)
        np.testing.assert_array_equal(m.evaluate(np.array([1.0, 2.0])), full[[2]])

    def test_batch_respects_selection(self):
        params = np.array([[1.0, 2.0], [1.5, 1.5]])
        m = self._three_output_model()
        assert m.evaluate_batch(params).shape == (2, 3, 5)
        m.set_output(1)
        assert m.evaluate_batch(params).shape == (2, 1, 5)

    def test_vectorised_batch_respects_selection(self):
        # The vectorised path bypasses the per-row loop, so it applies the selection itself.
        m = UserFunctionModel(
            func=lambda P, t: np.stack([np.vstack([p[0] * t, p[1] * t**2, p[0] + p[1] + 0 * t])
                                        for p in np.atleast_2d(P)], axis=0),
            names=["a", "b"], range=np.array([[1.0, 2.0], [1.0, 2.0]]),
            domain=np.linspace(0, 1, 5), output_names=["lin", "quad", "const"],
            vectorized=True,
        )
        params = np.array([[1.0, 2.0], [1.5, 1.5]])
        assert m.evaluate_batch(params).shape == (2, 3, 5)
        m.set_output([0, 1])
        assert m.evaluate_batch(params).shape == (2, 2, 5)

    def test_single_output_model_is_unaffected(self):
        m = UserFunctionModel(func=lambda p, t: p[0] * t, names=["a"],
                              range=np.array([[1.0, 2.0]]), domain=np.linspace(0, 1, 5))
        m.set_output(0)
        assert m.evaluate(np.array([1.5])).shape == (5,)

    def test_bad_selection_rejected(self):
        m = self._three_output_model()
        with pytest.raises(ValueError, match="unknown output"):
            m.set_output("nope")
        with pytest.raises(ValueError, match="out of range"):
            m.set_output(7)

    def test_legacy_subclass_overriding_evaluate_still_works(self):
        # A subclass written before `output` existed overrides evaluate() directly. It must
        # keep working -- bypassing selection, exactly as it did before.
        class Legacy(Model):
            names = ["a"]
            range = np.array([[1.0, 2.0]])
            nominal = np.array([1.5])
            output_names = ["out"]
            domain = None
            log_scale = np.array([False])

            def evaluate(self, params, xdata=None):
                return np.array([params[0] * 2])

        assert Legacy().evaluate(np.array([2.0])) == np.array([4.0])
