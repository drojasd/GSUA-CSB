"""Tests for ``design_matrix(method="joint")`` -- correlation-preserving sampling.

Mirrors the MATLAB suite (``tests/gsuaDmatrixJoint*Test.m``). The fixture pool lies exactly on
the hyperbola ``a*b = 2``, the manifold that ``y = a*b*exp(-k*t)`` actually identifies, so
"did the sampler leave the manifold" is a measurement rather than a judgement call.
"""

import warnings

import numpy as np
import pytest

from gsua_csb import UserFunctionModel, design_matrix

MANIFOLD = 2.0
XDATA = np.linspace(0.0, 5.0, 20)


def _model(fix_k=False):
    rng = np.array([[0.3, 3.0], [0.3, 3.0], [0.1, 1.2]])
    if fix_k:
        rng[2, :] = 0.5
    return UserFunctionModel(
        func=lambda p, xd: (p[0] * p[1] * np.exp(-p[2] * xd))[None, :],
        names=["a", "b", "k"],
        range=rng,
        nominal=np.array([1.0, 2.0, 0.5]),
        domain=XDATA,
        output_names=["y"],
    )


def _pool(n_pool=40):
    """Accepted-estimate ensemble lying exactly on a*b = 2.

    Deterministic (sin, not a Generator) so the fixture cannot perturb seeded draws.
    """
    a = np.linspace(0.5, 2.5, n_pool)
    return np.column_stack([a, MANIFOLD / a, 0.5 + 0.02 * np.sin(np.arange(1, n_pool + 1))])


def _leak(M):
    """Distance of each draw from the identified manifold."""
    return np.abs(M[:, 0] * M[:, 1] - MANIFOLD)


class TestContract:
    """Joint sampling must honor design_matrix's existing output contract."""

    @pytest.mark.parametrize("jt", ["bootstrap", "smooth_bootstrap", "gaussian"])
    def test_shape_and_finiteness(self, jt):
        M = design_matrix(_model(), 250, "joint", 7, joint_type=jt, pool=_pool())
        assert M.shape == (250, 3)
        assert np.all(np.isfinite(M))

    def test_every_bootstrap_draw_is_a_pool_row(self):
        # Drawing vectors intact is what preserves EVERY dependence in the pool,
        # not merely its pairwise correlations.
        pool = _pool()
        M = design_matrix(_model(), 200, "joint", 3, pool=pool)
        worst = max(np.min(np.linalg.norm(pool - row, axis=1)) for row in M)
        assert worst < 1e-12

    def test_fixed_parameters_stay_fixed(self):
        M = design_matrix(_model(fix_k=True), 100, "joint", 11, pool=_pool())
        np.testing.assert_allclose(M[:, 2], 0.5)
        assert np.ptp(M[:, 0]) > 0

    def test_missing_pool_raises(self):
        with pytest.raises(ValueError, match="needs an ensemble"):
            design_matrix(_model(), 10, "joint", 0)

    def test_wrong_pool_shape_raises(self):
        with pytest.raises(ValueError, match=r"pool must be \(n_pool, 3\)"):
            design_matrix(_model(), 10, "joint", 0, pool=np.ones((20, 2)))

    def test_unknown_joint_type_raises(self):
        with pytest.raises(ValueError, match="Unknown joint_type"):
            design_matrix(_model(), 10, "joint", 0, joint_type="mvn", pool=_pool())

    def test_joint_options_rejected_without_joint_method(self):
        with pytest.raises(ValueError, match="only appl"):
            design_matrix(_model(), 10, "latin_hypercube", 0, pool=_pool())

    @pytest.mark.parametrize("method", ["latin_hypercube", "uniform", "sobol"])
    def test_marginal_methods_unaffected(self, method):
        model = _model()
        M = design_matrix(model, 32, method=method, seed=0)
        assert M.shape == (32, 3)
        assert np.all(M >= model.range[:, 0]) and np.all(M <= model.range[:, 1])


class TestCorrelation:
    """The structure marginal sampling destroys must survive joint sampling."""

    def test_bootstrap_preserves_pool_correlation(self):
        pool = _pool()
        M = design_matrix(_model(), 4000, "joint", 21, pool=pool)
        assert np.corrcoef(M[:, 0], M[:, 1])[0, 1] == pytest.approx(
            np.corrcoef(pool[:, 0], pool[:, 1])[0, 1], abs=0.05
        )

    def test_marginal_sampling_destroys_correlation(self):
        M = design_matrix(_model(), 4000, "latin_hypercube", 21)
        assert abs(np.corrcoef(M[:, 0], M[:, 1])[0, 1]) < 0.15

    def test_gaussian_also_preserves_linear_correlation(self):
        # The boundary of correlation as a criterion: the Gaussian mode matches the
        # coefficient closely and is still off-manifold (see TestManifold).
        pool = _pool()
        M = design_matrix(_model(), 4000, "joint", 21, joint_type="gaussian", pool=pool)
        assert np.corrcoef(M[:, 0], M[:, 1])[0, 1] == pytest.approx(
            np.corrcoef(pool[:, 0], pool[:, 1])[0, 1], abs=0.05
        )

    def test_draws_inherit_pool_spread_not_model_range(self):
        # identifiability_analysis leaves `range` as a CI of the MEDIAN, narrower than the
        # pool and narrowing further as the pool grows. Joint draws must reflect the pool's
        # own spread -- which only holds because nothing clips them back by default.
        pool = _pool()
        model = _model()
        model.range[0, :] = [1.4, 1.6]  # a deliberately too-narrow "CI of the median"
        M = design_matrix(model, 2000, "joint", 33, pool=pool)
        assert M[:, 0].std() == pytest.approx(pool[:, 0].std(), rel=0.15)
        assert M[:, 0].max() > 1.6
        assert M[:, 0].min() < 1.4

    def test_clip_applied_when_requested(self):
        bounds = np.array([[1.0, 2.0], [0.3, 3.0], [0.1, 1.2]])
        M = design_matrix(_model(), 500, "joint", 33, pool=_pool(), clip=bounds)
        assert M[:, 0].min() >= 1.0
        assert M[:, 0].max() <= 2.0


class TestManifold:
    """The analytic control: the manifold is known in closed form, so leakage is measurable.

    This is what justifies bootstrap being the default. The manifold is CURVED, so a Gaussian
    fitted to it must place mass beside it -- preserving pairwise linear correlation is
    demonstrably not sufficient to stay on it.
    """

    def test_bootstrap_stays_on_the_manifold_exactly(self):
        M = design_matrix(_model(), 1000, "joint", 13, pool=_pool())
        assert _leak(M).max() < 1e-10

    def test_gaussian_leaves_the_manifold(self):
        Mg = design_matrix(_model(), 1000, "joint", 13, joint_type="gaussian", pool=_pool())
        Mb = design_matrix(_model(), 1000, "joint", 13, pool=_pool())
        assert _leak(Mg).mean() > 1e-3
        assert _leak(Mg).mean() > 1e3 * max(_leak(Mb).mean(), np.finfo(float).eps)

    def test_marginal_leaks_most_of_all(self):
        Ml = design_matrix(_model(), 1000, "latin_hypercube", 13)
        Mg = design_matrix(_model(), 1000, "joint", 13, joint_type="gaussian", pool=_pool())
        assert _leak(Ml).mean() > _leak(Mg).mean()

    def test_smooth_bootstrap_stays_closer_than_gaussian(self):
        Ms = design_matrix(
            _model(), 1000, "joint", 13, joint_type="smooth_bootstrap", pool=_pool()
        )
        Mg = design_matrix(_model(), 1000, "joint", 13, joint_type="gaussian", pool=_pool())
        assert _leak(Ms).mean() < _leak(Mg).mean()

    def test_smooth_bootstrap_reproduces_pool_spread(self):
        # The variance correction exists so smoothing reproduces the ensemble covariance
        # instead of inflating it by (1+h^2).
        pool = _pool()
        Ms = design_matrix(
            _model(), 20000, "joint", 99, joint_type="smooth_bootstrap", pool=pool
        )
        assert Ms[:, 0].std() / pool[:, 0].std() == pytest.approx(1.0, abs=0.06)


class TestGuard:
    """A pool too small to trust must be refused, not used.

    Two points always correlate at exactly +-1, so a joint band from such a pool is
    confidently wrong rather than approximately right -- a worse failure than the one joint
    sampling fixes, because it looks like success.
    """

    def test_tiny_pool_warns_and_falls_back(self):
        with pytest.warns(UserWarning, match="fewer than min_pool_n"):
            M = design_matrix(_model(), 400, "joint", 4, pool=_pool()[:3])
        assert M.shape == (400, 3)
        # Degrades honestly: no pretended correlation, and back inside the model range.
        assert abs(np.corrcoef(M[:, 0], M[:, 1])[0, 1]) < 0.2
        assert M[:, 0].min() >= 0.3 and M[:, 0].max() <= 3.0

    def test_adequate_pool_does_not_warn(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            design_matrix(_model(), 50, "joint", 4, pool=_pool()[:5])

    def test_min_pool_n_is_configurable(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            design_matrix(_model(), 50, "joint", 4, pool=_pool()[:3], min_pool_n=3)
        with pytest.warns(UserWarning, match="fewer than min_pool_n"):
            design_matrix(_model(), 50, "joint", 4, pool=_pool()[:8], min_pool_n=20)

    def test_non_finite_pool_rows_are_dropped(self):
        pool = _pool()
        pool[3, 1] = np.nan
        with pytest.warns(UserWarning, match="non-finite"):
            M = design_matrix(_model(), 200, "joint", 4, pool=pool)
        assert np.all(np.isfinite(M))

    def test_all_non_finite_pool_raises(self):
        pool = _pool()
        pool[:, 0] = np.nan
        with pytest.raises(ValueError, match="every pool row"):
            design_matrix(_model(), 50, "joint", 4, pool=pool)

    def test_rank_deficient_covariance_is_usable(self):
        # Fewer pool rows than free parameters is routine for an expensive model; the
        # eigendecomposition must sample the non-null subspace rather than fail.
        pool = _pool(6)
        M = design_matrix(_model(), 100, "joint", 4, joint_type="gaussian", pool=pool)
        assert np.all(np.isfinite(M))


class TestSeeding:
    def test_same_seed_reproduces_exactly(self):
        for jt in ["bootstrap", "smooth_bootstrap", "gaussian"]:
            a = design_matrix(_model(), 200, "joint", 2024, joint_type=jt, pool=_pool())
            b = design_matrix(_model(), 200, "joint", 2024, joint_type=jt, pool=_pool())
            np.testing.assert_array_equal(a, b)

    def test_different_seeds_differ(self):
        a = design_matrix(_model(), 200, "joint", 1, pool=_pool())
        b = design_matrix(_model(), 200, "joint", 2, pool=_pool())
        assert not np.array_equal(a, b)

    def test_accepts_a_generator(self):
        rng = np.random.default_rng(5)
        M = design_matrix(_model(), 50, "joint", rng, pool=_pool())
        assert M.shape == (50, 3)
