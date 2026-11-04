"""Smoke tests for the per-cluster / silhouette / truth-overlay plots (#6, #7, #10, Python side)."""
import numpy as np
import pytest

pytest.importorskip("matplotlib")
import matplotlib

matplotlib.use("Agg")

import gsua_csb as g


def _model():
    return g.UserFunctionModel(func=lambda p: p, names=["a", "b"],
                               range=np.array([[0.0, 10.0], [0.0, 10.0]]))


def _clustered():
    rng = np.random.default_rng(5)
    e = np.vstack([rng.normal([2.0, 2.0], 0.1, size=(60, 2)),
                   rng.normal([8.0, 8.0], 0.1, size=(60, 2))])
    return g.identifiability_analysis(_model(), e, cluster=True, seed=0)


def test_plot_silhouette_by_k():
    ax = g.plot_silhouette_by_k(_clustered(), sil_threshold=0.5)
    assert ax is not None


def test_plot_cluster_pairs_with_truth():
    axes = g.plot_cluster_pairs(_clustered(), truth=[2.0, 2.0])
    assert np.asarray(axes).shape == (2, 2)


def test_plot_parameter_cdfs_with_truth():
    axes = g.plot_parameter_cdfs(_clustered(), truth=[2.0, 2.0])
    assert len(np.atleast_1d(axes)) >= 2


def test_plot_cluster_curves():
    def decay(p, x):
        return (p[0] * np.exp(-p[1] * np.asarray(x, dtype=float)))[None, :]

    model = g.UserFunctionModel(decay, names=["a", "b"],
                                range=np.array([[0.1, 6.0], [0.05, 3.0]]),
                                domain=np.linspace(0, 2, 10))
    rng = np.random.default_rng(0)
    e = np.vstack([rng.normal([3.0, 1.0], 0.05, size=(40, 2)),
                   rng.normal([5.0, 2.0], 0.05, size=(40, 2))])
    ia = g.identifiability_analysis(model, e, cluster=True, seed=0)
    x = np.linspace(0, 2, 10)
    y = 3.0 * np.exp(-1.0 * x)
    ax = g.plot_cluster_curves(ia, model, x, y)
    assert ax is not None


def test_silhouette_plot_requires_clustering():
    rng = np.random.default_rng(1)
    ia = g.identifiability_analysis(_model(), rng.uniform(0, 10, size=(30, 2)))
    with pytest.raises(ValueError, match="silhouette-by-K"):
        g.plot_silhouette_by_k(ia)
