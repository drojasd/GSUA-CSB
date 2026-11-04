"""Per-cluster summary with identifiability verdict (cluster_summary / gsua_clustersummary, #8)."""
import numpy as np

import gsua_csb as g


def _model():
    return g.UserFunctionModel(func=lambda p: p, names=["a", "b"],
                               range=np.array([[0.0, 10.0], [0.0, 10.0]]))


def test_two_cluster_summary_rows_and_verdicts():
    rng = np.random.default_rng(5)
    e = np.vstack([rng.normal([2.0, 2.0], 0.05, size=(60, 2)),
                   rng.normal([8.0, 8.0], 0.05, size=(60, 2))])
    ia = g.identifiability_analysis(_model(), e, cluster=True, seed=0)
    rows = g.cluster_summary(ia)

    assert len(rows) == 2 * 2                       # 2 clusters x 2 params
    assert all(r.verdict == "identifiable" for r in rows)   # tight basins
    r = rows[0]
    assert r.n == 60
    assert r.median_ci_low <= r.median <= r.median_ci_high
    assert r.p2_5 <= r.median <= r.p97_5


def test_flat_parameter_is_not_identifiable():
    rng = np.random.default_rng(1)
    e = np.column_stack([5.0 + 0.01 * rng.standard_normal(200), rng.uniform(0, 10, 200)])
    ia = g.identifiability_analysis(_model(), e)
    by = {r.parameter: r for r in g.cluster_summary(ia)}

    assert by["a"].verdict == "identifiable"
    assert by["b"].verdict == "not identifiable"
    assert by["b"].pct_width_rel > 0.5
