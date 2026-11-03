import time
from types import SimpleNamespace

import numpy as np
import pytest

import admesh3d as am

NODES = np.array([[0, 0], [1, 0], [1, 1], [0, 1]], float)
TRIS = np.array([[0, 1, 2], [0, 2, 3]])


def test_counts_and_volume():
    m = am.extrude(nodes=NODES, triangles=TRIS, depth=np.full(4, 10.0), layers=5)
    assert m.n_nodes == 4 * 6 and m.n_elements == 2 * 5
    assert np.all(m.volumes > 0)
    assert m.volumes.sum() == pytest.approx(10.0)
    assert m.nodes[:, 2].min() == pytest.approx(-10.0) and m.nodes[:, 2].max() == 0.0
    assert m.element_type == "wedge"


def test_cw_triangles_are_reoriented():
    m = am.extrude(nodes=NODES, triangles=TRIS[:, ::-1], depth=np.full(4, 2.0), layers=2)
    assert np.all(m.volumes > 0)


def test_variable_depth_volume():
    depth = np.array([1.0, 2.0, 3.0, 2.0])
    m = am.extrude(nodes=NODES, triangles=TRIS, depth=depth, layers=3)
    assert m.volumes.sum() == pytest.approx(2.0)  # mean depth over unit square


def test_errors():
    with pytest.raises(ValueError, match="node 2 has depth"):
        am.extrude(nodes=NODES, triangles=TRIS, depth=np.array([1, 1, 0, -1.0]), layers=2)
    for bad in (0, 1.5, "3", True):
        with pytest.raises(ValueError, match="layers"):
            am.extrude(nodes=NODES, triangles=TRIS, depth=np.ones(4), layers=bad)
    with pytest.raises(ValueError, match="shape"):
        am.extrude(nodes=NODES, triangles=TRIS, depth=np.ones(3), layers=2)
    with pytest.raises(ValueError, match="depth"):
        am.extrude(nodes=NODES, triangles=TRIS, layers=2)


def test_duck_typed_mesh():
    n3 = np.c_[NODES, np.zeros(4)]
    m = am.extrude(SimpleNamespace(nodes=n3, elements=TRIS, depth=np.full(4, 4.0)), layers=2)
    assert m.volumes.sum() == pytest.approx(4.0)
    m = am.extrude(SimpleNamespace(nodes=NODES, elements=TRIS, bathymetry=np.full(4, -4.0)),
                   layers=2)
    assert m.volumes.sum() == pytest.approx(4.0)
    with pytest.raises(ValueError):
        am.extrude(SimpleNamespace(nodes=NODES, elements=TRIS, depth=None), layers=2)
    with pytest.raises(ValueError):
        am.extrude(SimpleNamespace(nodes=NODES, elements=TRIS), layers=2)


def test_sc006_large_extrusion():
    n = 75
    xs, ys = np.meshgrid(np.arange(n, dtype=float), np.arange(n, dtype=float))
    nodes = np.c_[xs.ravel(), ys.ravel()]
    idx = np.arange(n * n).reshape(n, n)
    a, b, c, d = idx[:-1, :-1], idx[:-1, 1:], idx[1:, 1:], idx[1:, :-1]
    tris = np.concatenate([np.stack([a, b, c], -1).reshape(-1, 3),
                           np.stack([a, c, d], -1).reshape(-1, 3)])
    assert len(tris) >= 10_000
    t = time.time()
    m = am.extrude(nodes=nodes, triangles=tris, depth=np.full(len(nodes), 5.0), layers=10)
    assert time.time() - t < 5
    assert m.volumes.sum() == pytest.approx((n - 1) ** 2 * 5.0, rel=1e-3)
