import numpy as np
import pytest

import admesh3d as am
from admesh3d import Box, FunctionDomain, Sphere, tetrahedralize
from admesh3d.api import estimate_nodes, _sample_grid

U = Box((0, 0, 0), (1, 1, 1))


@pytest.mark.parametrize("kw,msg", [
    (dict(h_max=0), "h_max"), (dict(h_max=-1), "h_max"),
    (dict(h_max=float("nan")), "h_max"), (dict(h_max=float("inf")), "h_max"),
    (dict(h_max=0.1, h_min=0), "h_min"), (dict(h_max=0.1, h_min=0.2), "h_min"),
    (dict(h_max=0.1, max_iter=0), "max_iter"),
    (dict(h_max=0.1, ttol=0), "ttol"), (dict(h_max=0.1, dptol=-1), "dptol"),
])
def test_bad_scalars(kw, msg):
    with pytest.raises(ValueError, match=msg):
        tetrahedralize(U, **kw)


@pytest.mark.parametrize("fn", [
    lambda p: np.full(len(p), np.nan), lambda p: np.full(len(p), -1.0),
    lambda p: np.ones(len(p) + 1), lambda p: np.ones((len(p), 2)),
])
def test_bad_size_fn(fn):
    with pytest.raises(ValueError, match="size_fn"):
        tetrahedralize(U, h_min=0.1, h_max=0.2, size_fn=fn)


def test_bad_domains():
    box = ((-1, -1, -1), (1, 1, 1))
    with pytest.raises(ValueError, match="no interior"):
        tetrahedralize(FunctionDomain(lambda p: np.ones(len(p)), box), h_max=0.3)
    with pytest.raises(ValueError, match="non-finite"):
        tetrahedralize(FunctionDomain(lambda p: np.full(len(p), np.nan), box), h_max=0.3)
    with pytest.raises(ValueError, match="non-finite"):
        tetrahedralize(FunctionDomain(lambda p: np.full(len(p), -np.inf), box), h_max=0.3)
    with pytest.raises(ValueError, match="degenerate"):
        tetrahedralize(FunctionDomain(lambda p: -np.ones(len(p)),
                                      ((0, 0, 0), (1, 1, 0))), h_max=0.3)
    with pytest.raises(ValueError, match="bbox"):
        FunctionDomain(lambda p: -np.ones(len(p)), ((0, 0, 0), (1, 1, np.nan)))


def test_node_cap_message_has_estimate():
    with pytest.raises(ValueError, match=r"estimated \d+ nodes"):
        tetrahedralize(U, h_max=0.1, max_nodes=100)


def test_node_estimate_close():
    _, _, pts, d = _sample_grid(U)
    est = estimate_nodes(U, pts, d, 0.1, 0.1, None)
    actual = tetrahedralize(U, h_max=0.1).n_nodes
    assert 0.5 * actual <= est <= 2 * actual


def test_huge_h_max_gives_coarse_mesh():
    with pytest.warns(UserWarning):
        m = tetrahedralize(U, h_max=5.0)
    assert m.n_nodes >= 8 and np.all(m.volumes > 0)
    assert abs(m.volumes.sum() - 1.0) < 0.1


def test_keyword_only():
    with pytest.raises(TypeError):
        tetrahedralize(U, 0.1)
    with pytest.raises(TypeError):
        am.extrude(np.zeros((3, 2)), np.array([[0, 1, 2]]), np.ones(3), 5)


def test_non_domain_rejected():
    with pytest.raises(TypeError):
        tetrahedralize("box", h_max=0.1)
