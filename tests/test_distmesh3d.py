import warnings

import numpy as np
import pytest

import admesh3d as am
from admesh3d import Box, ConvergenceWarning, Sphere, tetrahedralize
from admesh3d.quality import tet_quality


def _edges(m):
    e = m.elements[:, [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]].reshape(-1, 2)
    e = np.unique(np.sort(e, axis=1), axis=0)
    return np.linalg.norm(m.nodes[e[:, 0]] - m.nodes[e[:, 1]], axis=1)


@pytest.fixture(scope="module")
def cube():
    with warnings.catch_warnings():
        warnings.simplefilter("error", ConvergenceWarning)
        return tetrahedralize(Box((0, 0, 0), (1, 1, 1)), h_max=0.1)


@pytest.fixture(scope="module")
def sphere():
    with warnings.catch_warnings():
        warnings.simplefilter("error", ConvergenceWarning)
        return tetrahedralize(Sphere((0, 0, 0), 1.0), h_max=0.15)


def test_cube_volume_and_positive(cube):
    assert np.all(cube.volumes > 0)
    assert abs(cube.volumes.sum() - 1.0) < 0.01


def test_cube_quality(cube):
    q = cube.quality
    assert q.min() > 0.05 and q.mean() >= 0.6


def test_sphere_quality_volume_boundary(sphere):
    q = sphere.quality
    assert q.min() > 0.05 and q.mean() >= 0.6
    exact = 4 / 3 * np.pi
    assert abs(sphere.volumes.sum() - exact) / exact < 0.03
    # boundary nodes: nodes on boundary faces lie within 1% of h of the surface
    f = np.sort(sphere.elements[:, [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]]].reshape(-1, 3), axis=1)
    u, cnt = np.unique(f, axis=0, return_counts=True)
    b = np.unique(u[cnt == 1])
    r = np.linalg.norm(sphere.nodes[b], axis=1)
    assert np.all(np.abs(r - 1.0) < 0.01 * 0.15)


def test_hmin_hmax_median_edge():
    m = tetrahedralize(Box((0, 0, 0), (1, 1, 1)), h_min=0.05, h_max=0.1)
    assert abs(np.median(_edges(m)) - 0.1) / 0.1 < 0.25


def test_size_fn_graded():
    fn = lambda p: np.where(p[:, 0] < 0.5, 0.1, 0.25)  # noqa: E731
    m = tetrahedralize(Box((0, 0, 0), (1, 1, 1)), h_min=0.1, h_max=0.25, size_fn=fn)
    fine = np.sum(m.nodes[:, 0] < 0.5)
    coarse = np.sum(m.nodes[:, 0] >= 0.5)
    assert fine > coarse
    assert np.all(m.volumes > 0)


def test_size_fn_is_clipped():
    # size_fn asks for a size far below h_min; result must equal the clipped one
    big = tetrahedralize(Box((0, 0, 0), (1, 1, 1)), h_min=0.2, h_max=0.2,
                         size_fn=lambda p: np.full(len(p), 0.001))
    assert big.n_nodes < 500


def test_same_seed_identical():
    d = Box((0, 0, 0), (1, 1, 1))
    a = tetrahedralize(d, h_max=0.25, seed=3)
    b = tetrahedralize(d, h_max=0.25, seed=3)
    assert np.array_equal(a.nodes, b.nodes) and np.array_equal(a.elements, b.elements)


def test_max_iter_one_warns_and_stays_valid():
    with pytest.warns(ConvergenceWarning):
        m = tetrahedralize(Sphere((0, 0, 0), 1.0), h_max=0.3, max_iter=1)
    assert np.all(m.volumes > 0)
    assert len(np.unique(m.elements)) == m.n_nodes


@pytest.mark.slow
def test_sc001_cube_under_60s():
    import time
    t = time.time()
    tetrahedralize(Box((0, 0, 0), (1, 1, 1)), h_max=0.1)
    assert time.time() - t < 60


def _banded_edge_medians(m, fn, h_min, h_max, centers, bands):
    """Per band: (lo, hi, median edge length, median clipped target), all unique edges."""
    e = m.elements[:, [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]].reshape(-1, 2)
    e = np.unique(np.sort(e, axis=1), axis=0)
    a, b = m.nodes[e[:, 0]], m.nodes[e[:, 1]]
    length = np.linalg.norm(a - b, axis=1)
    mid = 0.5 * (a + b)
    d = np.min(np.linalg.norm(mid[:, None, :] - centers[None], axis=2), axis=1)
    target = np.clip(fn(mid), h_min, h_max)
    out = []
    for lo, hi in bands:
        s = (d >= lo) & (d < hi)
        assert s.sum() > 20, (lo, hi, int(s.sum()))
        out.append((lo, hi, np.median(length[s]), np.median(target[s])))
    return out


def test_size_fn_absolute_edge_lengths_cube_corners():
    corners = np.array([[x, y, z] for x in (0, 1) for y in (0, 1) for z in (0, 1)], float)

    def fn(p):
        d = np.min(np.linalg.norm(p[:, None, :] - corners[None], axis=2), axis=1)
        return 0.035 + 0.35 * d

    m = tetrahedralize(Box((0, 0, 0), (1, 1, 1)), h_min=0.035, h_max=0.15, size_fn=fn)
    bands = [(0, 0.15), (0.15, 0.3), (0.3, 0.5), (0.5, 0.9)]
    for lo, hi, med, tgt in _banded_edge_medians(m, fn, 0.035, 0.15, corners, bands):
        assert abs(med - tgt) / tgt < 0.25, (lo, hi, med, tgt)


def test_size_fn_absolute_edge_lengths_sphere_surface():
    def fn(p):
        return 0.06 + 0.3 * (1.0 - np.linalg.norm(p, axis=1))

    m = tetrahedralize(Sphere((0, 0, 0), 1.0), h_min=0.06, h_max=0.25, size_fn=fn)
    # Bands by edge-midpoint distance from the centre.
    centers = np.zeros((1, 3))
    bands = [(0, 0.3), (0.3, 0.6), (0.6, 0.85)]
    for lo, hi, med, tgt in _banded_edge_medians(m, fn, 0.06, 0.25, centers, bands):
        assert abs(med - tgt) / tgt < 0.25, (lo, hi, med, tgt)


def test_callback_snapshots():
    import admesh3d
    snaps = []
    admesh3d.tetrahedralize(admesh3d.Box((0, 0, 0), (1, 1, 1)), h_max=0.3, max_iter=40,
                            callback=lambda it, p: snaps.append((it, p.shape)),
                            callback_every=5)
    assert snaps[0][0] == 0
    assert all(it % 5 == 0 for it, _ in snaps)
    assert all(s[1] == 3 for _, s in snaps)
