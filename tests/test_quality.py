import numpy as np
import pytest

from admesh3d.quality import (tet_quality, tet_signed_volumes, wedge_quality,
                              wedge_signed_volumes)

UNIT_TET = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
UNIT_WEDGE = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0],
                       [0, 0, 1], [1, 0, 1], [0, 1, 1]], float)


def test_tet_signed_volume_sign():
    e = np.array([[0, 1, 2, 3]])
    assert tet_signed_volumes(UNIT_TET, e)[0] == pytest.approx(1 / 6)
    mirror = UNIT_TET * np.array([1, 1, -1])
    assert tet_signed_volumes(mirror, e)[0] == pytest.approx(-1 / 6)


def test_wedge_signed_volume_sign():
    e = np.array([[0, 1, 2, 3, 4, 5]])
    assert wedge_signed_volumes(UNIT_WEDGE, e)[0] == pytest.approx(0.5)
    mirror = UNIT_WEDGE * np.array([1, 1, -1])
    assert wedge_signed_volumes(mirror, e)[0] == pytest.approx(-0.5)


def regular_tet():
    return np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], float)


def test_regular_tet_quality_is_one():
    p = regular_tet()
    if tet_signed_volumes(p, np.array([[0, 1, 2, 3]]))[0] < 0:
        p = p[[0, 2, 1, 3]]
    assert tet_quality(p, np.array([[0, 1, 2, 3]]))[0] == pytest.approx(1.0)


def test_sliver_tet_quality_near_zero():
    p = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0.5, 0.5, 1e-4]], float)
    assert tet_quality(p, np.array([[0, 1, 2, 3]]))[0] < 0.01


def test_regular_prism_quality_is_one():
    s = 1.0
    tri = np.array([[0, 0], [s, 0], [s / 2, s * np.sqrt(3) / 2]])
    p = np.vstack([np.c_[tri, np.zeros(3)], np.c_[tri, np.full(3, s)]])
    q = wedge_quality(p, np.array([[0, 1, 2, 3, 4, 5]]))
    assert q[0] == pytest.approx(1.0)


def test_inverted_wedge_quality_zero():
    mirror = UNIT_WEDGE * np.array([1, 1, -1])
    assert wedge_quality(mirror, np.array([[0, 1, 2, 3, 4, 5]]))[0] == 0.0


def _cleaned(domain_h, domain):
    import warnings
    from admesh3d import ConvergenceWarning, tetrahedralize
    with warnings.catch_warnings():
        warnings.simplefilter("error", ConvergenceWarning)
        return tetrahedralize(domain, h_max=domain_h)


def test_report_fields_and_str():
    import admesh3d as am
    m = am.tetrahedralize(am.Box((0, 0, 0), (1, 1, 1)), h_max=0.5)
    r = m.quality_report()
    q = np.sort(m.quality)
    assert r.q_min == pytest.approx(q.min()) and r.q_mean == pytest.approx(q.mean())
    k = max(1, int(np.ceil(0.05 * len(q))))
    assert r.q_worst5 == pytest.approx(q[:k].mean())
    assert r.total_volume == pytest.approx(1.0)
    assert 0 <= r.q_mean <= 1 and r.element_type == "tetra"
    assert "q_mean" in str(r)
    with pytest.raises(Exception):
        r.q_min = 0.0


def test_sc003_cube_and_sphere():
    import admesh3d as am
    for m in (_cleaned(0.1, am.Box((0, 0, 0), (1, 1, 1))),
              _cleaned(0.15, am.Sphere((0, 0, 0), 1.0))):
        r = m.quality_report()
        assert r.q_min > 0.05 and r.q_mean >= 0.6


def test_wedge_report_sc007():
    import admesh3d as am
    m = am.extrude(nodes=np.array([[0, 0], [1, 0], [1, 1], [0, 1]], float),
                   triangles=np.array([[0, 1, 2], [0, 2, 3]]),
                   depth=np.full(4, 10.0), layers=5)
    r = m.quality_report()
    assert r.element_type == "wedge" and r.q_min >= 0.1
