import numpy as np
import pytest

from admesh3d.domain import (Box, Cylinder, FunctionDomain, Sphere)


def test_box_sdf_sign_and_bbox():
    b = Box((0, 0, 0), (1, 2, 3))
    d = b.sdf(np.array([[0.5, 1, 1.5], [2, 1, 1], [1, 2, 3]]))
    assert d[0] < 0 and d[1] > 0 and abs(d[2]) < 1e-12
    assert b.bbox == ((0, 0, 0), (1, 2, 3))
    assert b.volume == pytest.approx(6.0)


def test_box_fixed_points_are_corners():
    fp = Box((0, 0, 0), (1, 1, 1)).fixed_points()
    assert fp.shape == (8, 3)
    assert len({tuple(r) for r in fp}) == 8
    assert Sphere((0, 0, 0), 1).fixed_points().shape == (0, 3)
    assert Cylinder((0, 0, 0), 1, 2).fixed_points().shape == (0, 3)


def test_sphere_and_cylinder():
    s = Sphere((1, 1, 1), 2.0)
    assert s.sdf(np.array([[1, 1, 1]]))[0] == pytest.approx(-2.0)
    assert s.sdf(np.array([[3, 1, 1]]))[0] == pytest.approx(0.0)
    assert s.volume == pytest.approx(4 / 3 * np.pi * 8)
    c = Cylinder((0, 0, 0), 1.0, 2.0)
    assert c.sdf(np.array([[0, 0, 0]]))[0] < 0
    assert c.sdf(np.array([[0, 0, 1.5]]))[0] > 0
    assert c.sdf(np.array([[2, 0, 0]]))[0] > 0
    assert c.volume == pytest.approx(2 * np.pi)
    assert c.bbox == ((-1, -1, -1), (1, 1, 1))
    cx = Cylinder((0, 0, 0), 1.0, 4.0, axis="x")
    assert cx.bbox == ((-2, -1, -1), (2, 1, 1))


def test_set_operations():
    a = Sphere((0, 0, 0), 1.0)
    b = Sphere((1, 0, 0), 1.0)
    p = np.array([[-0.9, 0, 0], [0.5, 0, 0], [1.9, 0, 0], [5, 0, 0]])
    assert ((a | b).sdf(p) < 0).tolist() == [True, True, True, False]
    assert ((a & b).sdf(p) < 0).tolist() == [False, True, False, False]
    assert ((a - b).sdf(p) < 0).tolist() == [True, False, False, False]
    assert (a | b).bbox == ((-1, -1, -1), (2, 1, 1))
    assert (a & b).bbox == ((0, -1, -1), (1, 1, 1))
    assert (a - b).bbox == a.bbox
    assert (a - b).fixed_points().shape == (0, 3)


def test_function_domain():
    d = FunctionDomain(lambda p: np.linalg.norm(p, axis=1) - 1,
                       ((-1, -1, -1), (1, 1, 1)))
    assert d.sdf(np.zeros((1, 3)))[0] == pytest.approx(-1)
    assert d.volume is None


def test_validation():
    with pytest.raises(ValueError, match="radius"):
        Sphere((0, 0, 0), 0)
    with pytest.raises(ValueError, match="radius"):
        Cylinder((0, 0, 0), -1, 1)
    with pytest.raises(ValueError, match="height"):
        Cylinder((0, 0, 0), 1, 0)
    with pytest.raises(ValueError, match="hi"):
        Box((0, 0, 0), (1, 0, 1))
    with pytest.raises(ValueError, match="finite"):
        Box((0, 0, 0), (1, np.inf, 1))
    with pytest.raises(ValueError, match="bbox"):
        FunctionDomain(lambda p: p[:, 0], ((0, 0, 0), (np.inf, 1, 1)))
    with pytest.raises(ValueError, match="axis"):
        Cylinder((0, 0, 0), 1, 1, axis="w")
