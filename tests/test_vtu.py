import numpy as np
import pytest

import admesh3d as am


def _check(m, tmp_path):
    f = tmp_path / "m.vtu"
    m.to_vtu(f)
    r = am.read_vtu(f)
    assert r.element_type == m.element_type
    assert r.n_nodes == m.n_nodes and r.n_elements == m.n_elements
    assert np.array_equal(r.elements, m.elements)
    assert np.allclose(r.nodes, m.nodes, rtol=0, atol=1e-15)
    assert 'Name="quality"' in f.read_text()
    return r


def test_tetra_roundtrip(tmp_path):
    _check(am.tetrahedralize(am.Box((0, 0, 0), (1, 1, 1)), h_max=0.3), tmp_path)


def test_wedge_roundtrip_and_vtk_orientation(tmp_path):
    nodes = np.array([[0, 0], [1, 0], [1, 1], [0, 1]], float)
    m = am.extrude(nodes=nodes, triangles=np.array([[0, 1, 2], [0, 2, 3]]),
                   depth=np.full(4, 2.0), layers=2)
    r = _check(m, tmp_path)
    p = r.nodes[r.elements]
    # VTK parametric Jacobian det[p1-p0, p2-p0, p3-p0] must be positive
    det = np.einsum("ij,ij->i", np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]),
                    p[:, 3] - p[:, 0])
    assert np.all(det > 0)
    n = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
    assert np.all(np.einsum("ij,ij->i", n, p[:, 3] - p[:, 0]) > 0)


def test_bad_file(tmp_path):
    f = tmp_path / "x.vtu"
    f.write_text("<VTKFile type='PolyData'/>")
    with pytest.raises(ValueError):
        am.read_vtu(f)
