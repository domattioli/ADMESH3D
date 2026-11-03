import numpy as np
import pytest

import admesh3d as am

admesh = pytest.importorskip("admesh")


def test_admesh2d_mesh_extrudes():
    # Duck typing only: any object with .nodes/.elements/.depth works.  Use the
    # admesh2D package only to prove it is importable alongside admesh3d.
    from types import SimpleNamespace
    nodes = np.array([[0, 0], [1, 0], [1, 1], [0, 1]], float)
    m2 = SimpleNamespace(nodes=nodes, elements=np.array([[0, 1, 2], [0, 2, 3]]),
                         bathymetry=np.full(4, -3.0))
    m = am.extrude(m2, layers=3)
    assert m.n_nodes == 16 and m.n_elements == 6 and np.all(m.volumes > 0)
