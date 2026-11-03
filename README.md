# admesh3D

> **EXPERIMENTAL PROTOTYPE (version 0.1.0, below 1.0).** The API may change
> without notice. Not for production use.

Tetrahedral meshing of signed-distance-function domains (3D DistMesh), sigma-layer
prism extrusion of 2D triangle meshes, a quality report, and ASCII `.vtu` output.
It is a sibling of the 2D mesh generator
[ADMESH](https://github.com/domattioli/ADMESH) (PyPI: `admesh2D`).

## Quickstart

```bash
pip install -e ".[test]"
pytest -q
```

```python
import admesh3d as am

m = am.tetrahedralize(am.Box((0,0,0),(1,1,1)), h_max=0.1)
print(m.quality_report())
m.to_vtu("cube.vtu")               # open in ParaView

ball = am.Sphere((0,0,0), 1.0) - am.Cylinder((0,0,0), 0.3, 3.0)
am.tetrahedralize(ball, h_min=0.1, h_max=0.15).to_vtu("ball.vtu")

import numpy as np
nodes = np.array([[0,0],[1,0],[1,1],[0,1]], float)
tris  = np.array([[0,1,2],[0,2,3]])
p = am.extrude(nodes=nodes, triangles=tris, depth=np.full(4,10.0), layers=5)
assert p.n_elements == 10 and p.n_nodes == 24
assert abs(p.volumes.sum() - 10.0) < 1e-9
```

## License

Apache-2.0. See [LICENSE](LICENSE).
