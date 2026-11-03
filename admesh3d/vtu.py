"""ASCII VTK XML unstructured grid (.vtu) writer and reader.

Only tetrahedra (VTK_TETRA = 10) and wedges (VTK_WEDGE = 13) are supported.
Wedges are written in the internal node order, which already matches the VTK
wedge convention, so no permutation is applied.
"""
from __future__ import annotations

import io
import xml.etree.ElementTree as ET

import numpy as np

VTK_TETRA = 10
VTK_WEDGE = 13
_TYPE_TO_VTK = {"tetra": VTK_TETRA, "wedge": VTK_WEDGE}
_VTK_TO_TYPE = {v: k for k, v in _TYPE_TO_VTK.items()}
_NODES_PER_CELL = {"tetra": 4, "wedge": 6}


def _text(arr: np.ndarray, fmt: str) -> str:
    buf = io.StringIO()
    np.savetxt(buf, np.atleast_2d(arr), fmt=fmt)
    return "\n" + buf.getvalue()


def _array(parent, name, dtype, values, fmt, components=None):
    attrs = {"type": dtype, "format": "ascii"}
    if name:
        attrs["Name"] = name
    if components:
        attrs["NumberOfComponents"] = str(components)
    el = ET.SubElement(parent, "DataArray", attrs)
    el.text = _text(values, fmt)
    return el


def write_vtu(path, nodes, elements, element_type, quality=None) -> None:
    """Write a mesh as an ASCII .vtu file."""
    if element_type not in _TYPE_TO_VTK:
        raise ValueError(f"unsupported element_type {element_type!r}")
    nodes = np.asarray(nodes, dtype=float)
    elements = np.asarray(elements, dtype=np.int64)
    k = _NODES_PER_CELL[element_type]
    root = ET.Element("VTKFile", {"type": "UnstructuredGrid", "version": "0.1",
                                  "byte_order": "LittleEndian"})
    grid = ET.SubElement(root, "UnstructuredGrid")
    piece = ET.SubElement(grid, "Piece", {"NumberOfPoints": str(len(nodes)),
                                          "NumberOfCells": str(len(elements))})
    points = ET.SubElement(piece, "Points")
    _array(points, None, "Float64", nodes, "%.17g", components=3)
    cells = ET.SubElement(piece, "Cells")
    _array(cells, "connectivity", "Int64", elements.reshape(-1), "%d")
    offsets = np.arange(1, len(elements) + 1, dtype=np.int64) * k
    _array(cells, "offsets", "Int64", offsets, "%d")
    types = np.full(len(elements), _TYPE_TO_VTK[element_type], dtype=np.int64)
    _array(cells, "types", "UInt8", types, "%d")
    if quality is not None:
        cell_data = ET.SubElement(piece, "CellData", {"Scalars": "quality"})
        _array(cell_data, "quality", "Float64", np.asarray(quality, dtype=float), "%.17g")
    ET.indent(root)
    ET.ElementTree(root).write(path, encoding="utf-8", xml_declaration=True)


def _read_array(parent, name, dtype):
    for el in parent.iter("DataArray"):
        if el.get("Name") == name or (name is None and el.get("Name") is None):
            text = (el.text or "").split()
            return np.array(text, dtype=dtype)
    raise ValueError(f"vtu file has no DataArray {name!r}")


def read_vtu(path):
    """Read an ASCII .vtu file; returns (nodes, elements, element_type, quality)."""
    root = ET.parse(path).getroot()
    piece = root.find("./UnstructuredGrid/Piece")
    if piece is None:
        raise ValueError("not an UnstructuredGrid .vtu file")
    n_cells = int(piece.get("NumberOfCells"))
    nodes = _read_array(piece.find("Points"), None, float).reshape(-1, 3)
    cells = piece.find("Cells")
    conn = _read_array(cells, "connectivity", np.int64)
    types = _read_array(cells, "types", np.int64)
    if n_cells == 0 or len(types) == 0:
        raise ValueError("vtu file has no cells")
    unique = set(types.tolist())
    if len(unique) != 1 or next(iter(unique)) not in _VTK_TO_TYPE:
        raise ValueError(f"unsupported or mixed cell types {sorted(unique)}")
    element_type = _VTK_TO_TYPE[next(iter(unique))]
    elements = conn.reshape(n_cells, _NODES_PER_CELL[element_type])
    quality = None
    cell_data = piece.find("CellData")
    if cell_data is not None:
        try:
            quality = _read_array(cell_data, "quality", float)
        except ValueError:
            quality = None
    return nodes, elements, element_type, quality
