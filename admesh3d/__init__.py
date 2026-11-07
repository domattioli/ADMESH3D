"""admesh3D: experimental prototype 3D mesher in the ADMESH family."""

__version__ = "0.1.1"

from .errors import ConvergenceWarning  # noqa: E402,F401
from .domain import (Box, Cylinder, Domain, FunctionDomain, Sphere)  # noqa: E402,F401
from .quality import QualityReport  # noqa: E402,F401
from .api import Mesh, extrude, read_vtu, tetrahedralize  # noqa: E402,F401
