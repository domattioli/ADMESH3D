import re
from pathlib import Path

import admesh3d


def test_version_below_one_and_matches_pyproject():
    text = (Path(__file__).parents[1] / "pyproject.toml").read_text()
    pv = re.search(r'^version\s*=\s*"([^"]+)"', text, re.M).group(1)
    assert admesh3d.__version__ == pv
    assert int(pv.split(".")[0]) < 1
