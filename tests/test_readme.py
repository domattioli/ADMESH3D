import re
from pathlib import Path

import pytest


@pytest.mark.slow
def test_readme_first_python_block_runs(tmp_path, monkeypatch):
    text = (Path(__file__).parents[1] / "README.md").read_text()
    code = re.search(r"```python\n(.*?)```", text, re.S).group(1)
    monkeypatch.chdir(tmp_path)
    exec(compile(code, "README.md", "exec"), {})
    assert (tmp_path / "cube.vtu").exists()
