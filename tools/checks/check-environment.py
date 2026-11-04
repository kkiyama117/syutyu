"""Run the environment check notebook using this project's virtual environment.

Run from the repository root:

    uv run --project library --locked library/tools/checks/check-environment.py

The executed notebook is written under library/local/validation/ (Git-ignored), never
into the notebook's source directory. No credentials are read and no server is
started.
"""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys


LIBRARY = Path(__file__).resolve().parents[2]
NOTEBOOK = "tools/checks/environment-check.ipynb"
OUTPUT_DIR = "local/validation"


def main() -> None:
    argparse.ArgumentParser(description=__doc__).parse_args()
    environment = Path(os.environ.get("UV_PROJECT_ENVIRONMENT", LIBRARY / ".venv"))
    if not environment.is_absolute():
        environment = LIBRARY / environment
    if Path(sys.prefix).resolve() != environment.resolve():
        raise RuntimeError("Run from the repository root: uv run --project library --locked library/tools/checks/check-environment.py")
    output_name = (
        "environment-check.executed.ipynb"
        if environment.resolve() == LIBRARY / ".venv"
        else "environment-check.container.executed.ipynb"
    )

    subprocess.run(["uv", "lock", "--check"], cwd=LIBRARY, check=True)
    subprocess.run(
        ["uv", "pip", "check", "--python", sys.executable], cwd=LIBRARY, check=True
    )
    subprocess.run(
        [sys.executable, "-m", "jupyterlab", "--version"], cwd=LIBRARY, check=True
    )
    (LIBRARY / OUTPUT_DIR).mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            sys.executable,
            "-m",
            "nbconvert",
            "--to", "notebook",
            "--execute", NOTEBOOK,
            "--output", output_name,
            "--output-dir", OUTPUT_DIR,
            "--ExecutePreprocessor.timeout=120",
        ],
        cwd=LIBRARY,
        check=True,
    )

    output_path = LIBRARY / OUTPUT_DIR / output_name
    notebook = json.loads(output_path.read_text(encoding="utf-8"))
    cells = [cell for cell in notebook["cells"] if cell["cell_type"] == "code"]
    if not cells or any(cell["execution_count"] is None for cell in cells):
        raise RuntimeError("Some notebook cells were not executed")
    outputs = [output for cell in cells for output in cell["outputs"]]
    if any(output["output_type"] == "error" for output in outputs):
        raise RuntimeError("Notebook contains an execution error")
    if not any("Calculation OK" in "".join(output.get("text", [])) for output in outputs):
        raise RuntimeError("Numerical check did not complete")
    if not any("image/png" in output.get("data", {}) for output in outputs):
        raise RuntimeError("Notebook did not produce a PNG plot")
    print(f"PASS: {sys.executable}; calculation and plot verified")
    print(f"Executed notebook: {output_path}")


if __name__ == "__main__":
    main()
