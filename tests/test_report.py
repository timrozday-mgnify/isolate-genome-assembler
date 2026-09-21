"""Render the report template against the fixture run_summary/, as QUARTO_REPORT does.

Needs Quarto plus the report image's Python packages, so it skips on a bare checkout. CI
runs it inside the report image.
"""

import importlib.util
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "data" / "run_summary"
SAMPLES = ("iso_pass", "iso_warn", "iso_fail")

pytestmark = pytest.mark.skipif(
    shutil.which("quarto") is None
    or any(
        importlib.util.find_spec(m) is None for m in ("papermill", "plotly", "pandas")
    ),
    reason="needs quarto, papermill, plotly and pandas (the report image)",
)


def render(tmp_path: Path, *params: str) -> str:
    template = tmp_path / "render"
    shutil.copytree(ROOT / "assets" / "report", template, dirs_exist_ok=True)
    args = ["-P", f"summary_dir:{FIXTURE}"]
    for param in params:
        args += ["-P", param]
    subprocess.run(
        ["quarto", "render", str(template / "isolate_assembly_report.qmd"), *args],
        check=True,
        capture_output=True,
    )
    return (template / "isolate_assembly_report.html").read_text()


def test_run_report_names_every_sample(tmp_path: Path) -> None:
    report = render(tmp_path)
    for sample in SAMPLES:
        assert f'id="sample-{sample}"' in report
    assert "status-fail" in report and "status-warn" in report


def test_per_sample_report_holds_only_that_sample(tmp_path: Path) -> None:
    report = render(tmp_path, "sample:iso_warn")
    assert 'id="sample-iso_warn"' in report
    assert "iso_pass" not in report and "iso_fail" not in report
