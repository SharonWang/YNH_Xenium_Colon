import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NOTEBOOKS = ROOT / "notebooks"
REGIONS = [f"Region_{i}" for i in range(1, 7)]


def source_text(notebook):
    return "\n".join("".join(cell.get("source", [])) for cell in notebook["cells"])


def test_region_notebooks():
    for index, region in enumerate(REGIONS, start=1):
        path = NOTEBOOKS / f"01_QC_Region{index}.ipynb"
        assert path.exists(), path
        notebook = json.loads(path.read_text(encoding="utf-8"))
        text = source_text(notebook)
        assert notebook["metadata"]["kernelspec"]["name"] == "ir"
        assert f'REGION_ID <- "{region}"' in text
        assert "FULL_HPC" in text and "LOCAL_SUBSET" in text
        assert "fixed_cell_qc_thresholds.tsv" in text
        assert "primary_include" in text
        for heading in ("## Goal", "## Setup", "## Inputs and integrity", "## Cell QC", "## Outputs and checks", "## Next steps"):
            assert heading in text, (path, heading)
        assert "scWAT" not in text
        assert "adipose_analysis" not in text
        assert all(cell.get("execution_count") is None for cell in notebook["cells"] if cell["cell_type"] == "code")
        assert all(cell.get("outputs", []) == [] for cell in notebook["cells"] if cell["cell_type"] == "code")


def test_summary_notebook():
    path = NOTEBOOKS / "02_slide_QC_summary.ipynb"
    notebook = json.loads(path.read_text(encoding="utf-8"))
    text = source_text(notebook)
    assert "expected_colon_regions()" in text
    assert "summarise_colon_mouse_position" in text
    assert "Mouse_1" in text and "Mouse_2" in text
    assert "top" in text and "middle" in text and "bottom" in text
    assert "scWAT" not in text and "adipose_analysis" not in text
    assert all(cell.get("execution_count") is None for cell in notebook["cells"] if cell["cell_type"] == "code")
    assert all(cell.get("outputs", []) == [] for cell in notebook["cells"] if cell["cell_type"] == "code")


if __name__ == "__main__":
    test_region_notebooks()
    test_summary_notebook()
    print("Seven colon notebook contracts passed.")
