#!/usr/bin/env python3
"""Validate colon notebook JSON structure, with optional Jupyter nbformat validation."""

import argparse
import json
import re
from pathlib import Path

CELL_ID = re.compile(r"^[A-Za-z0-9_-]{1,64}$")


def validate_structure(path: Path) -> None:
    with path.open("r", encoding="utf-8") as handle:
        notebook = json.load(handle)
    if notebook.get("nbformat") != 4:
        raise ValueError(f"{path}: nbformat must be 4")
    if not isinstance(notebook.get("metadata"), dict):
        raise TypeError(f"{path}: notebook metadata must be a JSON object")
    cells = notebook.get("cells")
    if not isinstance(cells, list) or not cells:
        raise TypeError(f"{path}: cells must be a non-empty JSON array")
    cell_ids = []
    for index, cell in enumerate(cells, start=1):
        if not isinstance(cell, dict):
            raise TypeError(f"{path}: cell {index} must be a JSON object")
        if not isinstance(cell.get("metadata"), dict):
            raise TypeError(f"{path}: cell {index} metadata must be a JSON object, not a list")
        cell_id = cell.get("id")
        if not isinstance(cell_id, str) or not CELL_ID.fullmatch(cell_id):
            raise ValueError(f"{path}: cell {index} has an invalid or missing nbformat 4.5 id")
        cell_ids.append(cell_id)
        if cell.get("cell_type") == "code" and not isinstance(cell.get("outputs"), list):
            raise TypeError(f"{path}: code cell {index} outputs must be a JSON array")
        if not isinstance(cell.get("source"), (str, list)):
            raise TypeError(f"{path}: cell {index} source must be a string or string array")
    if len(cell_ids) != len(set(cell_ids)):
        raise ValueError(f"{path}: cell ids must be unique within the notebook")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--require-nbformat", action="store_true")
    parser.add_argument("notebooks", nargs="+", type=Path)
    args = parser.parse_args()
    for path in args.notebooks:
        validate_structure(path)
    try:
        import nbformat  # type: ignore
    except ImportError:
        if args.require_nbformat:
            raise SystemExit("Python package 'nbformat' is required for HPC/Jupyter schema validation.")
        schema_status = "SKIP_NOT_INSTALLED"
    else:
        for path in args.notebooks:
            notebook = nbformat.read(path, as_version=4)
            nbformat.validate(notebook)
        schema_status = "PASS"
    print(f"NOTEBOOKS_VALID={len(args.notebooks)} STRUCTURE=PASS NBFORMAT_SCHEMA={schema_status}")


if __name__ == "__main__":
    main()
