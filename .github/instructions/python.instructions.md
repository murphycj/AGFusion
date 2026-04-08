---
description: "Use when writing, editing, or reviewing Python code in the AGFusion package. Covers code style, module responsibilities, testing patterns, and Python-specific conventions."
applyTo: "**/*.py"
---

# AGFusion Python Guidelines

## Code Style

- Format all Python code with **black** (`line-length = 100`)
- Target Python 3.7+
- Run `black .` before committing

## Module Responsibilities

| Module                   | Purpose                                                                              |
| ------------------------ | ------------------------------------------------------------------------------------ |
| `agfusion/model.py`      | Core domain classes: `_Gene`, `Fusion`. Holds exon/protein architecture logic.       |
| `agfusion/database.py`   | `AGFusionDB` — SQLite-backed store of protein domain annotations (Pfam, TMHMM, etc.) |
| `agfusion/parsers.py`    | Parsers for fusion-finding algorithm output files (one parser per algorithm)         |
| `agfusion/plot.py`       | Matplotlib-based visualization of protein domain and exon structures                 |
| `agfusion/cli.py`        | Click-based CLI entry points: `annotate`, `batch`, `download`                        |
| `agfusion/utils.py`      | Shared utility functions                                                             |
| `agfusion/exceptions.py` | Custom exception classes                                                             |

## Python Conventions

- Protein database annotation keys: `"pfam"`, `"tmhmm"`, etc.
- `pyensembl.EnsemblRelease` objects provide transcript/exon lookups
- Add new fusion-finding parsers to `parsers.py` — one parser per algorithm

## Testing

- Tests use `unittest.TestCase` in `test/`
- Tests require local AGFusion `.db` files and pyensembl data installed
- Add new test cases to `test_base.py`, `test_parsers.py`, or `test_plots.py` as appropriate
- Run tests with: `python -m pytest test/`
