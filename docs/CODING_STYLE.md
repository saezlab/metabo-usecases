# Coding style

The repository follows the Saez lab guidelines in
[saezverse](https://github.com/saezlab/saezverse):
`guidelines/r-coding-style.md` for R and
`guidelines/python-coding-style.md` for Python. This file summarises them
and adds the rules specific to this repository. Where they disagree, the
saezverse guidelines win; update this file.

## Both languages

- Paths are resolved from the repository root: `here::here()` in R
  analyses, `Path(__file__)` in Python. No personal absolute paths, no
  dependence on the working directory.
- Names describe content (`resource_overlap`, not `fig2_panel_c`). Figure
  and table numbers live only in `manuscript.yaml`.
- Analyses write to `data/derived/<name>/`, figures and tables to their
  `out/`. Nothing writes into `data/raw/` or `data/vendored/`.
- Read OmniPath data through the backends in `docs/DB_ACCESS.md`; open
  DuckDB files read-only.
- Comments explain why, above the code they refer to. Notebooks explain
  each decision in text next to the code and plots that support it.
- Prefer the simplest standard method that answers the question.

## R

- Style as in `r-coding-style.md`; `.lintr` checks the essentials:
  snake_case, 80 characters per line, no absolute paths, no `T`/`F`,
  `logger` instead of `print`/`cat`/`message` in the package, `purrr`
  instead of `sapply`/`vapply`, `dplyr::if_else` instead of `ifelse`.
- Functions of the package `metabo.figures` (`R/`) have roxygen2
  documentation; regenerate `man/` and `NAMESPACE` with
  `devtools::document()` when possible. If roxygen2 is not available, edit
  them by hand in roxygen's format.
- Analyses are R Markdown notebooks, run by the folder's entry script in
  fresh sessions; shared helpers go to `analyses/<name>/R/`.

## Python

From `python-coding-style.md` (PEP 8 and the Google style guide, with these
differences):

- Names: `snake_case` functions and variables, `PascalCase` classes.
  Resource names are written as one word (`proteinatlas_annotations`).
  Leave out redundant prefixes (`get_`, `load_`, `parse_`, `create_`):
  `uniprot_locations()`, not `get_uniprot_locations()`.
- Spaces around every operator, also around `=` in keyword arguments and
  defaults: `function(a = 1)`.
- Single quotes, unless there is a reason for double quotes.
- Generous blank lines: two around functions and classes, one before and
  after blocks, and between short groups of related lines inside blocks.
- Long argument lists: one argument per line, one level deeper, trailing
  comma, closing parenthesis on its own line. Long comprehensions: value,
  each `for` and each `if` on their own lines.
- Docstrings: Google (Napoleon) style; the text starts on the line after
  the opening quotes, closing quotes on their own line. Types go into type
  hints, not into the docstring.
- Imports in groups separated by a blank line: standard library, typing,
  third party, this repository. Import heavy dependencies inside the
  function that needs them.
- Return the data structure the caller needs: `dict` for lookups,
  `list[tuple]` for pairs, a data frame for tables.
- Python 3.11+, type hints on public functions.

Specific to this repository:

- One environment for the whole repository: `python/pyproject.toml` and
  `python/uv.lock`. Run with `uv run --project python ...`; add packages
  with `uv add --project python ...`. Do not install into a global Python.
- Lint with `uv run --project python ruff check python analyses`. Do not
  run `ruff format` or black: they rewrite the style above.
- Scripts that the pipeline runs call
  `utils.log.configure_pipeline_log('<component>')` first in their
  `if __name__ == '__main__':` block.
- Notebooks (`.ipynb`) are fine for exploration and for presenting an
  analysis. Code that writes to `data/derived/` lives in a `.py` module
  next to the notebook, which the notebook imports, so it can run without
  Jupyter. Commit notebooks with their outputs cleared.
- The OmniPath clients (`omnipath_client`) are installed from the monorepo
  checkout, not from PyPI (`docs/DB_ACCESS.md`).
