# python

The Python environment of the repository, managed with
[uv](https://docs.astral.sh/uv/): `pyproject.toml` declares the
dependencies, `uv.lock` pins them. It serves the pipeline helpers in this
folder and the Python analyses in `analyses/`.

## Setup

```sh
uv sync --project python            # creates python/.venv from uv.lock
uv run --project python python analyses/<name>/<script>.py
```

Add a dependency with `uv add --project python <package>` and commit both
`pyproject.toml` and `uv.lock`. Check code with
`uv run --project python ruff check python analyses`.

## Modules

- `utils/log.py`: `configure_pipeline_log(component)` sends a script's log
  records to the pipeline log (`METABO_FIGURES_LOG`, same line format as
  `R/utils-log.R`). Call it in the `if __name__ == '__main__':` block of
  every Python step run by the pipeline.

Coding style: `docs/CODING_STYLE.md`.
