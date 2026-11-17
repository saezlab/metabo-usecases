"""CORNETO graph loader for the bundled Mouse-GEM model."""

from __future__ import annotations

from importlib import resources
from pathlib import Path

import corneto as cn
from corneto.io import import_sbml_model


_MOUSE_GEM_FILENAME = "Mouse-GEM-v1.8.0.xml"


def _default_mouse_gem_path() -> Path:
    """Find the checked-out model asset or its copy embedded in an installed wheel."""
    repository_asset = (
        Path(__file__).resolve().parents[2] / "assets" / "models" / _MOUSE_GEM_FILENAME
    )
    if repository_asset.is_file():
        return repository_asset

    packaged_asset = resources.files("corneto_mcmt").joinpath(
        "assets", "models", _MOUSE_GEM_FILENAME
    )
    if not packaged_asset.is_file():
        raise FileNotFoundError(
            f"Bundled Mouse-GEM model not found: {repository_asset}. "
            "Pass the path to a Mouse-GEM SBML file explicitly."
        )
    return Path(str(packaged_asset))


def load_mouse_gem(path: str | Path | None = None) -> cn.Graph:
    """Load Mouse-GEM SBML as a CORNETO metabolic hypergraph.

    Reaction IDs, bounds, net stoichiometry, names, and GPR strings are stored
    on graph edges. Boundary and constant species are retained in graph
    metadata rather than added as FBA balance rows.
    """
    model_path = Path(path) if path is not None else _default_mouse_gem_path()
    if not model_path.is_file():
        raise FileNotFoundError(f"SBML model file does not exist: {model_path}")
    return import_sbml_model(model_path)
