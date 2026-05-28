"""SVG-level panel composition for the metabo.figures pipeline.

Invoked from ``R/compose/mixed_source.R`` when SVG-native placement is
required (mostly: composites mixing pre-rendered SVG panels with
manual SVG assets where exact alignment matters).

Composition spec is YAML, passed on stdin or as a file argument:

.. code-block:: yaml

    output: figures/figXX/out/composite.svg
    page:
        width_mm: 180
        height_mm: 220
    panels:
        - id: A
          path: figures/figXX/out/panelA.svg
          x_mm: 0
          y_mm: 0
          width_mm: 90
        - id: B
          path: figures/figXX/out/panelB.svg
          x_mm: 90
          y_mm: 0
          width_mm: 90
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

import svgutils.transform as sg  # type: ignore[import-untyped]
import yaml

from python.utils.log import configure_pipeline_log

# 1 mm = 3.7795275591 px at 96 DPI; SVG default unit is the user unit ==
# px when no viewport scaling is set. We standardize on px throughout.
PX_PER_MM = 96.0 / 25.4


def _to_px(value_mm: float) -> float:
    return value_mm * PX_PER_MM


def assemble(spec: dict[str, Any]) -> Path:
    """Assemble one composite SVG from the spec dict."""

    out_path = Path(spec["output"])
    page = spec["page"]
    panels = spec["panels"]

    page_w = _to_px(page["width_mm"])
    page_h = _to_px(page["height_mm"])

    figure = sg.SVGFigure()
    figure.set_size((f"{page_w}px", f"{page_h}px"))

    elements = []
    for panel in panels:
        panel_path = Path(panel["path"])
        if not panel_path.exists():
            raise FileNotFoundError(
                f"Panel SVG missing: {panel_path}"
            )

        svg = sg.fromfile(str(panel_path))
        root = svg.getroot()

        # Move into the requested (x, y) in px.
        root.moveto(_to_px(panel["x_mm"]), _to_px(panel["y_mm"]))

        # Optional uniform scale to the requested width.
        if "width_mm" in panel:
            target_w = _to_px(panel["width_mm"])
            current_w_str = svg.width if svg.width else None
            if current_w_str:
                current_w = float(
                    str(current_w_str).rstrip("pxptmmcin")
                )
                scale = target_w / current_w
                root.scale(scale)

        elements.append(root)

    figure.append(elements)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    figure.save(str(out_path))

    return out_path


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__.splitlines()[0]
    )
    parser.add_argument(
        "spec",
        nargs="?",
        help="Path to YAML composition spec; '-' or omitted reads from stdin",
    )
    parser.add_argument(
        "--component",
        default="svg_assemble",
        help="Component name for unified-log lines",
    )
    args = parser.parse_args()

    logger = configure_pipeline_log(args.component)

    if args.spec in (None, "-"):
        spec = yaml.safe_load(sys.stdin)
    else:
        with open(args.spec) as fh:
            spec = yaml.safe_load(fh)

    logger.info(
        "Composing %d panels → %s", len(spec["panels"]), spec["output"]
    )
    out = assemble(spec)
    logger.info("Wrote %s", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
