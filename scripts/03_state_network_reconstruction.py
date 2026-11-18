"""Reconstruct the active metabolic network at each of the 16 states.

For each (genotype, timepoint) state, take the tasks found feasible by
02_temporal_task_feasibility and the reactions its proteomics closed, then run
the same two-stage solve as the base network (script 01): the smallest
reaction set that supports all the feasible tasks, then a parsimonious flux on
that fixed-size network. The result is a state-specific minimal network that
shows how the active metabolism rewires over the starvation time course --
WT expands as fasting modules switch on, ob/ob stays smaller.

Infeasible tasks are left out: a joint solve that included one would itself be
infeasible, so only the feasible subset is reconstructed.

Reads 02's feasibility_status_matrix.csv (which tasks are feasible per state)
and recomputes each state's reaction closures from the gene scores.

Outputs (result/state_network_reconstruction/):
  selected_reactions_by_state.csv  reaction x state membership in each network
  reconstruction_summary.csv       per state: feasible tasks, reaction count
  build_summary.txt                per-genotype reaction-count trajectory
  figures/<state>.{svg,png}        union network at that state
  figures/<state>/<task>.{svg,png} flux-carrying reactions per task
  figures/<state>_gallery.html     the combined gallery (tutorial's own view)
"""

from __future__ import annotations

import importlib
import re
import subprocess
import sys
from pathlib import Path

PKG_SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(PKG_SRC))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from corneto_mcmt.model_example import plot_task_flux_gallery  # noqa: E402

network = importlib.import_module("01_metabolic_task_network")  # noqa: E402
gating = importlib.import_module("02_temporal_task_feasibility")  # noqa: E402

PKG_ROOT = Path(__file__).resolve().parents[1]
OUTDIR = PKG_ROOT / "result" / "state_network_reconstruction"
FEASIBILITY_MATRIX = PKG_ROOT / "result" / "temporal_task_feasibility" / "feasibility_status_matrix.csv"


def _sanitize(name: str) -> str:
    """Filesystem-safe task name, e.g. 'Gluconeogenesis from lactate' -> 'Gluconeogenesis_from_lactate'."""
    return re.sub(r"[^0-9A-Za-z]+", "_", name).strip("_")


def _split_gallery(html: str) -> list[tuple[str, str]]:
    """Pull (title, svg) pairs from a plot_task_flux_gallery HTML, in order.
    The first is the union network; the rest are one per task."""
    pairs = re.findall(r"<h4[^>]*>(.*?)</h4>.*?(<svg.*?</svg>)", html, re.DOTALL)
    return [(re.sub(r"<[^>]+>", "", title).strip(), svg) for title, svg in pairs]


def _render(svg: str, out_base: Path) -> None:
    """Write an .svg next to a .png rendered with rsvg-convert."""
    out_base.parent.mkdir(parents=True, exist_ok=True)
    if "xmlns" not in svg[:200]:
        svg = svg.replace(
            "<svg ", '<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" ', 1
        )
    svg_path = out_base.with_suffix(".svg")
    svg_path.write_text(svg)
    subprocess.run(
        ["rsvg-convert", "--background-color=white", "-o", str(out_base.with_suffix(".png")), str(svg_path)],
        check=True,
    )


def main() -> None:
    if not FEASIBILITY_MATRIX.exists():
        raise RuntimeError(
            f"{FEASIBILITY_MATRIX} not found -- run 02_temporal_task_feasibility.py first."
        )
    figures_dir = OUTDIR / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    graph = network.build_full_metabolic_network()
    task_bounds_all = network.build_task_bounds()
    task_names = list(task_bounds_all)
    reaction_ids = list(graph.get_attr_from_edges("id"))
    print(f"Network: {graph.num_edges} reactions, {len(task_names)} tasks")

    status_matrix = pd.read_csv(FEASIBILITY_MATRIX, index_col=0)

    network_genes = gating.network_mouse_genes(graph)
    wide_scores = pd.read_csv(gating.GENE_SCORES_PATH, index_col=0)
    low_by_state = gating.low_expression_genes_per_state(wide_scores, network_genes)

    summary_rows = []
    membership = pd.DataFrame(index=reaction_ids, dtype=bool)

    for genotype in gating.GENOTYPES:
        for timepoint in gating.STATE_SUFFIXES:
            state = f"{genotype}_{timepoint}"
            feasible = [t for t in task_names if bool(status_matrix.loc[t, state])]
            if not feasible:
                print(f"  {state}: no feasible task -- skipped")
                continue

            availability = gating.reaction_availability(graph, low_by_state[state])
            closure = {
                rid: (0.0, 0.0) for rid, info in availability.items() if not info["available"]
            }
            sub_bounds = {t: {**task_bounds_all[t], **closure} for t in feasible}
            selected, flux, n_reactions = network.solve_minimum_network(graph, sub_bounds)

            membership[state] = selected
            summary_rows.append({
                "state": state, "genotype": genotype, "timepoint": timepoint,
                "n_feasible_tasks": len(feasible), "n_reactions": n_reactions,
            })

            html = plot_task_flux_gallery(graph, flux, feasible, selected)
            (figures_dir / f"{state}_gallery.html").write_text(f"<html><body>{html}</body></html>")
            panels = _split_gallery(html)
            if panels:
                _render(panels[0][1], figures_dir / state)              # union network
                for title, svg in panels[1:]:
                    _render(svg, figures_dir / state / _sanitize(title))  # per task

            print(f"  {state}: {len(feasible)} feasible tasks -> {n_reactions}/{graph.num_edges} reactions")

    membership.fillna(False).astype(int).to_csv(OUTDIR / "selected_reactions_by_state.csv")
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(OUTDIR / "reconstruction_summary.csv", index=False)

    lines = [f"Network: {graph.num_edges} reactions, {len(task_names)} tasks\n"]
    for genotype in gating.GENOTYPES:
        sub = summary_df[summary_df["genotype"] == genotype]
        traj = ", ".join(f"{r.timepoint}:{r.n_reactions}" for r in sub.itertuples())
        lo, hi = int(sub["n_reactions"].min()), int(sub["n_reactions"].max())
        lines.append(f"{genotype}: {lo}-{hi} of {graph.num_edges} reactions (0h->24h: {traj})")
    summary = "\n".join(lines) + "\n"
    (OUTDIR / "build_summary.txt").write_text(summary)
    print("\n" + summary)


if __name__ == "__main__":
    main()
