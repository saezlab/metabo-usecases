"""Presentation-ready summary of the core-network-size / solver-time
parameter testing done in scripts 02-04, for Christina's requested
"what was tried and why, what worked, what's the final decision" slides.

This is a hand-curated record of actual runs (not re-derived from result/
files, several of which were overwritten before the labeled-output-folder
convention existed) -- see conversation history / git log for the
underlying logs this summarizes.

Outputs (corneto_starvation/result/QC/):
  parameter_testing_network_size.png     lambda_reg sweep: size vs convergence vs pathway completeness
  parameter_testing_temporal_fit.png     time-budget sweep: convergence outcome per attempt
  parameter_testing_decision.png         final decision + rationale
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUTDIR = Path(__file__).resolve().parents[1] / "result" / "QC"

STATUS_COLOR = {
    "ok": "#1a7a4c",
    "partial": "#b8860b",
    "fail": "#a73935",
    "na": "#8a8f98",
}


def draw_table(
    title: str,
    subtitle: str,
    headers: list[str],
    rows: list[list[str]],
    status_col: int,
    status_map: dict[int, str],
    outpath: Path,
    col_widths: list[float] | None = None,
    figsize: tuple[float, float] = (11, 5),
) -> None:
    fig, ax = plt.subplots(figsize=figsize)
    ax.axis("off")

    n_cols = len(headers)
    n_rows = len(rows)
    col_widths = col_widths or [1.0 / n_cols] * n_cols

    ax.text(0, 1.08, title, fontsize=15, fontweight="bold", transform=ax.transAxes, va="bottom")
    ax.text(0, 1.015, subtitle, fontsize=9.5, color="#555555", transform=ax.transAxes, va="bottom")

    table = ax.table(
        cellText=rows,
        colLabels=headers,
        cellLoc="left",
        colWidths=col_widths,
        bbox=[0, 0, 1, 1],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9.5)

    for (r, c), cell in table.get_celld().items():
        cell.set_edgecolor("#dddddd")
        if r == 0:
            cell.set_facecolor("#eef1ef")
            cell.set_text_props(fontweight="bold")
            cell.set_height(cell.get_height() * 1.3)
        else:
            cell.set_facecolor("#ffffff" if r % 2 else "#f7f8f7")
            if c == status_col:
                status = status_map.get(r - 1, "na")
                cell.set_text_props(color=STATUS_COLOR[status], fontweight="bold")

    fig.tight_layout()
    fig.savefig(outpath, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {outpath}")


def build_network_size_table() -> None:
    headers = ["lambda_reg", "reactions", "solve gap\n(1200s budget)", "status", "key pathways present?"]
    rows = [
        ["0.01", "4,585", "not converged", "partial", "(not individually checked; near-full model)"],
        ["0.10", "4,559", "not converged", "partial", "(not individually checked; near-full model)"],
        ["0.15", "4,466", "not converged", "partial", "(not individually checked; near-full model)"],
        ["0.175", "4,418", "1.33%", "ok", "(not individually checked; near-full model)"],
        ["0.20", "4,267", "1.64%", "ok", "(not individually checked; near-full model)"],
        ["0.50", "1,492", "23.0%", "partial", "glycolysis/PPP/TCA/beta-ox YES — ketogenesis MISSING"],
        ["1.00", "485", "29.5%", "fail", "glycolysis near-absent, PPP/TCA/ketogenesis ALL MISSING"],
    ]
    status_map = {0: "partial", 1: "partial", 2: "partial", 3: "ok", 4: "ok", 5: "partial", 6: "fail"}
    draw_table(
        "Core network extraction — lambda_reg sweep",
        "Script 02 (MultiSampleIMAT, sparsity penalty on reaction union) · full Mouse-GEM = 12,987 reactions · Gurobi, 1200s budget per run",
        headers, rows, status_col=3, status_map=status_map,
        outpath=OUTDIR / "parameter_testing_network_size.png",
        col_widths=[0.12, 0.12, 0.16, 0.12, 0.48],
        figsize=(12, 4.2),
    )


def build_temporal_fit_table() -> None:
    headers = ["network", "genotype(s)", "time budget", "outcome", "result"]
    rows = [
        ["1,492 reactions", "WT only", "300s (5 min)", "fail", "mass-balance check failed on time-limit incumbent"],
        ["1,492 reactions", "WT only", "3,000s (50 min)", "ok", "valid solution — gap 25.5%, 434 switches"],
        ["1,492 reactions", "WT only", "10,800s (3h)", "fail", "flux-bound check failed, despite gap down to ~11%"],
        ["1,492 reactions", "ob/ob only", "3,000s (50 min)", "ok", "valid solution — gap 2.5%, 65 switches"],
        ["485 reactions", "WT + ob/ob joint", "60s (!)", "ok", "proven OPTIMAL — gap 0.01%, 39 / 12 switches"],
    ]
    status_map = {0: "fail", 1: "ok", 2: "fail", 3: "ok", 4: "ok"}
    draw_table(
        "Temporal joint iMAT fit — convergence by network size & time budget",
        "Scripts 03/04 (TemporalMultiSampleIMAT) · same hyperparameters throughout (epsilon_fit=0.03, lambda_reg=0.08, lambda_temporal=0.16, lambda_flux=0.05)",
        headers, rows, status_col=3, status_map=status_map,
        outpath=OUTDIR / "parameter_testing_temporal_fit.png",
        col_widths=[0.17, 0.17, 0.15, 0.1, 0.41],
        figsize=(12, 3.4),
    )


def build_decision_summary() -> None:
    fig, ax = plt.subplots(figsize=(11, 6.5))
    ax.axis("off")

    ax.text(0, 1.0, "Final decision", fontsize=16, fontweight="bold", transform=ax.transAxes, va="top")

    body = (
        "The trade-off is structural, not a matter of more compute:\n\n"
        "  •  Networks large enough to keep every relevant pathway (glycolysis, PPP, TCA,\n"
        "     beta-oxidation, ketogenesis) take hours to converge, sometimes fail validation\n"
        "     even after hours, and WT consistently converges worse than ob/ob at the same size.\n\n"
        "  •  Networks small enough to converge cleanly and fast (485 reactions, 60s, proven\n"
        "     optimal) do so by discarding whole pathways the biology needs — ketogenesis and\n"
        "     most of glycolysis were gone entirely at that size.\n\n"
        "  •  lambda_reg sparsity has no way to know which reactions matter biologically — it\n"
        "     only minimizes count. There is no lambda_reg value that is reliably both small\n"
        "     and complete; it depends on which reactions the sparsity penalty happens to keep.\n\n"
        "Decision: stop tuning lambda_reg on the full Mouse-GEM. Switch to Pablo's\n"
        "metabolomics + proteomics-anchored network construction (MultiSampleFBA on the full\n"
        "model, constrained by real qualitative exchange directions from our own metabolomics\n"
        "data — glucose/lactate/ketone body trends already computed) — this guarantees by\n"
        "construction that the metabolites and pathways we care about stay in the network,\n"
        "rather than hoping a sparsity penalty leaves them in."
    )
    ax.text(0, 0.88, body, fontsize=11, transform=ax.transAxes, va="top", linespacing=1.6, family="monospace")

    fig.tight_layout()
    fig.savefig(OUTDIR / "parameter_testing_decision.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUTDIR / 'parameter_testing_decision.png'}")


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    build_network_size_table()
    build_temporal_fit_table()
    build_decision_summary()


if __name__ == "__main__":
    main()
