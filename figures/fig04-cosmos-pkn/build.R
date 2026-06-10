#' Build script for Figure 4 — COSMOS PKN comparison
#'
#' Produces panelB (edge-count grouped bar) and panelC (proportional
#' stacked bar) via \code{\link{cosmos_old_vs_new}}. Manual schematics
#' (panels A, D, E) must exist under \code{manual/} before the
#' xelatex composite step runs.
#'
#' @return Invisible NULL; outputs written to \code{out/}.
#' @importFrom logger log_info log_error
#' @importFrom here here
#' @importFrom ggplot2 ggsave
#' @export
build_fig04 <- function() {

    out_dir <- here::here("figures", "fig04-cosmos-pkn", "out")
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

    logger::log_info("[fig04] loading old COSMOS PKN")
    old_pkn <- cosmos_old_pkn()

    logger::log_info("[fig04] loading new OmniPath Metabo PKN edge counts")
    conn <- connect_deployment("dev3")
    on.exit(DBI::dbDisconnect(conn), add = TRUE)
    new_pkn <- cosmos_new_pkn_counts(conn)

    logger::log_info("[fig04] rendering panelB (grouped bar)")
    panels <- cosmos_old_vs_new(old_pkn, new_pkn)

    panel_b_pdf <- file.path(out_dir, "panelB.pdf")
    panel_b_svg <- file.path(out_dir, "panelB.svg")
    ggplot2::ggsave(panel_b_pdf, panels$panel_b,
        width = 89, height = 70, units = "mm", device = "pdf")
    ggplot2::ggsave(panel_b_svg, panels$panel_b,
        width = 89, height = 70, units = "mm", device = "svg")

    logger::log_info("[fig04] rendering panelC (proportional stacked bar)")
    panel_c_pdf <- file.path(out_dir, "panelC.pdf")
    panel_c_svg <- file.path(out_dir, "panelC.svg")
    ggplot2::ggsave(panel_c_pdf, panels$panel_c,
        width = 89, height = 70, units = "mm", device = "pdf")
    ggplot2::ggsave(panel_c_svg, panels$panel_c,
        width = 89, height = 70, units = "mm", device = "svg")

    logger::log_info("[fig04] writing provenance sidecar")
    sidecar_path <- file.path(out_dir, "fig04-cosmos-pkn.pdf.provenance.json")
    write_sidecar(
        artifact_id  = "fig04-cosmos-pkn",
        script_path  = "figures/fig04-cosmos-pkn/build.R",
        deployments  = list(list(name = "dev3")),
        parameters   = list(
            old_pkn_source = attr(old_pkn, "source_path"),
            old_pkn_md5    = attr(old_pkn, "fingerprint")
        ),
        sidecar_path = sidecar_path
    )

    logger::log_info("[fig04] done — outputs in {out_dir}")
    invisible(NULL)
}
