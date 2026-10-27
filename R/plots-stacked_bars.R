#' Canonical COSMOS-PKN resource order for Figure 6 Panels E / F
#'
#' Reproduces the resource ordering of the legacy
#' \code{analyses/cancer-cell-lines-cosmos/03_visualization.R} so the
#' stacked-bar fills render in a publication-stable sequence.
#'
#' @return Character vector of resource names in canonical order.
#'
#' @keywords internal
#' @noRd
case_study_resource_order <- function() {
    c(
        "STITCH",
        "BRENDA",
        "GEM:Human-GEM",
        "GEM:Recon3D",
        "GEM:Human-GEM;GEM:Recon3D",
        "KEGG",
        "GEM:Human-GEM;KEGG",
        "GEM:Recon3D;KEGG",
        "GEM:Human-GEM;GEM:Recon3D;KEGG"
    )
}


#' Strip the redundant "GEM:" prefix for display
#'
#' The underlying PKN edges carry the prefixed labels
#' (\code{"GEM:Human-GEM"}, \code{"GEM:Recon3D"}, plus composite
#' memberships like \code{"GEM:Human-GEM;KEGG"}), so the data column
#' and the category-colour registry keys stay prefixed for
#' downstream consumers. Legends and axis ticks display the cleaner
#' \code{"Human-GEM"} / \code{"Recon3D"} forms via this helper.
#'
#' @param x Character vector.
#' @return Character vector with every \code{"GEM:"} occurrence
#'     removed.
#'
#' @keywords internal
#' @noRd
case_study_resource_display <- function(x) {
    gsub("GEM:", "", x, fixed = TRUE)
}


#' Canonical resource → hex palette for Panels E / F
#'
#' Returns a named character vector mapping every resource in the
#' canonical order to its registered colour. Used in
#' \code{scale_fill_manual()} so the same scale object appears in
#' every E / F sub-panel — that's what lets
#' \code{patchwork::plot_layout(guides = "collect")} merge the
#' Resource legend reliably (\code{scale_fill_identity} does not
#' collect across plots even with identical breaks / labels).
#'
#' @return Named character vector: names are canonical resource
#'     names (in canonical order), values are hex codes from the
#'     \code{"resources"} category registry.
#'
#' @keywords internal
#' @noRd
case_study_resource_palette <- function() {
    resources <- case_study_resource_order()
    stats::setNames(
        unname(category_colour("resources", resources)),
        resources
    )
}


#' Single-letter to full subcellular-location labels used in COSMOS PKN
#'
#' @return Named character vector keyed by single-letter code.
#'
#' @keywords internal
#' @noRd
case_study_location_labels <- function() {
    c(
        c = "Cytoplasm",
        m = "Mitochondria",
        n = "Nucleus",
        e = "Extracellular",
        r = "ER",
        x = "Peroxisome",
        g = "Golgi",
        v = "Vesicle",
        l = "Lysosome"
    )
}


#' Register the COSMOS-PKN resource palette in the category registry
#'
#' Looks up the canonical resource order and assigns the first N
#' colours of the lead palette (FR-020 / FR-022 / SC-004), then
#' registers them under the \code{"resources"} category so Panels E
#' and F render the same resource in the same colour across both
#' plots, and across any other figure that touches the COSMOS PKN.
#'
#' @return Invisibly the registered mapping (named character vector).
#'
#' @examples
#' \dontrun{
#' register_case_study_resource_colours()
#' category_colour("resources", "GEM:Human-GEM")
#' }
#'
#' @importFrom stats setNames
#' @export
register_case_study_resource_colours <- function() {

    resources <- case_study_resource_order()
    fills <- stats::setNames(
        palette_n(as.integer(length(resources)), unknown = FALSE),
        resources
    )
    register_category_colours("resources", fills)
    invisible(fills)
}


#' Filter a COSMOS-PKN tibble to edges touching a DEM ChEBI ID
#'
#' Reproduces the legacy \code{screen_pkn()} helper from
#' \code{02_connect_dem_pkn.R}: keeps PKN rows where either
#' \code{source} or \code{target} matches a ChEBI in
#' \code{chebi_ids}. Adds a \code{matched_chebi} column carrying the
#' matched ID, used by the location-unnest step in Panel F.
#'
#' @param pkn Tibble: a single COSMOS-PKN slice (allosteric or
#'     enzyme-metabolite).
#' @param chebi_ids Character vector of unique trimmed ChEBI IDs
#'     produced by \code{\link{case_study_extract_chebi}}.
#'
#' @return Tibble: the PKN rows whose source or target matches a
#'     ChEBI, with an added \code{matched_chebi} column.
#'
#' @importFrom dplyr filter mutate case_when
#' @importFrom rlang .data
#' @export
case_study_screen_pkn <- function(pkn, chebi_ids) {

    # NSE vs. R CMD check workaround
    source <- target <- NULL

    pkn |>
        dplyr::filter(
            .data$source %in% chebi_ids |
                .data$target %in% chebi_ids
        ) |>
        dplyr::mutate(
            matched_chebi = dplyr::case_when(
                .data$source %in% chebi_ids ~ .data$source,
                .data$target %in% chebi_ids ~ .data$target
            )
        )
}


#' Summarise PKN edges by direction (up / down) and PKN slice
#'
#' Builds the four-group long table consumed by
#' \code{\link{gem_allosteric_panel}}. Groups are
#' \code{"up - Allos"}, \code{"down - Allos"}, \code{"up - GEM"},
#' \code{"down - GEM"} — preserving the legacy
#' \code{03_visualization.R} naming.
#'
#' @param dem_tibble Output of \code{\link{case_study_differential}}
#'     for a single contrast; must carry \code{logFC} and
#'     \code{chebi}.
#' @param pkn Output of \code{\link{case_study_cosmos_pkn}}: a named
#'     list with \code{allosteric} and \code{enzyme_metabolite}
#'     tibbles.
#' @param top_n Integer: top-N DEMs per direction (legacy default
#'     \code{10L}). Selection is by t-statistic.
#'
#' @return A tibble combining matched edges across all four groups
#'     with columns \code{group} (factor with levels in legacy
#'     order), \code{resource} (factor in canonical order), plus the
#'     original PKN columns (including \code{locations}).
#'
#' @importFrom dplyr filter mutate bind_rows
#' @importFrom rlang abort
#' @export
case_study_pkn_summary <- function(dem_tibble, pkn, top_n = 10L) {

    # NSE vs. R CMD check workaround
    logFC <- NULL

    if (!all(c("allosteric", "enzyme_metabolite") %in% names(pkn))) {
        rlang::abort(paste0(
            "`pkn` must be the output of case_study_cosmos_pkn() ",
            "(a named list with 'allosteric' and 'enzyme_metabolite')."
        ))
    }

    top <- case_study_top_dems(dem_tibble, n = top_n)
    up   <- dplyr::filter(top, logFC > 0)
    down <- dplyr::filter(top, logFC < 0)

    up_chebi   <- case_study_extract_chebi(up)
    down_chebi <- case_study_extract_chebi(down)

    group_levels <- c(
        "up - Allos", "down - Allos", "up - GEM", "down - GEM"
    )

    rows <- dplyr::bind_rows(
        dplyr::mutate(
            case_study_screen_pkn(pkn$allosteric, up_chebi),
            group = "up - Allos"
        ),
        dplyr::mutate(
            case_study_screen_pkn(pkn$allosteric, down_chebi),
            group = "down - Allos"
        ),
        dplyr::mutate(
            case_study_screen_pkn(pkn$enzyme_metabolite, up_chebi),
            group = "up - GEM"
        ),
        dplyr::mutate(
            case_study_screen_pkn(pkn$enzyme_metabolite, down_chebi),
            group = "down - GEM"
        )
    )

    dplyr::mutate(
        rows,
        group    = factor(rows$group, levels = group_levels),
        resource = factor(
            rows$resource, levels = case_study_resource_order()
        )
    )
}


#' Figure 6 Panel E -- GEM / allosteric edge counts by direction
#'
#' Stacked-bar plot: x = group (up / down × Allos / GEM), y = edge
#' count, fill = resource. Reproduces the legacy
#' \code{03_visualization.R} \code{bar_group()} composition with the
#' shared style module applied.
#'
#' @param summary_tibble Output of
#'     \code{\link{case_study_pkn_summary}}.
#' @param contrast_label Character: panel title (e.g.
#'     \code{"KRAS"} or \code{"EGFR"}).
#' @param width_mm Numeric: target physical panel width in mm.
#' @param font_scale Numeric: passed to \code{\link{theme_bw_metabo}}.
#' @param legend_scale Numeric or \code{NULL}: passed to
#'     \code{\link{theme_bw_metabo}}.
#'
#' @return A ggplot object.
#'
#' @importFrom dplyr count
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_manual labs
#' @importFrom ggplot2 theme element_text
#' @importFrom rlang .data
#' @export
gem_allosteric_panel <- function(
    summary_tibble,
    contrast_label,
    width_mm = 89L,
    font_scale = 1,
    legend_scale = NULL
) {

    # NSE vs. R CMD check workaround
    group <- resource <- n <- NULL

    counts <- dplyr::count(
        summary_tibble,
        .data$group,
        .data$resource,
        name   = "n",
        .drop  = FALSE
    )

    resources <- case_study_resource_order()

    ggplot2::ggplot(
        counts,
        ggplot2::aes(
            x    = .data$group,
            y    = .data$n,
            fill = .data$resource
        )
    ) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::scale_fill_manual(
            values = case_study_resource_palette(),
            breaks = resources,
            labels = case_study_resource_display(resources),
            name   = "Resource",
            drop   = FALSE
        ) +
        ggplot2::labs(
            title = contrast_label,
            x     = NULL,
            y     = "Number of edges"
        ) +
        theme_bw_metabo(
            width_mm     = width_mm,
            font_scale   = font_scale,
            legend_scale = legend_scale
        ) +
        ggplot2::theme(
            axis.text.x     = ggplot2::element_text(
                angle = 30, hjust = 1
            ),
            legend.position = "bottom"
        )
}


#' Figure 6 Panel F -- subcellular-location stacked bars
#'
#' Builds one stacked-bar panel per contrast × direction (KRAS up,
#' KRAS down, EGFR up, EGFR down). Each panel groups edges by
#' single-letter subcellular-location code (expanded to full names
#' on the x axis) and stacks by resource.
#'
#' Compartment letters are parsed from the COSMOS PKN
#' \code{locations} column, which carries a Python tuple
#' representation like \code{"('c', 'e', 'n')"}; the legacy
#' \code{03_visualization.R} regex \code{'([a-z])'} is reproduced.
#' Edges with an empty \code{locations} tuple contribute no rows
#' (legacy behaviour).
#'
#' @param summary_tibble Output of
#'     \code{\link{case_study_pkn_summary}}; must carry the
#'     \code{locations} column from the PKN.
#' @param contrast_label Character: panel title (e.g.
#'     \code{"KRAS"} or \code{"EGFR"}).
#' @param direction Character: one of \code{"up"} or \code{"down"};
#'     filters \code{summary_tibble} to edges with that direction.
#' @param width_mm Numeric: target physical panel width in mm.
#' @param font_scale Numeric: passed to \code{\link{theme_bw_metabo}}.
#' @param legend_scale Numeric or \code{NULL}: passed to
#'     \code{\link{theme_bw_metabo}}.
#'
#' @return A ggplot object.
#'
#' @importFrom dplyr filter mutate count
#' @importFrom tidyr unnest
#' @importFrom stringr str_match_all str_detect
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_manual
#' @importFrom ggplot2 scale_x_discrete labs theme element_text
#' @importFrom rlang .data abort
#' @export
subcellular_location_panel <- function(
    summary_tibble,
    contrast_label,
    direction = c("up", "down"),
    width_mm = 89L,
    font_scale = 1,
    legend_scale = NULL
) {

    # NSE vs. R CMD check workaround
    group <- locations <- loc <- resource <- n <- NULL

    direction <- match.arg(direction)

    if (!"locations" %in% names(summary_tibble)) {
        rlang::abort(
            "summary_tibble is missing the 'locations' column"
        )
    }

    sub <- dplyr::filter(
        summary_tibble,
        stringr::str_detect(as.character(.data$group), direction)
    )

    sub <- dplyr::mutate(
        sub,
        loc = lapply(
            stringr::str_match_all(.data$locations, "'([a-z])'"),
            function(m) if (nrow(m) == 0L) character(0L) else m[, 2L]
        )
    )
    sub <- tidyr::unnest(sub, cols = "loc")

    # Cast loc to a factor with every canonical level so each F
    # sub-panel renders an identical x-axis. Combined with
    # scale_x_discrete(drop = FALSE), localizations absent from a
    # given (contrast, direction) slice show as empty slots so bar
    # widths stay constant for cross-sub-panel comparison.
    loc_levels <- names(case_study_location_labels())
    sub$loc <- factor(sub$loc, levels = loc_levels)

    counts <- dplyr::count(
        sub,
        .data$loc,
        .data$resource,
        name  = "n",
        .drop = FALSE
    )

    resources <- case_study_resource_order()

    ggplot2::ggplot(
        counts,
        ggplot2::aes(
            x    = .data$loc,
            y    = .data$n,
            fill = .data$resource
        )
    ) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::scale_x_discrete(
            labels = case_study_location_labels(),
            drop   = FALSE
        ) +
        ggplot2::scale_fill_manual(
            values = case_study_resource_palette(),
            breaks = resources,
            labels = case_study_resource_display(resources),
            name   = "Resource",
            drop   = FALSE
        ) +
        ggplot2::labs(
            title = sprintf("%s %s", contrast_label, direction),
            x     = NULL,
            y     = "Number of edges"
        ) +
        theme_bw_metabo(
            width_mm     = width_mm,
            font_scale   = font_scale,
            legend_scale = legend_scale
        ) +
        ggplot2::theme(
            axis.text.x     = ggplot2::element_text(
                angle = 30, hjust = 1
            ),
            legend.position = "bottom"
        )
}
