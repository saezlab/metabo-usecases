library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(patchwork)

DATA_FILE <- "/Users/priscilla/Library/CloudStorage/OneDrive-UniversitätHeidelberg/00_project/CCLE_case_study/data/ExtendedDataTable_DifferentialAnalysis_Shorthouse_LungMutations.xlsx"
PKN_DIR   <- "/Users/priscilla/Library/CloudStorage/OneDrive-UniversitätHeidelberg/00_project/CCLE_case_study/data"
FIG_DIR   <- "/Users/priscilla/Library/CloudStorage/OneDrive-UniversitätHeidelberg/00_project/CCLE_case_study/figure"

dir.create(FIG_DIR, showWarnings = FALSE)

# ── Colour palette (in canonical order) ───────────────────────────────────────

palette_order <- c("#006384", "#9F0162", "#FEAF16", "#BBCC33", "#EA6572",
                   "#009E73", "#99DDFF", "#D03293", "#984EA3")
palette_na    <- "#BEBEBE"

# 9 distinct resources → first 9 colours in order
resource_order <- c("STITCH", "BRENDA",
                    "GEM:Human-GEM", "GEM:Recon3D", "GEM:Human-GEM;GEM:Recon3D",
                    "KEGG",
                    "GEM:Human-GEM;KEGG", "GEM:Recon3D;KEGG",
                    "GEM:Human-GEM;GEM:Recon3D;KEGG")
res_colors <- setNames(palette_order[seq_along(resource_order)], resource_order)

loc_labels <- c(
  c = "cytoplasm", m = "mitochondria", n = "nucleus",
  e = "extracellular", r = "ER", x = "peroxisome",
  g = "Golgi", v = "vesicle", l = "lysosome"
)

# ── Shared theme ──────────────────────────────────────────────────────────────

theme_fig <- function() {
  theme_bw() +
  theme(
    text                 = element_text(family = "Arial"),
    panel.grid.major     = element_blank(),
    panel.grid.minor     = element_blank(),
    axis.text            = element_text(family = "Arial", size = 9),
    axis.title           = element_text(family = "Arial", size = 10),
    legend.title         = element_text(family = "Arial", size = 10, face = "bold"),
    legend.text          = element_text(family = "Arial", size = 10),
    plot.title           = element_text(family = "Arial", size = 10, face = "bold"),
    strip.text           = element_text(family = "Arial", size = 9),
    axis.text.x          = element_text(family = "Arial", size = 9, angle = 30, hjust = 1),
    legend.position      = "bottom",
    legend.key.size      = unit(0.4, "cm")
  )
}

# ── Figure sizing helpers ─────────────────────────────────────────────────────
# y_cm = 8 (fixed plot-area height), x_cm = 0.5 × n_bars (plot-area width)
# margins: left ~2.5 cm (y-axis title + labels), right 0.3 cm,
#           top 0.8 cm (title), bottom 1.5 cm (x labels)

margin_left   <- 2.5
margin_right  <- 0.3
margin_top    <- 0.8
margin_bottom <- 1.5
y_cm          <- 8

fig_dims <- function(n_bars, n_panels_v = 1, legend_h = 1.2) {
  w <- n_bars * 0.5 + margin_left + margin_right
  h <- n_panels_v * (y_cm + margin_top + margin_bottom) + legend_h
  list(w = w / 2.54, h = h / 2.54)   # convert cm → inches for pdf()
}

# For faceted plots: total bars = sum across all facets
fig_dims_facet <- function(bars_per_facet, n_facets, n_panels_v = 1, legend_h = 1.2) {
  w <- (bars_per_facet * 0.5 + 0.6) * n_facets + margin_left + margin_right
  h <- n_panels_v * (y_cm + margin_top + margin_bottom) + legend_h
  list(w = w / 2.54, h = h / 2.54)
}

# ── Data helpers ──────────────────────────────────────────────────────────────

load_dem <- function(sheet) {
  raw <- read_excel(DATA_FILE, sheet = sheet) |>
    select(metabolite_name, logFC, t, chebi)
  bind_rows(
    slice_max(raw, t, n = 10),
    slice_min(raw, t, n = 10)
  )
}

chebi_map <- function(dem_df) {
  dem_df |>
    filter(!is.na(chebi)) |>
    mutate(chebi_single = strsplit(chebi, ";")) |>
    unnest(chebi_single) |>
    mutate(chebi_single = trimws(chebi_single)) |>
    select(chebi_single) |>
    distinct()
}

match_pkn <- function(pkn, cmap) {
  pkn |>
    mutate(matched_chebi = case_when(
      source %in% cmap$chebi_single ~ source,
      target %in% cmap$chebi_single ~ target
    )) |>
    filter(!is.na(matched_chebi))
}

build_summary <- function(up_map, down_map) {
  group_levels <- c("up - Allos", "down - Allos", "up - GEM", "down - GEM")
  bind_rows(
    match_pkn(pkn_allos, up_map)   |> mutate(group = "up - Allos"),
    match_pkn(pkn_allos, down_map) |> mutate(group = "down - Allos"),
    match_pkn(pkn_gem,   up_map)   |> mutate(group = "up - GEM"),
    match_pkn(pkn_gem,   down_map) |> mutate(group = "down - GEM")
  ) |>
    mutate(group = factor(group, levels = group_levels),
           resource = factor(resource, levels = resource_order))
}

count_by_group <- function(df) {
  df |> count(group, resource, name = "n", .drop = FALSE)
}

count_by_location <- function(df, direction) {
  df |>
    filter(str_detect(group, direction)) |>
    mutate(loc = str_match_all(locations, "'([a-z])'") |>
             lapply(\(m) if (nrow(m) == 0) character(0) else m[, 2])) |>
    unnest(loc) |>
    count(loc, resource, name = "n")
}

# ── Plot functions ────────────────────────────────────────────────────────────

bar_group <- function(df, title) {
  ggplot(df, aes(x = group, y = n, fill = resource)) +
    geom_col() +
    scale_fill_manual(values = res_colors, drop = FALSE, name = "Resource") +
    labs(title = title, x = NULL, y = "Number of edges") +
    theme_fig()
}

bar_location <- function(df, title) {
  ggplot(df, aes(x = loc, y = n, fill = resource)) +
    geom_col() +
    scale_x_discrete(labels = loc_labels) +
    scale_fill_manual(values = res_colors, drop = FALSE, name = "Resource") +
    labs(title = title, x = NULL, y = "Number of edges") +
    theme_fig()
}

# ── Load PKN ──────────────────────────────────────────────────────────────────

pkn_allos <- read.csv(file.path(PKN_DIR, "pkn_allosteric.csv"))
pkn_gem   <- read.csv(file.path(PKN_DIR, "pkn_enzyme_metabolite.csv"))

# ── KRAS ─────────────────────────────────────────────────────────────────────

kras_dem  <- load_dem("KRAS_filt_limma")
kras_up   <- chebi_map(filter(kras_dem, logFC > 0))
kras_down <- chebi_map(filter(kras_dem, logFC < 0))
kras_all  <- build_summary(kras_up, kras_down)

# ── EGFR ─────────────────────────────────────────────────────────────────────

egfr_dem  <- load_dem("EGFR_filt_limma")
egfr_up   <- chebi_map(filter(egfr_dem, logFC > 0))
egfr_down <- chebi_map(filter(egfr_dem, logFC < 0))
egfr_all  <- build_summary(egfr_up, egfr_down)

# ── Figure 1: edge count by group ─────────────────────────────────────────────

p_kras <- bar_group(count_by_group(kras_all), "KRAS")
p_egfr <- bar_group(count_by_group(egfr_all), "EGFR")

fig1 <- p_kras / p_egfr +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom") &
  guides(fill = guide_legend(nrow = 2))

dims1 <- fig_dims(n_bars = 4, n_panels_v = 2)
quartz(type = "pdf", file = file.path(FIG_DIR, "fig1_DEM_PKN_groups.pdf"),
       width = 10, height = dims1$h, family = "Arial")
print(fig1)
dev.off()
cat(sprintf("Saved: fig1_DEM_PKN_groups.pdf  (%.1f x %.1f cm)\n",
            10 * 2.54, dims1$h * 2.54))

# ── Figure 2: edge count by location ─────────────────────────────────────────

l_kras_up   <- bar_location(count_by_location(kras_all, "up"),   "KRAS up")
l_kras_down <- bar_location(count_by_location(kras_all, "down"), "KRAS down")
l_egfr_up   <- bar_location(count_by_location(egfr_all, "up"),   "EGFR up")
l_egfr_down <- bar_location(count_by_location(egfr_all, "down"), "EGFR down")

fig2 <- l_kras_up | l_kras_down | l_egfr_up | l_egfr_down
fig2 <- fig2 +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.key.size = unit(0.25, "cm"),
        legend.text     = element_text(family = "Arial", size = 7),
        legend.title    = element_text(family = "Arial", size = 7, face = "bold")) &
  guides(fill = guide_legend(nrow = 1))

max_locs <- bind_rows(
  count_by_location(kras_all, "up"),   count_by_location(kras_all, "down"),
  count_by_location(egfr_all, "up"),   count_by_location(egfr_all, "down")
) |> pull(loc) |> n_distinct()

dims2 <- fig_dims_facet(bars_per_facet = max_locs, n_facets = 4, n_panels_v = 1)
quartz(type = "pdf", file = file.path(FIG_DIR, "fig2_DEM_PKN_locations.pdf"),
       width = dims2$w, height = dims2$h, family = "Arial")
print(fig2)
dev.off()
cat(sprintf("Saved: fig2_DEM_PKN_locations.pdf  (%.1f x %.1f cm)\n",
            dims2$w * 2.54, dims2$h * 2.54))

