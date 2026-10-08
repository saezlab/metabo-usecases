# User Story 3 capstone: summary statistics and key-driver/key-metabolite
# tables across the 8-timepoint pruned mechanistic network series, for a
# first-pass read of the reconstructed starvation response.
#
# Run from omnipath_metabo_case2/, after scripts/06_footprint_moon.R.

moon_results <- readRDS("result/moon/all_timepoints.rds")
measured_features <- readRDS("result/pk_retrieval/measured_features.rds")

timepoints <- as.numeric(names(moon_results))

## ---------------------------------------------------------------------
## Name lookups (CHEBI -> metabolite name, UniProt -> gene symbol) from
## measured_features$feature_id, which encodes "<layer>;<name>;<id>...".
## ---------------------------------------------------------------------

name_lookup <- function(omics_layers) {
    rows <- measured_features[measured_features$omics_layer %in% omics_layers & !is.na(measured_features$pk_node_id), ]
    parts <- strsplit(rows$feature_id, ";", fixed = TRUE)
    name <- vapply(parts, function(p) if (length(p) >= 2) p[2] else NA_character_, character(1))
    lut <- stats::setNames(name, rows$pk_node_id)
    lut[!duplicated(names(lut))]
}
metab_name <- name_lookup(c("metabolome", "plasma_metabolome"))
gene_symbol <- name_lookup(c("transcriptome", "proteome"))

strip_compartment <- function(metab_node_id) sub("^Metab__", "", sub("_[a-z]+$", "", metab_node_id))

## ---------------------------------------------------------------------
## 1. Network size per timepoint
## ---------------------------------------------------------------------

size_table <- do.call(rbind, lapply(names(moon_results), function(tp) {
    r <- moon_results[[tp]]
    data.frame(
        timepoint_h = as.numeric(tp),
        pruned_nodes = nrow(r$nodes),
        direct_edges = nrow(r$edges),
        gem_edges_reattached = nrow(r$gem_edges),
        metabolite_nodes = sum(grepl("^Metab__", r$nodes$source)),
        tf_kinase_upstream_nodes = sum(r$nodes$type == "upstream_input")
    )
}))
size_table <- size_table[order(size_table$timepoint_h), ]

cat("=== Network size per timepoint ===\n")
print(size_table, row.names = FALSE)
write.csv(size_table, "result/moon/network_size_per_timepoint.csv", row.names = FALSE)

## ---------------------------------------------------------------------
## 2. Key metabolites: trajectory of MOON score across timepoints
## (level0 nodes = directly measured metabolite t-stats, the ones that
## actually "change" in the input data -- not propagated scores).
## ---------------------------------------------------------------------

metab_long <- do.call(rbind, lapply(names(moon_results), function(tp) {
    r <- moon_results[[tp]]
    met <- r$nodes[grepl("^Metab__", r$nodes$source) & r$nodes$type == "level0", ]
    if (nrow(met) == 0) return(NULL)
    data.frame(
        timepoint_h = as.numeric(tp),
        chebi_id = strip_compartment(met$source),
        score = met$score,
        stringsAsFactors = FALSE
    )
}))
metab_long <- unique(metab_long)  # multiple compartments collapse to the same (timepoint, chebi, score)

metab_wide <- reshape(metab_long, idvar = "chebi_id", timevar = "timepoint_h", direction = "wide")
names(metab_wide) <- sub("^score\\.", "t", names(metab_wide))
metab_wide$name <- metab_name[metab_wide$chebi_id]
score_cols <- grep("^t[0-9]", names(metab_wide), value = TRUE)
metab_wide$n_timepoints <- rowSums(!is.na(metab_wide[score_cols]))
metab_wide$range <- apply(metab_wide[score_cols], 1, function(x) diff(range(x, na.rm = TRUE)))
metab_wide$max_abs <- apply(metab_wide[score_cols], 1, function(x) max(abs(x), na.rm = TRUE))
metab_wide <- metab_wide[order(-metab_wide$n_timepoints, -metab_wide$range), c("chebi_id", "name", score_cols, "n_timepoints", "range", "max_abs")]

cat(sprintf(
    "\n=== %d metabolites pass pruning in >=1 timepoint (of %d measured); shown by # timepoints present, then swing ===\n",
    nrow(metab_wide), length(unique(metab_long$chebi_id))
))
print(metab_wide, row.names = FALSE)
write.csv(metab_wide, "result/moon/metabolite_trajectories.csv", row.names = FALSE)

## ---------------------------------------------------------------------
## 3. Key drivers: TF/kinase upstream_input nodes, trajectory across
## timepoints (footprint activity score, not raw measurement).
## ---------------------------------------------------------------------

driver_long <- do.call(rbind, lapply(names(moon_results), function(tp) {
    r <- moon_results[[tp]]
    drv <- r$nodes[r$nodes$type == "upstream_input", ]
    if (nrow(drv) == 0) return(NULL)
    data.frame(timepoint_h = as.numeric(tp), uniprot = drv$source, score = drv$score, stringsAsFactors = FALSE)
}))

driver_wide <- reshape(driver_long, idvar = "uniprot", timevar = "timepoint_h", direction = "wide")
names(driver_wide) <- sub("^score\\.", "t", names(driver_wide))
driver_wide$symbol <- gene_symbol[driver_wide$uniprot]
score_cols <- grep("^t[0-9]", names(driver_wide), value = TRUE)
driver_wide$range <- apply(driver_wide[score_cols], 1, function(x) diff(range(x, na.rm = TRUE)))
driver_wide <- driver_wide[order(-driver_wide$range), c("uniprot", "symbol", score_cols, "range")]

cat("\n=== Top 15 TF/kinase drivers by largest swing in activity score across timepoints ===\n")
print(head(driver_wide, 15), row.names = FALSE)
write.csv(driver_wide, "result/moon/driver_trajectories.csv", row.names = FALSE)

cat("\nSaved result/moon/{network_size_per_timepoint,metabolite_trajectories,driver_trajectories}.csv\n")
