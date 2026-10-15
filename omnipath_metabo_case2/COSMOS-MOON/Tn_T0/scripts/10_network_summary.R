# User Story 3 capstone (spec 003-case-study-2-temporal-trajectory FR-013):
# summary statistics and key-driver/key-metabolite tables across the
# 14-(genotype,timepoint)-pair pruned mechanistic network series.
#
# Mirrors 002's 10_network_summary.R, with two structural changes (research.md
# R4): moon_results is keyed by "WT_2".."ob_24" (14 keys), not bare
# timepoints (8 keys) -- parsed into explicit genotype + timepoint_h columns
# everywhere `as.numeric(tp)` used to work directly; and the metabolite/
# driver trajectory tables keep an explicit `genotype` column, pivoting wide
# on timepoint_h WITHIN each genotype (not a combined genotype_timepoint
# axis) -- so one row is "this metabolite's trajectory within this genotype."
#
# Run from omnipath_metabo_case2/, after scripts/06_footprint_moon.R.

moon_results <- readRDS("Tn_T0/result/moon/all_pairs.rds")
measured_features <- readRDS("Tn_T0/result/pk_retrieval/measured_features.rds")

parse_key <- function(keys) {
    data.frame(
        genotype = sub("_.*", "", keys),
        timepoint_h = as.numeric(sub(".*_", "", keys)),
        stringsAsFactors = FALSE
    )
}

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
## 1. Network size per (genotype, timepoint) pair
## ---------------------------------------------------------------------

size_table <- do.call(rbind, lapply(names(moon_results), function(key) {
    r <- moon_results[[key]]
    gt <- parse_key(key)
    data.frame(
        genotype = gt$genotype,
        timepoint_h = gt$timepoint_h,
        pruned_nodes = nrow(r$nodes),
        direct_edges = nrow(r$edges),
        gem_edges_reattached = nrow(r$gem_edges),
        metabolite_nodes = sum(grepl("^Metab__", r$nodes$source)),
        tf_kinase_upstream_nodes = sum(r$nodes$type == "upstream_input")
    )
}))
size_table <- size_table[order(size_table$genotype, size_table$timepoint_h), ]

cat("=== Network size per (genotype, timepoint) pair ===\n")
print(size_table, row.names = FALSE)
write.csv(size_table, "Tn_T0/result/moon/network_size_per_pair.csv", row.names = FALSE)

## ---------------------------------------------------------------------
## 2. Key metabolites: trajectory of MOON score across timepoints, within
## each genotype (level0 nodes = directly measured metabolite t-stats).
## ---------------------------------------------------------------------

metab_long <- do.call(rbind, lapply(names(moon_results), function(key) {
    r <- moon_results[[key]]
    gt <- parse_key(key)
    met <- r$nodes[grepl("^Metab__", r$nodes$source) & r$nodes$type == "level0", ]
    if (nrow(met) == 0) return(NULL)
    data.frame(
        genotype = gt$genotype,
        timepoint_h = gt$timepoint_h,
        chebi_id = strip_compartment(met$source),
        score = met$score,
        stringsAsFactors = FALSE
    )
}))
metab_long <- unique(metab_long)  # multiple compartments collapse to the same (genotype, timepoint, chebi, score)

metab_wide <- reshape(metab_long, idvar = c("genotype", "chebi_id"), timevar = "timepoint_h", direction = "wide")
names(metab_wide) <- sub("^score\\.", "t", names(metab_wide))
metab_wide$name <- metab_name[metab_wide$chebi_id]
score_cols <- grep("^t[0-9]", names(metab_wide), value = TRUE)
metab_wide$n_timepoints <- rowSums(!is.na(metab_wide[score_cols]))
metab_wide$range <- apply(metab_wide[score_cols], 1, function(x) diff(range(x, na.rm = TRUE)))
metab_wide$max_abs <- apply(metab_wide[score_cols], 1, function(x) max(abs(x), na.rm = TRUE))
metab_wide <- metab_wide[order(metab_wide$genotype, -metab_wide$n_timepoints, -metab_wide$range),
                          c("genotype", "chebi_id", "name", score_cols, "n_timepoints", "range", "max_abs")]

cat(sprintf(
    "\n=== %d (genotype, metabolite) rows pass pruning in >=1 timepoint; shown by genotype, then # timepoints present, then swing ===\n",
    nrow(metab_wide)
))
print(metab_wide, row.names = FALSE)
write.csv(metab_wide, "Tn_T0/result/moon/metabolite_trajectories.csv", row.names = FALSE)

## ---------------------------------------------------------------------
## 3. Key drivers: TF/kinase upstream_input nodes, trajectory across
## timepoints within each genotype (footprint activity score).
## ---------------------------------------------------------------------

driver_long <- do.call(rbind, lapply(names(moon_results), function(key) {
    r <- moon_results[[key]]
    gt <- parse_key(key)
    drv <- r$nodes[r$nodes$type == "upstream_input", ]
    if (nrow(drv) == 0) return(NULL)
    data.frame(genotype = gt$genotype, timepoint_h = gt$timepoint_h, uniprot = drv$source, score = drv$score, stringsAsFactors = FALSE)
}))

driver_wide <- reshape(driver_long, idvar = c("genotype", "uniprot"), timevar = "timepoint_h", direction = "wide")
names(driver_wide) <- sub("^score\\.", "t", names(driver_wide))
driver_wide$symbol <- gene_symbol[driver_wide$uniprot]
score_cols <- grep("^t[0-9]", names(driver_wide), value = TRUE)
driver_wide$range <- apply(driver_wide[score_cols], 1, function(x) diff(range(x, na.rm = TRUE)))
driver_wide <- driver_wide[order(driver_wide$genotype, -driver_wide$range), c("genotype", "uniprot", "symbol", score_cols, "range")]

cat("\n=== Top 15 TF/kinase drivers by largest swing in activity score, per genotype ===\n")
for (gt in sort(unique(driver_wide$genotype))) {
    cat(sprintf("\n-- %s --\n", gt))
    print(head(driver_wide[driver_wide$genotype == gt, ], 15), row.names = FALSE)
}
write.csv(driver_wide, "Tn_T0/result/moon/driver_trajectories.csv", row.names = FALSE)

cat("\nSaved Tn_T0/result/moon/{network_size_per_pair,metabolite_trajectories,driver_trajectories}.csv\n")
