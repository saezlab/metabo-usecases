# Compare the refactored pre-processing (shorthouse_v1) with the previous
# pipeline (analyses/cancer-cell-lines/0[1-3]_*.Rmd, outputs in
# data/derived/cancer-cell-lines/*.RData). Every difference is either
# attributed to an intended change or counted as unexplained.
#
# Run from the repository root after both pipelines have been run:
#   Rscript analyses/cancer-cell-lines/preprocessing/verify_against_previous.R

suppressMessages({
    library(dplyr)
    library(tidyr)
    library(stringr)
})
source(here::here("analyses/cancer-cell-lines/R/use_case.R"))
config <- load_config()
new_dir <- uc_dir(config, "derived")
old_dir <- here::here("data/derived/cancer-cell-lines")
out_dir <- uc_dir(config, "output", "verify_against_previous")

load_one <- function(f) { e <- new.env(); load(f, envir = e); as.list(e)[[1]] }
old <- load_one(file.path(old_dir, "Shorthouse_2022_extended.RData"))
new_lines <- readRDS(file.path(new_dir, "01_lines.rds"))
new_cells <- readr::read_tsv(file.path(new_dir, "02_cell_lines.tsv"), show_col_types = FALSE)
new_genes <- readr::read_tsv(file.path(new_dir, "02_gene_status.tsv"), show_col_types = FALSE)
summary <- list()

# 1. Lines -------------------------------------------------------------------------
old_names <- rownames(old$data_matrix)
member <- new_cells %>% separate_rows(member_names, sep = "; ") %>% select(line, name = member_names)
line_diff <- tibble(name = union(old_names, member$name)) %>%
    mutate(in_old = name %in% old_names, in_new_as = member$line[match(name, member$name)]) %>%
    mutate(explanation = case_when(
        in_old & is.na(in_new_as) & name %in% config$shorthouse$exclude_lines ~ "B14 excluded line",
        in_old & !is.na(in_new_as) & name != in_new_as ~ "B3 merged into another line",
        in_old & !is.na(in_new_as) ~ "same",
        TRUE ~ "UNEXPLAINED"
    ))
summary$lines <- count(line_diff, explanation)

# 2. Intensities --------------------------------------------------------------------
# Old: mean over the injections of each name. New: mean over the injections of
# each line, so a merged line equals the injection-weighted mean of its names.
old_mat <- as.matrix(old$data_matrix)
old_mat <- old_mat[, new_lines$ions$feature_id]
n_inj <- readRDS(file.path(new_dir, "01_injections.rds"))$injections %>% count(name)
expected <- sapply(new_lines$lines$line, function(l) {
    names_l <- member$name[member$line == l]
    w <- n_inj$n[match(names_l, n_inj$name)]
    colSums(old_mat[names_l, , drop = FALSE] * w) / sum(w)
})
expected <- t(expected)
rel <- abs(new_lines$intensities - expected) / pmax(abs(expected), 1)
summary$intensities <- tibble(
    lines = nrow(rel), ions = ncol(rel),
    max_relative_difference = max(rel),
    lines_over_1e_9 = sum(apply(rel, 1, max) > 1e-9)
)

# 3. Mutations: genes counted as mutated before vs. statuses now -----------------------
old_mut <- old$samples %>%
    transmute(name, gene = str_split(coalesce(mutated_genes, ""), ";\\s*")) %>%
    unnest(gene) %>%
    filter(gene != "") %>%
    left_join(member, by = "name") %>%
    distinct(line, gene)
mut_cmp <- full_join(
    old_mut %>% mutate(old_mutated = TRUE),
    new_genes %>% select(line, gene, status),
    by = c("line", "gene")
) %>%
    mutate(
        old_mutated = coalesce(old_mutated, FALSE),
        new_mutated = coalesce(status == "mutation", FALSE),
        accession_changed = line %in% c("977", "NCIH441"),
        explanation = case_when(
            old_mutated == new_mutated ~ "same",
            is.na(line) ~ "B14 excluded line",
            line == "NCIH460DNP53" ~ "B3 engineered line, no parental variations",
            accession_changed ~ "B2/B15 corrected accession",
            old_mutated & status == "wild_type_reported" ~ "B1 reported wild type",
            # the old parser kept the raw fusion field, e.g. "ABL1 + HGNC"
            old_mutated & str_detect(gene, " \\+ HGNC$") ~ "B1 fusion field taken as gene name",
            old_mutated & status %in% c("deletion", "amplification", "fusion", "repeat_expansion", "mutation_engineered") ~
                paste("B1 not a mutation:", status),
            TRUE ~ "UNEXPLAINED"
        )
    )
summary$mutations <- count(mut_cmp, explanation)

# 4. Tissue ---------------------------------------------------------------------------
tissue_cmp <- old$samples %>%
    select(name, old_site = general_site) %>%
    inner_join(member, by = "name") %>%
    left_join(new_cells %>% select(line, new_site = general_site, accession_source), by = "line") %>%
    mutate(explanation = case_when(
        coalesce(old_site, "") == coalesce(new_site, "") ~ "same",
        line %in% c("977", "NCIH441") ~ "B2/B15 corrected accession",
        TRUE ~ "UNEXPLAINED"
    ))
summary$tissue <- count(tissue_cmp, explanation)

# 5. Feature identifiers (only once 03 has run) ---------------------------------------
features_file <- file.path(new_dir, "03_features.rds")
if (file.exists(features_file)) {
    new_feat <- readRDS(features_file)
    id_cols <- c("hmdb", "kegg", "chebi", "pubchem", "metabolite_class")
    norm <- function(x) vapply(str_split(coalesce(as.character(x), ""), ";\\s*"),
                               function(v) paste(sort(unique(v[v != ""])), collapse = ";"), "")
    feat_cmp <- inner_join(
        old$features %>% select(ionIdx, all_of(id_cols)),
        new_feat %>% select(ionIdx, all_of(id_cols)),
        by = "ionIdx", suffix = c("_old", "_new")
    )
    feat_diff <- bind_rows(lapply(id_cols, function(col) {
        tibble(ionIdx = feat_cmp$ionIdx, column = col,
               same = norm(feat_cmp[[paste0(col, "_old")]]) == norm(feat_cmp[[paste0(col, "_new")]]))
    })) %>%
        mutate(explanation = case_when(
            same ~ "same",
            ionIdx %in% c(797, 815) ~ "B6 ions sharing a name no longer exchange identifiers",
            TRUE ~ "UNEXPLAINED"
        ))
    summary$features <- count(feat_diff, column, explanation)
    readr::write_csv(filter(feat_diff, !same), file.path(out_dir, "feature_differences.csv"))
}

readr::write_csv(line_diff, file.path(out_dir, "lines.csv"))
readr::write_csv(filter(mut_cmp, explanation != "same"), file.path(out_dir, "mutation_differences.csv"))
readr::write_csv(filter(tissue_cmp, explanation != "same"), file.path(out_dir, "tissue_differences.csv"))
for (n in names(summary)) { cat("\n##", n, "\n"); print(as.data.frame(summary[[n]]), row.names = FALSE) }
