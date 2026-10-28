# Reading cell line annotation from the local Cellosaurus flat file.
#
# The flat file (https://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt,
# header removed) is used instead of the Cellosaurus API so that results are
# tied to one release and do not change between runs.

#' Read the parts of the Cellosaurus flat file needed here
#'
#' @param path Path to the prepared flat file.
#' @return A list with
#'   `entries`: one row per entry (ac, id, synonyms, human, disease,
#'   derived_from_site, category, parent);
#'   `variations`: one row per "Sequence variation" comment (ac, text).
read_cellosaurus <- function(path) {
    lines <- readLines(path, encoding = "UTF-8", warn = FALSE)
    # an entry ends with "//"; lines of one entry share the same index
    entry <- cumsum(c(0L, head(lines == "//", -1)))
    tag <- substr(lines, 1, 2)
    value <- substring(lines, 6)

    keep <- tag %in% c("ID", "AC", "SY", "OX", "DI", "CA", "HI") |
        startsWith(lines, "CC   Sequence variation:") |
        startsWith(lines, "CC   Derived from site:")
    d <- data.frame(entry = entry[keep], tag = tag[keep], value = value[keep])

    first_of <- function(t) {
        x <- d[d$tag == t, ]
        x$value[match(sort(unique(d$entry)), x$entry)]
    }
    all_of <- function(t, sep = " | ", pattern = NULL) {
        x <- d[d$tag == t, ]
        if (!is.null(pattern)) x <- x[startsWith(x$value, pattern), ]
        v <- tapply(x$value, x$entry, paste, collapse = sep)
        unname(v[as.character(sort(unique(d$entry)))])
    }

    entries <- data.frame(
        entry = sort(unique(d$entry)),
        id = first_of("ID"),
        # the first AC is the primary accession
        ac = sub(";.*$", "", first_of("AC")),
        synonyms = all_of("SY", sep = "; "),
        human = grepl("NCBI_TaxID=9606;", all_of("OX")),
        disease = all_of("DI"),
        derived_from_site = all_of("CC", pattern = "Derived from site:"),
        category = first_of("CA"),
        parent = all_of("HI")
    )
    entries$derived_from_site <- sub("^Derived from site: ", "", entries$derived_from_site)

    variations <- d[d$tag == "CC" & startsWith(d$value, "Sequence variation:"), ]
    variations <- data.frame(
        ac = entries$ac[match(variations$entry, entries$entry)],
        text = sub("^Sequence variation: ", "", variations$value)
    )

    list(entries = entries[, -1], variations = variations)
}

#' Split Cellosaurus sequence variations into fields and classify them
#'
#' Formats (Cellosaurus release 53):
#'   Mutation; HGNC; HGNC:id; GENE; Kind; Variant; [ClinVar=..;] Zygosity=..[; Note=..] (Source).
#'   Gene deletion; HGNC; HGNC:id; GENE; Zygosity=.. (Source).
#'   Gene amplification; HGNC; HGNC:id; GENE; Kind; Zygosity=.. (Source).
#'   Gene fusion; HGNC; HGNC:id; GENE + HGNC; HGNC:id; PARTNER; Name(s)=.. (Source).
#' The kind "None_reported" records that no mutation was found in the gene,
#' i.e. reported wild type; "_edited" kinds are engineered mutations.
#'
#' @param variations Data frame with columns `ac` and `text`.
#' @return The input with columns variation_type, gene, partner_gene, kind,
#'   variant, zygosity, note, source and status.
parse_sequence_variations <- function(variations) {
    text <- sub("\\.$", "", variations$text)
    # trailing "(PubMed=..; ..)" holds the sources and may contain "; "
    source <- ifelse(grepl("\\(([^()]*)\\)$", text), sub("^.*\\(([^()]*)\\)$", "\\1", text), NA)
    body <- sub("\\s*\\([^()]*\\)$", "", text)
    parts <- strsplit(body, "; ", fixed = TRUE)
    field <- function(i) vapply(parts, function(p) if (length(p) >= i) p[[i]] else NA_character_, "")
    keyval <- function(key) {
        m <- regmatches(body, regexec(paste0(key, "=([^;]+)"), body))
        vapply(m, function(x) if (length(x) == 2) x[[2]] else NA_character_, "")
    }

    type <- field(1)
    is_fusion <- type == "Gene fusion"
    is_mutation <- type == "Mutation"
    out <- data.frame(
        ac = variations$ac,
        variation_type = type,
        gene = ifelse(is_fusion, sub(" \\+ .*$", "", field(4)), field(4)),
        partner_gene = ifelse(is_fusion, field(6), NA_character_),
        kind = ifelse(is_mutation | type == "Gene amplification", field(5), NA_character_),
        variant = ifelse(is_mutation, field(6), NA_character_),
        zygosity = keyval("Zygosity"),
        note = keyval("Note"),
        source = source
    )
    out$status <- dplyr::case_when(
        is_mutation & startsWith(out$kind, "None_reported") ~ "wild_type_reported",
        is_mutation & grepl("_edited$", out$kind) ~ "mutation_engineered",
        is_mutation & startsWith(out$kind, "Repeat_expansion") ~ "repeat_expansion",
        is_mutation ~ "mutation",
        type == "Gene deletion" ~ "deletion",
        type == "Gene amplification" ~ "amplification",
        is_fusion ~ "fusion",
        TRUE ~ "other"
    )
    out
}
