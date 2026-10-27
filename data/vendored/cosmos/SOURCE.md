# Vendored asset: meta_network.RData

Source:  R/Bioconductor package, 
cosmosR version: 1.18.1
Retrieval date: 2026-06-10
File fingerprint (MD5): 

## What this file is

The old (pre-OmniPath Metabo) COSMOS prior-knowledge network serialised as
an R data object. Used by Figure 4 (FR-011, FR-011a) as the baseline to
compare against the new COSMOS PKN built by the  package.

The old PKN contains only edges derived from metabolic reactions. The new
PKN additionally includes receptors, inhibitors, allosteric regulators,
transport, and other metabolite-protein interaction types — the comparison
panel makes this scope extension explicit.

## Refresh policy

Re-vendor when cosmosR releases a new version that updates this object, or
when the upstream GitHub file changes
(https://github.com/saezlab/cosmosR/blob/master/data/meta_network.RData).
After refresh, recompute MD5 above and re-run ./rebuild.sh cosmos-pkn.
