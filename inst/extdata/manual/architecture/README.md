# OmniPath architecture diagram — vendored manual asset (Figure 1 Panel A)

Per FR-005 / FR-005a / SC-010 (2026-06-09 follow-up), Figure 1 Panel A
is a **vendored manual asset** authored in Inkscape outside this
repository. The pipeline embeds the PDF verbatim at the composite's
Panel A position; it never regenerates the artwork.

A swap, silent re-edit, or content drift since the last pin causes the
rebuild to fail (FR-005a). To refresh the asset:

1. Re-export the PDF + SVG from Inkscape.
2. Copy both files over the vendored copies in this folder.
3. Recompute the PDF SHA-256 (`sha256sum *.pdf`) and update the
   `pdf_sha256` line below.
4. Bump the `last_edit_date` field.
5. Re-run `./rebuild.sh fig01-overview` — the rebuild fails fast if
   the PDF on disk does not match the pinned `pdf_sha256`.

## Files

- `omnipath-architecture-new2026.pdf` — composite-embedded PDF
- `omnipath-architecture-new2026.svg` — Inkscape source (vendored for
  reference / re-edit; not consumed by the composite)

## Provenance

- `source_path`: `/home/denes/omnipath_figures/omnipath-architecture-new2026.{pdf,svg}`
- `contributor`: denes
- `last_edit_date`: 2026-06-08
- `reference_style`: `/home/denes/omnipath_figures/omnipath-architecture-legacy2025-v1.svg`
- `pdf_sha256`: `4fd43e054435c4e16353c2c32a746becf323bc59d4bf60d59fd8c1f183615638`
- `svg_sha256`: `bf9807e60b2888fb2abda7b98c086d53a8544569acd578d9818b7ab2bd900312`

## Why a manual asset and not pipeline-generated?

The 2026-06-09 follow-up retired the prior FR-005..FR-006a TikZ
authoring track. The architecture diagram is now hand-arranged in
Inkscape — it gives the contributor full control over visual identity
(box typography, logo placement, downstream-triangle geometry) without
the constraints of programmatic layout. The pipeline's role is to
guarantee the embedded PDF is byte-equivalent to the version the
contributor signed off on, via the SHA-256 fingerprint above.
