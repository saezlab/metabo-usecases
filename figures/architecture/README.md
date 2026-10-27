# architecture — Architecture diagram + statistics digest

**Mode**: manual (architecture asset vendored from Inkscape) + pipeline
side-output (FR-043 statistics digest).

Post-2026-06-14 six-figure renumbering: Figure 1 is the standalone
workflow-architecture diagram. The previous quantitative panels that used
to live next to it are now Figure 2 (under
`figures/database-content/`).

## Inputs

- `figures/architecture/manual/omnipath-architecture-new2026.{pdf,svg}`
  — vendored Inkscape source, SHA-256-pinned in the folder's README.md
  (FR-005, FR-005a).
- `digest-config.yaml` under `panel-a-stats/` — definition variants for the
  FR-043 digest (transporter/receptor curated resource lists, pathway/
  reaction preferred-vs-fallback switches, slot-3 fallback chain).

## Outputs

- `out/architecture.pdf` — verbatim copy of the vendored Inkscape PDF.
- `out/architecture.svg` — verbatim copy of the vendored Inkscape SVG.
- `out/architecture-with-caption.pdf` — typeset caption-and-figure PDF
  (FR-041).
- `out/caption.txt` — plain-text caption (FR-041a).
- `out/architecture.pdf.provenance.json` — sidecar carrying the dev5
  build_id, the architecture asset fingerprint, and a pointer at the digest
  sidecar.
- `panel-a-stats/{stats.json, stats.csv, stats.md, stats.pdf,
  stats.provenance.json}` — the FR-043 statistics digest (FR-043 family).
  The author hand-transcribes these numbers into the Inkscape diagram
  before re-vendoring.

## How to rebuild

```sh
./rebuild.sh architecture
```

This verifies the architecture asset fingerprint, re-runs the digest
against dev5, recompiles the caption, and writes the sidecar.

## Updating the architecture diagram

When the Inkscape source is re-edited and a new PDF/SVG exported:

1. Copy the new files to
   `figures/architecture/manual/omnipath-architecture-new2026.{pdf,svg}`.
2. Recompute the SHA-256 fingerprint and update the value in
   `figures/architecture/manual/README.md` (also bump the
   `last_edit_date` field).
3. Re-run `./rebuild.sh architecture`. The rebuild fails fast if the
   PDF on disk does not match the README fingerprint.
