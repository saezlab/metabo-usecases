# architecture

The OmniPath Metabo architecture diagram (hand-made in Inkscape), plus a
digest of headline database statistics whose numbers are copied into the
diagram by hand.

**Status (2026-11-24): fails.** The diagram is copied, but the statistics
digest stops because dev5 lacks `metabo_entity_structural_specificity`
(Structures section); see `docs/DB_ACCESS.md` for what the backends
provide.

## Inputs

| Input | What |
|---|---|
| `manual/omnipath-architecture-new2026.{pdf,svg}` | The diagram; SHA-256 pinned in `manual/README.md` and checked by `build.R` |
| `digest-config.yaml` | Definitions of the digest's counts (e.g. which resources count as transporter or receptor sources) |
| dev5 Postgres | Counts for the digest (`R/stats-*.R`) |

## Outputs (`out/`, gitignored)

| File | What |
|---|---|
| `architecture.{pdf,svg}` | Copy of the diagram |
| `architecture-with-caption.pdf`, `caption.txt` | Figure with typeset caption; plain-text caption |
| `panel-a-stats/stats.{json,csv,md,pdf}` | The statistics digest (schema `inst/extdata/contracts/panel-a-stats.schema.json`) |
| `*.provenance.json` | Provenance sidecars |

## Run

```sh
./rebuild.sh architecture
```

## Updating the diagram

1. Edit the Inkscape source; transfer the numbers from
   `out/panel-a-stats/stats.md`.
2. Export to `manual/omnipath-architecture-new2026.{pdf,svg}`.
3. Update the SHA-256 in `manual/README.md`; the build fails if it does
   not match.
4. `./rebuild.sh architecture`.
