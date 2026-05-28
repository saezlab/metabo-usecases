# Palette sources

Per FR-020a, every palette consumed from `saezlab/graphics` is copied
into and committed to this repo so the pipeline is self-contained.

| Palette | Upstream | Copied on | Notes |
|---|---|---|---|
| `metabo.gpl` | `saezlab/graphics:palettes/metabo.gpl` | 2026-05-28 | Lead categorical palette authored for the OmniPath Metabo manuscript. Hex order: `#006384`, `#9F0162`, `#FEAF16`, `#BBCC33`, `#EA6572`, `#009E73`, `#99DDFF`, `#D03293`, `#BEBEBE` (grey reserved for None/Unknown/NA per spec FR-020). |
| `rwth.gpl` | `saezlab/graphics:palettes/rwth.gpl` | 2026-05-28 | Extended palette for categorical needs beyond the lead palette (FR-021). |

When the upstream `saezlab/graphics` palettes change, re-copy with:

```sh
cp ~/graphics/palettes/{metabo,rwth}.gpl inst/extdata/palettes/
```

and update the "Copied on" date here.
