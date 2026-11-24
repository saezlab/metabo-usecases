# cancer-cell-lines

Figure of the cancer cell lines use case. **Not built yet.**

The figure will be rebuilt from the results of `analyses/cancer-cell-lines/`
(see that folder's README and `docs/`). Until then `manuscript.yaml` lists it
as `null`, so `rebuild.R` skips it.

The previous version (KRAS/EGFR volcanos, azelate evidence, COSMOS PKN
panels) was removed in the 2026-11 cleanup; see `RESTRUCTURING.md`. The
figure of the first preprint is at the tag `preprint_v0`.

When adding the figure, follow the other folders in `figures/`: `build.R`,
`caption.tex`, this README, `manual/` for hand-made assets, outputs in `out/`
(gitignored). Then set its number in `manuscript.yaml`.
