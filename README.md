# Omnipath Metabo

This repository includes all notebooks used to create figures and data tables for the manuscript. \

## Structure
- The `Omnipath Metabo` database overview and each use-case has their own folder with its own environment management with `uv`, `nix` or `renv`.
- The final supplementary tables should be saved as `ExtendedDataTable_Number_ShortName.xlsx` with a sheet overview in the first sheet:

    | Sheet | Content |
    |---------|---------|
    | Sheet 1 | Rawdata. |
    | Sheet 2 | Normalised data. |
    | Sheet 3 | ... |

- The final figures should be saved as `.svg` or anything else compatible with Inkscape or Adobe Illustrator.
- Ensure that any papers/methods you mention/use in your notebooks are cited (e.g. using `.bib` file)

