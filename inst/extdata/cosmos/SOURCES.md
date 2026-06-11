# Vendored assets: COSMOS+ CSVs

## cosmos_plus_human.csv

- Organism: Homo sapiens (NCBI taxonomy 9606)
- CLI command: `cosmos-pkn export --all-columns --organism 9606 --output cosmos_plus_human.csv`
- Package: omnipath-metabo (version: TBD — fill after generation on beauty)
- Git commit: TBD
- Generation date: TBD
- File fingerprint (SHA-256): TBD

## cosmos_plus_mouse.csv

- Organism: Mus musculus (NCBI taxonomy 10090)
- CLI command: `cosmos-pkn export --all-columns --organism 10090 --output cosmos_plus_mouse.csv`
- Package: omnipath-metabo (version: TBD — fill after generation on beauty)
- Git commit: TBD
- Generation date: TBD
- File fingerprint (SHA-256): TBD

## Columns (--all-columns)

`source`, `target`, `sign`, `interaction_type`, `resource`, `source_type`, `target_type`, `locations`

## Refresh policy

Re-vendor whenever `omnipath-metabo` releases a new COSMOS PKN build.
After re-generation, update the version, commit hash, date, and fingerprints above,
then re-run `./rebuild.sh fig04-cosmos-pkn`.

## Generation instructions

On beauty (or any host with omnipath-metabo installed and the DB reachable):

```bash
cd ~/omnipath-metabo/usecases
cosmos-pkn export --all-columns --organism 9606  --output inst/extdata/cosmos/cosmos_plus_human.csv
cosmos-pkn export --all-columns --organism 10090 --output inst/extdata/cosmos/cosmos_plus_mouse.csv
sha256sum inst/extdata/cosmos/cosmos_plus_human.csv inst/extdata/cosmos/cosmos_plus_mouse.csv
# Paste fingerprints above and commit.
```
