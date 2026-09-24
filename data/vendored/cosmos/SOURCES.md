# Vendored assets: COSMOS+ CSVs

## cosmos_plus_human.csv

- Organism: Homo sapiens (NCBI taxonomy 9606)
- CLI command: `cosmos-pkn export --all-columns --organism 9606 --no-stitch --output data/vendored/cosmos/cosmos_plus_human.csv`
- Package: omnipath-metabo git commit `bc3ae3f7ecdf27279176b5f1b6bbc390e2260310`
- Generation date: 2026-06-11
- File fingerprint (SHA-256): `49fd4948cd7ecf12f94a4e270b9ae8841b9db5e53caf4538ec1546f9821f8d31`

## cosmos_plus_mouse.csv

- Organism: Mus musculus (NCBI taxonomy 10090)
- CLI command: `cosmos-pkn export --all-columns --organism 10090 --no-stitch --output data/vendored/cosmos/cosmos_plus_mouse.csv`
- Package: omnipath-metabo git commit `bc3ae3f7ecdf27279176b5f1b6bbc390e2260310`
- Generation date: 2026-06-11
- File fingerprint (SHA-256): `aa94c9336ec901ed97bd5116ff93cb40dabe9a6834c56e932a6213c435a71f3c`

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
cosmos-pkn export --all-columns --organism 9606  --output data/vendored/cosmos/cosmos_plus_human.csv
cosmos-pkn export --all-columns --organism 10090 --output data/vendored/cosmos/cosmos_plus_mouse.csv
sha256sum data/vendored/cosmos/cosmos_plus_human.csv data/vendored/cosmos/cosmos_plus_mouse.csv
# Paste fingerprints above and commit.
```
