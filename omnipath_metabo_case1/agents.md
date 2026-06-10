# Agent Notes

## Beauty DB access for the figures pipeline

When running the `metabo-usecases` figures pipeline against OmniPath deployments, prefer running on `beauty` because the Postgres instances are exposed only on loopback there.

Documented source of truth:
- `DB_ACCESS.md`
- `CONFIGURATION.md`

Working beauty setup:
- host: `localhost`
- database: `omnipath`
- user: `omnipath`
- password: `omnipath`
- ports:
  - `dev3`: `5403`
  - `dev4`: `5404`
  - `dev5`: `5405`

Reusable beauty-side wrapper:
```bash
~/.local/bin/with-metabo-figures-env <command>
```

That wrapper sources:
```bash
~/.config/metabo-figures/env.sh
```

So future noninteractive runs should use, for example:
```bash
cd /home/omnipath/dev/metabo-usecases-fig3
~/.local/bin/with-metabo-figures-env ./rebuild.sh fig03-metalinks-versions
```

If you need direct R execution during development:
```bash
cd /home/omnipath/dev/metabo-usecases-fig3
~/.local/bin/with-metabo-figures-env Rscript -e 'pkgload::load_all(".", export_all = FALSE, quiet = TRUE); source("figures/fig03-metalinks-versions/build.R")'
```

Do not use public hostnames like `dev3.omnipathdb.org:5403` directly from outside beauty. Use beauty local execution or an SSH tunnel as described in `DB_ACCESS.md`.
