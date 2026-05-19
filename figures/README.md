# OmniPath Database Figures

This folder contains a small data-driven D3 dashboard for summarizing the live Postgres database configured in `.env`.

## Rebuild the figure data

```sh
npm run data
```

The export script reads `DATABASE_URL` from `.env`, runs [sql/figure_data.sql](sql/figure_data.sql), and writes [data/figures.json](data/figures.json).

## View the dashboard

```sh
npm run dev
```

Then open [http://127.0.0.1:8080](http://127.0.0.1:8080).

The visualization is static after export: the browser only reads `data/figures.json`. When the database changes, rerun `npm run data` and refresh the page.
