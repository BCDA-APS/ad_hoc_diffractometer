# gh-pages root assets

Static files in this directory are published to the **root** of the
`gh-pages` branch (not under `latest/` or any version directory) by the
"Deploy gh-pages root redirect" step in `.github/workflows/docs.yml` on
every push to `main`.

Currently this directory contains:

- `index.html` &mdash; meta-refresh redirect from
  <https://bcda-aps.github.io/ad_hoc_diffractometer/> to
  <https://bcda-aps.github.io/ad_hoc_diffractometer/latest/>, so any
  reference to the bare repo URL lands on the latest published docs.

Add additional files here only if they should appear at the gh-pages
site root (e.g. a `robots.txt` or a `CNAME`).  Per-version assets
belong in `docs/source/_static/` instead.
