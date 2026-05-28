# shell.nix — project-local Nix environment for metabo.figures pipeline.
#
# Beauty is NixOS; the pipeline depends on:
#   - libpq for RPostgres
#   - V8 + nodejs for jsonvalidate / juicyjuice / gt
#   - TeX Live (xelatex + pdfpages + fontspec + tikz + standalone)
#   - pdftk for FR-029 bundle concatenation
#   - graphics tooling already declared elsewhere (uv, R)
#
# Usage on beauty:
#
#   cd ~/dev/metabo-usecases
#   nix-shell
#   # Inside the shell, the env has libpq, V8, xelatex, pdftk on PATH.
#   Rscript -e 'devtools::install_deps(dependencies = TRUE)'
#   ./rebuild.sh --dry-run
#
# Laptops without nix can still work — only the system deps fail to
# install in that case; the R package logic itself remains portable.

{ pkgs ? import <nixpkgs> {} }:

let
  texEnv = pkgs.texlive.combine {
    inherit (pkgs.texlive)
      scheme-medium
      pdfpages
      eso-pic        # required by pdfpages for page overlays
      fontspec
      standalone
      pgf            # tikz
      xcolor
      geometry
      l3packages
      ;
  };
in

pkgs.mkShell {
  name = "metabo-figures";

  buildInputs = with pkgs; [
    # R is installed system-wide on beauty already, but pin here for
    # laptops that lack it.
    R

    # Python helpers via uv (kept out of nix — uv manages its own venv).
    uv
    python311

    # System libs for R packages that won't build without them.
    postgresql.dev   # libpq for RPostgres
    nodejs           # provides V8 for jsonvalidate / juicyjuice / gt
    icu              # text utilities used by stringi / R
    libxml2
    openssl
    curl
    pkg-config

    # Composition + manuscript-bundle tooling.
    texEnv
    pdftk
    poppler-utils    # pdftoppm for --png PNG companions

    # Quality-of-life tools used by quickstart.md / docs.
    git
    rsync
    jq
    yq-go
  ];

  shellHook = ''
    echo "[shell.nix] metabo.figures dev environment ready"
    echo "  R:        $(R --version 2>&1 | head -1)"
    echo "  xelatex:  $(xelatex --version 2>&1 | head -1)"
    echo "  pdftk:    $(pdftk --version 2>&1 | head -1)"
    echo "  Postgres include path for RPostgres: ${pkgs.postgresql.dev}/include"
    export PKG_CONFIG_PATH="${pkgs.postgresql.dev}/lib/pkgconfig:$PKG_CONFIG_PATH"
  '';
}
