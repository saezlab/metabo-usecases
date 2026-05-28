#!/usr/bin/env bash
set -euo pipefail
exec Rscript "$(dirname "$0")/rebuild.R" "$@"
