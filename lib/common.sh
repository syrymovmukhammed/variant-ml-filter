#!/usr/bin/env bash
# Shared helpers, sourced by every step.
set -euo pipefail

_here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=../config.sh
source "${_here}/config.sh"

log()  { printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*" >&2; }
warn() { printf '[%s] WARNING: %s\n' "$(date '+%H:%M:%S')" "$*" >&2; }
die()  { printf '[%s] ERROR: %s\n' "$(date '+%H:%M:%S')" "$*" >&2; exit 1; }

need_cmd() {
  local c
  for c in "$@"; do
    command -v "$c" >/dev/null 2>&1 || die "'$c' not found in PATH (see README → Requirements)"
  done
}

need_file() {
  local f
  for f in "$@"; do
    [[ -s "$f" ]] || die "missing input: $f — run the previous step first"
  done
}


# Iterate over the cohort table: prints "sample coverage error contam paralogs seed split"
cohort() { awk -F'\t' 'NR>1 && $1!~/^#/ {print $1, $2, $3, $4, $5, $6, $7}' "${SAMPLES_TSV}"; }

# Print a TSV as an aligned table on stderr (no dependency on `column`).
pretty_tsv() { awk -F'\t' '{for(i=1;i<=NF;i++) printf "%-14s", $i; print ""}' "$@" >&2; }
