#!/bin/bash
# Record benchmark results tagged with the current git commit.
# Usage: tools/record_bench.sh [note]
set -e
cd "$(dirname "$0")/.."

HASH=$(git rev-parse --short HEAD)
DIRTY=$(git diff --quiet && git diff --cached --quiet || echo "+dirty")
NOTE="${1:-}"
TS=$(date -u +%Y-%m-%dT%H:%M:%SZ)

OUT=$(cargo bench --features lvl1,sign,bench-internals 2>/dev/null | grep -E '^(verify|keygen|sign)\b')
# Parse "min=   X.XXXms" — minimum is most stable under system load.
VERIFY=$(echo "$OUT" | sed -n 's/^verify.*min= *\([^ ]*\).*/\1/p')
KEYGEN=$(echo "$OUT" | sed -n 's/^keygen.*min= *\([^ ]*\).*/\1/p')
SIGN=$(echo "$OUT" | sed -n 's/^sign.*min= *\([^ ]*\).*/\1/p')

LINE="$HASH$DIRTY,$TS,verify=$VERIFY,keygen=$KEYGEN,sign=$SIGN,$NOTE"
echo "$LINE" | tee -a BENCH_HISTORY.csv
