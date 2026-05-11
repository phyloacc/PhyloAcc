#!/usr/bin/env bash
set -e
TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT
ALN="$TMPDIR/malformed.fa"
BED="$TMPDIR/segments.bed"
printf '>sp1\nAAAA\n>sp2\nAAA\n>sp3\nAAAA\n' > "$ALN"
printf 'chr1\t0\t4\tlocus1\n' > "$BED"
OUTPUT=$(EXPECT_FAIL=1 MALFORMED_ALN="$ALN" MALFORMED_BED="$BED" ./tests/cpp/phyloacc_cpp_fail_st_profile 2>&1 || true)
printf '%s\n' "$OUTPUT"
if [[ "$OUTPUT" != *"Sequence length mismatch in phylogenetic profile input file"* ]] || [[ "$OUTPUT" != *"Species sp2 has 3 sites; expected 4"* ]]; then
  echo "Missing expected ST profile length-mismatch output" >&2
  exit 1
fi
echo "Observed expected ST profile failure."
