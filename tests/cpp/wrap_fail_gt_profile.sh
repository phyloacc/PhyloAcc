#!/usr/bin/env bash
set -e
TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT
ALN="$TMPDIR/malformed.fa"
BED="$TMPDIR/segments.bed"
printf '>sp1\nAAAA\n>sp2\nAAA\n>sp3\nAAAA\n' > "$ALN"
printf 'chr1\t0\t4\tlocus1\n' > "$BED"
OUTPUT=$(EXPECT_FAIL=1 MALFORMED_ALN="$ALN" MALFORMED_BED="$BED" ./tests/cpp/phyloacc_cpp_fail_gt_profile 2>&1 || true)
printf '%s\n' "$OUTPUT"
if [[ "$OUTPUT" != *"wholeline.length() == prof.G"* ]]; then
  echo "Missing expected GT profile assertion output" >&2
  exit 1
fi
echo "Observed expected GT profile failure."
