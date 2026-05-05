#!/usr/bin/env bash
set -e
OUTPUT=$(EXPECT_FAIL=1 MALFORMED_ALN=tests/data/minimal/does_not_exist.fa MALFORMED_BED=tests/data/minimal/bed.bed ./tests/cpp/phyloacc_cpp_fail_gt_profile 2>&1 || true)
printf '%s
' "$OUTPUT"
if [[ "$OUTPUT" != *"Cannot open the phylogenetic profile input file"* ]]; then
  echo "Missing expected GT missing-profile output" >&2
  exit 1
fi
echo "Observed expected GT missing-profile failure."
