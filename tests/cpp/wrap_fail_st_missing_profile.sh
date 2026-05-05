#!/usr/bin/env bash
set -e
OUTPUT=$(EXPECT_FAIL=1 MALFORMED_ALN=tests/data/minimal/does_not_exist.fa MALFORMED_BED=tests/data/minimal/bed.bed ./tests/cpp/phyloacc_cpp_fail_st_profile 2>&1 || true)
printf '%s
' "$OUTPUT"
if [[ "$OUTPUT" != *"Cannot open the phylogenetic profile input file"* ]]; then
  echo "Missing expected ST missing-profile output" >&2
  exit 1
fi
echo "Observed expected ST missing-profile failure."
