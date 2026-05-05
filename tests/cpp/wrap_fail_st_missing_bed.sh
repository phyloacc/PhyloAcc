#!/usr/bin/env bash
set -e
OUTPUT=$(EXPECT_FAIL=1 MALFORMED_ALN=tests/data/minimal/aln.fa MALFORMED_BED=tests/data/minimal/does_not_exist.bed ./tests/cpp/phyloacc_cpp_fail_st_profile 2>&1 || true)
printf '%s
' "$OUTPUT"
if [[ "$OUTPUT" != *"Cannot open the segment input file"* ]]; then
  echo "Missing expected ST missing-bed output" >&2
  exit 1
fi
echo "Observed expected ST missing-bed failure."
