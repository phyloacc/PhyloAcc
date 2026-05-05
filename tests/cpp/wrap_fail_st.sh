#!/usr/bin/env bash
set -e
OUTPUT=$(EXPECT_FAIL=1 ./tests/cpp/phyloacc_cpp_fail_st 2>&1 || true)
printf '%s\n' "$OUTPUT"
if [[ "$OUTPUT" != *"Cannot open the phylogenetic tree input file"* ]]; then
  echo "Missing expected ST error output" >&2
  exit 1
fi
echo "Observed expected ST failure."
