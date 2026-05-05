#!/usr/bin/env bash
set -e
OUTPUT=$(EXPECT_FAIL=1 ./tests/cpp/phyloacc_cpp_fail_gt 2>&1 || true)
printf '%s\n' "$OUTPUT"
if [[ "$OUTPUT" != *"Cannot open tree2 input file"* ]]; then
  echo "Missing expected GT error output" >&2
  exit 1
fi
echo "Observed expected GT failure."
