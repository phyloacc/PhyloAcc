#!/usr/bin/env bash
set -e
OUTPUT=$(EXPECT_FAIL=1 ./tests/cpp/phyloacc_cpp_fail_config_bool 2>&1 || true)
printf '%s\n' "$OUTPUT"
if [[ "$OUTPUT" != *"Invalid boolean value for parameter WL: maybe"* ]]; then
  echo "Missing expected config boolean error output" >&2
  exit 1
fi
echo "Observed expected config boolean failure."
