#!/usr/bin/env bash
set -euo pipefail
root="/n/holylfs05/LABS/informatics/Lab/projects/gwct/phyloacc/PhyloAcc/dev/phyloacc-buildcheck4"
export LIBRARY_PATH="$root/x86_64-conda-linux-gnu/sysroot/usr/lib64:$root/lib/gcc/x86_64-conda-linux-gnu/14.3.0:${LIBRARY_PATH:-}"
exec "$root/bin/x86_64-conda-linux-gnu-c++" \
  -B "$root/bin" \
  -B "$root/lib/gcc/x86_64-conda-linux-gnu/14.3.0" \
  -B "$root/libexec/gcc/x86_64-conda-linux-gnu/14.3.0" \
  -isystem "$root/lib/gcc/x86_64-conda-linux-gnu/14.3.0/include" \
  -isystem "$root/lib/gcc/x86_64-conda-linux-gnu/14.3.0/include/c++" \
  -isystem "$root/lib/gcc/x86_64-conda-linux-gnu/14.3.0/include/c++/x86_64-conda-linux-gnu" \
  -isystem "$root/lib/gcc/x86_64-conda-linux-gnu/14.3.0/include/c++/backward" \
  -isystem "$root/x86_64-conda-linux-gnu/sysroot/usr/include" \
  "$@"
