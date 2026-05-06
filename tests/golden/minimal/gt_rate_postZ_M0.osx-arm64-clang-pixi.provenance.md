GT Golden Provenance: macOS arm64 Pixi/Clang

Golden file
- `tests/golden/minimal/gt_rate_postZ_M0.osx-arm64-clang-pixi.txt`

Status
- Recorded on `2026-05-06`.
- Golden SHA256: `5be539b7e13db75218ca67cc801e82a336d5db38613252a828ab5fc751aed848`

What matches this golden
- Local macOS arm64 Pixi environment in this repository.
- `PhyloAcc-GT` built by `pixi run build`.
- Compiler: `clang version 19.1.7`, target `arm64-apple-darwin20.0.0`.

Build provenance
- Platform: `darwin`, `arm64`.
- Pixi lock SHA256: `8a2d224c7e4dc8601b9cc52f90a7467356b92ff111f40e16311562ef2c9de623`
- `PhyloAcc-GT` binary SHA256: `e6e21c92fb66d045a6a76bea13878757a05c3519afc24e45a8a294961b2a2636`
- Compiler executable: `arm64-apple-darwin20.0.0-clang++`
- Compiler version output:
  - `clang version 19.1.7`
  - `Target: arm64-apple-darwin20.0.0`

Repeatability checks
- Recording run: `PHYLOACC_RUN_GT=1 PHYLOACC_RECORD_GOLDEN=1 pytest -q tests/test_st_gt.py::test_gt_minimal_run`
  - Result: `1 passed in 104.53s`
- First comparison run: `pixi run test-gt`
  - Result: `1 passed in 101.90s`
- Second comparison run: `pixi run test-gt`
  - Result: `1 passed in 109.96s`

Selection behavior
- `tests/test_st_gt.py` selects this golden automatically when all of the following are true:
  - platform is Darwin
  - machine is arm64
  - `PIXI_PROJECT_NAME=phyloacc`
  - `CXX` contains `clang++`
- To force this baseline, set `PHYLOACC_GT_GOLDEN_ID=osx-arm64-clang-pixi`.
- To force the packaged GCC14 baseline, set `PHYLOACC_GT_GOLDEN_ID=packaged-gcc14`.
- To compare against any explicit file, set `PHYLOACC_GT_GOLDEN=/path/to/golden.txt`.

Relationship to packaged-build golden
- `tests/golden/minimal/gt_rate_postZ_M0.txt` remains the conda-forge/bioconda GCC14 packaged-build baseline.
- This file is a local development baseline for the macOS arm64 Pixi/Clang toolchain, not a replacement for the packaged-build golden.
