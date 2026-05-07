GT Golden Provenance: macOS arm64 Pixi/Clang

Golden file
- `tests/golden/minimal/gt_rate_postZ_M0.osx-arm64-clang-pixi.txt`

Status
- Recorded on `2026-05-07` after the one-public-`SEED` RNG contract change.
- Golden SHA256: `4d4fd45b61590089d6da4b6d60ca1d41d1ebe7fb089d7a0154744fb200f3038a`

What matches this golden
- Local macOS arm64 Pixi environment in this repository.
- `PhyloAcc-GT` built by `pixi run build`.
- Compiler: `clang version 19.1.7`, target `arm64-apple-darwin20.0.0`.

Build provenance
- Platform: `darwin`, `arm64`.
- Pixi lock SHA256: `8a2d224c7e4dc8601b9cc52f90a7467356b92ff111f40e16311562ef2c9de623`
- `PhyloAcc-GT` binary SHA256: `3a7c44e14b18d0ad2fdc39db855dee0b9860fb1ca6d408ab374b733faa4cb467`
- Compiler executable: `.pixi/envs/default/bin/clang++`
- Compiler version output:
  - `clang version 19.1.7`
  - `Target: arm64-apple-darwin24.6.0`

Repeatability checks
- Recording run: `PHYLOACC_RUN_GT=1 PHYLOACC_RUN_TESTDATA=1 PHYLOACC_RECORD_GOLDEN=1 pytest -q tests/test_st_gt.py::test_gt_minimal_run tests/test_optional_testdata.py::test_optional_testdata_st_golden`
  - Result: `2 passed in 107.22s`
- RNG contract run: `PHYLOACC_RUN_RNG=1 PHYLOACC_RUN_GT=1 pytest -q tests/test_rng_repro.py`
  - Result: `6 passed in 548.70s`

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
