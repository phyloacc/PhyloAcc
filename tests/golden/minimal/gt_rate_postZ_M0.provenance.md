GT Golden Provenance

Golden file
- `tests/golden/minimal/gt_rate_postZ_M0.txt`

Status
- Reconfirmed on `2026-03-26`.
- Golden SHA256: `dbcee113d6ffe33143d726844e2049ae7fbfd4b8b7f3626e6e8d3936c4542d05`

What matches this golden exactly
- Bioconda packaged `phyloacc 2.4.3` `PhyloAcc-GT`
- Bioconda packaged `phyloacc 2.4.5` `PhyloAcc-GT`
- Local build of `v2.4.5` source using the conda-forge GCC 14 toolchain reconstructed under `dev/phyloacc-buildcheck4`

What did not match
- Earlier local source rebuilds done through the system compiler path (`g++ 8.5`)

Conclusion
- This GT golden is a packaged-build baseline.
- It should be compared against binaries built with the conda-forge/bioconda compiler stack, not against arbitrary local compiler builds.

Relevant local artifacts from the investigation
- GCC 14 wrapper used for the successful local reproduction: `dev/cc14-wrap.sh`
- Matching local source build tree: `dev/phyloacc-2.4.5-src-gcc14/`
- Local conda-build recipe copy: `dev/recipe-phyloacc-245/`
- Local conda-build output package: `dev/conda-bld245/linux-64/phyloacc-2.4.5-py313h4c9e609_1.conda`

Key build-path notes
- Matching local source reproduction required:
  - conda-forge GCC 14 front-end
  - conda-forge sysroot
  - conda-forge host libraries from the repaired prefix env under `dev/phyloacc-buildcheck4`
- The decisive factor was the build toolchain/provenance, not source drift between `2.4.3`, `2.4.4`, and `2.4.5`.
