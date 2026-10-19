# Third-party notices

Athanor builds on open-source scientific software and public data. Each component keeps its own license. **Verified** means the license was read from the package's own registry metadata or license file; **to verify** means it must be confirmed from the component's license file when it is integrated.

Please cite the underlying engines, not just Athanor, when results contribute to published work.

## Application framework and visualization

| Component | Use | License | Status |
|---|---|---|---|
| Tauri | Desktop shell, installer, updater | Apache-2.0 OR MIT | Verified |
| three.js | 3D rendering | MIT | Verified |
| fflate | Unzipping report archives | MIT | Verified |
| igv.js | Genome browser | MIT | Verified |
| Mol* (molstar) | Protein structure viewer | MIT | Verified |
| 3Dmol.js | Molecular viewer | BSD-3-Clause | Verified |
| ideogram.js data | Cytogenetic bands and gene positions | CC0-1.0 | Verified |
| Chakra Petch, JetBrains Mono, Saira Semi Condensed (via Fontsource) | Bundled interface fonts | OFL-1.1 | Verified |
| conda-pack | Packs the science runtime into the installer | BSD-3-Clause | Verified |

## Chemistry engines

| Component | Use | License | Status |
|---|---|---|---|
| RDKit | MMFF94 force field, cheminformatics | BSD-3-Clause | Verified |
| OpenMM | Molecular dynamics (CUDA) | MIT, with LGPL parts | To verify |
| tblite (Grimme lab) | Independent GFN2-xTB implementation for CI cross-checks (not shipped) | LGPL-3.0-or-later | Verified |
| xtb (Grimme lab) | GFN2-xTB engine | LGPL-3.0-or-later | Verified |
| ASE | Independent vibrational-analysis cross-check | LGPL-2.1-or-later | Verified |
| Psi4 | DFT and Hartree-Fock | LGPL-3.0 | To verify |
| PySCF | DFT engine of the GPU tier (CPU path checked against Psi4 in CI) | Apache-2.0 | Verified |
| GPU4PySCF (ByteDance) | GPU-accelerated DFT (optional compute pack) | Apache-2.0 | Verified |
| CuPy | GPU arrays for GPU4PySCF (compute pack) | MIT | Verified |
| NVIDIA CUDA runtime libraries (cuBLAS, cuSOLVER, cuSPARSE, cuFFT, cuRAND, NVRTC, cuTENSOR; PyPI wheels) | GPU math libraries (compute pack) | NVIDIA proprietary license; redistribution terms to verify before public release | To verify |
| Ubuntu 24.04 base image | Root filesystem of the compute pack's WSL2 distribution | Various free licenses (Canonical) | To verify |
| NumPy | Numerics | BSD-3-Clause (with bundled MIT, Zlib, 0BSD, CC0 parts) | Verified |
| SciPy | Numerics | BSD-3-Clause | To verify |

## Genomics engines

| Component | Use | License | Status |
|---|---|---|---|
| biowasm / Aioli | Runs genomics tools as WebAssembly | MIT | Verified |
| samtools, bcftools, htslib | Read and variant processing | MIT/Expat | To verify |
| minimap2 | Sequence alignment | MIT | To verify |
| Bowtie2 | Ancient-DNA read mapping | GPL-3.0 | To verify |
| wgsim | Read simulation | MIT | To verify |
| bedtools | Interval operations | MIT | To verify |
| MAFFT | Multiple sequence alignment | BSD | To verify |
| MUSCLE | Multiple sequence alignment | Depends on version | To verify |
| FastTree | Phylogenetic trees | GPL-2.0-or-later | To verify |
| IQ-TREE | Phylogenetic trees with model selection | GPL-2.0 | To verify |
| BLAST+ | Sequence search | Public domain (US government work) | To verify |
| Primer3 (primer3-py) | Primer design | GPL-2.0-or-later | Verified |
| Biopython | Sequence utilities | Biopython License Agreement / BSD-3-Clause | Verified |
| parasail | Pairwise alignment | BSD | Verified |
| AutoDock Vina | Docking | Apache-2.0 | To verify |

## Data sources

| Source | Use | Terms | Status |
|---|---|---|---|
| GenBank / NCBI | Mitochondrial genome sequences | Public, no general restrictions; some individual records may carry their own | To verify |
| GetOrganelleDB | Curated mitochondrial sequence and gene-label database (prototype only) | GPL-3.0 | Verified |
| MMFF94 validation suite and parameter tables (T. A. Halgren) | MMFF94 correctness test; fitted-parameter coverage check | Distributed via CCL and RDKit (BSD-3-Clause) | Verified |
| NIST CCCBDB (SRD 101), NIST Chemistry WebBook (SRD 69) | Experimental reference values for benchmarks | US government data | Verified |
| Ensembl (release 110) | Gene coordinates and transcript models | No restrictions on use | To verify |

## Compute pack note

The optional GPU compute pack is a separate download that bundles NVIDIA's CUDA runtime libraries from PyPI. They are proprietary; NVIDIA's CUDA EULA lists these runtime components as redistributable with applications, but the exact terms must be checked against the current EULA before Athanor is released publicly.

## License compatibility note

Several engines are GPL-licensed (Bowtie2, FastTree, IQ-TREE, Primer3) and the curated mitochondrial data comes from GetOrganelleDB (GPL-3.0). While the repository is private and nothing is distributed, this has no effect. Before any public release, Athanor's own license and the way these components are bundled (as separate programs rather than linked libraries) must be reviewed together. This note is not legal advice.
