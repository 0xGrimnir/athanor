# Third-party notices

Athanor builds on open-source scientific software and public data. Each component keeps its own license. **Verified** means the license was read from the package's own registry metadata or license file; **to verify** means it must be confirmed from the component's license file when it is integrated.

Please cite the underlying engines, not just Athanor, when results contribute to published work.

## Application framework and visualization

| Component | Use | License | Status |
|---|---|---|---|
| Tauri | Desktop shell, installer, updater | Apache-2.0 OR MIT | Verified |
| tauri-plugin-dialog, tauri-plugin-fs | Native file open/save dialogs | Apache-2.0 OR MIT | Verified |
| rfd | Native dialogs on Windows (used by tauri-plugin-dialog) | MIT | Verified |
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
| gargammel 1.1.4 (with fragSim, deamSim, adptSim) | Ancient-DNA read simulation (Linux engines pack) | GPL-3.0-only | Verified (conda metadata) |
| ART 2016.06.05 | Illumina read simulation, used by gargammel | GPL-3.0-only | Verified (conda metadata) |
| AdapterRemoval 3.0.2 | Adapter trimming and read merging | GPL-3.0-or-later | Verified (conda metadata) |
| BWA 0.7.19 | Ancient-DNA read mapping (aln) | GPL-3.0-or-later | Verified (conda metadata) |
| samtools 1.24, htslib 1.24 | BAM processing | MIT | Verified (conda metadata) |
| mapDamage2 2.2.2 | Damage estimation (uses R 4.3.3, GPL-2.0-or-later) | MIT | Verified (conda metadata) |
| minimap2 2.31 | Aligning the truth genome to references for scoring | MIT | Verified (conda metadata) |
| pysam 0.24.1 | Reading BAM files in the science core | MIT | Verified (PyPI) |
| Bowtie2 | Read mapping (planned) | GPL-3.0 | To verify |
| wgsim | Read simulation | MIT | To verify |
| bedtools | Interval operations | MIT | To verify |
| MAFFT | Multiple sequence alignment | BSD | To verify |
| MUSCLE | Multiple sequence alignment | Depends on version | To verify |
| FastTree | Phylogenetic trees | GPL-2.0-or-later | To verify |
| IQ-TREE | Phylogenetic trees with model selection | GPL-2.0 | To verify |
| BLAST+ | Sequence search | Public domain (US government work) | To verify |
| Primer3 (primer3-py) | Primer design | GPL-2.0-or-later | Verified |
| OpenMM 8.6 (Eastman et al., J. Phys. Chem. B 128, 109, 2024) | DNA Structure: backbone relaxation and molecular dynamics, with the AMBER OL15 DNA force field and GBn2 implicit solvent files it ships | MIT (library and force-field files); the CUDA and OpenCL platform plugins are LGPL-3.0 | Verified (conda metadata) |
| Deep DNAshape (Li, Chiu & Rohs, Nature Commun. 15, 1243, 2024; commit fe12e24) | DNA Structure: trained network weights for DNA shape prediction, exported from the authors' checkpoints in CI and shipped with the licence (`core/athanor_core/data/dna/dnashape/LICENSE-DeepDNAshape.txt`); the network itself is reimplemented in NumPy | BSD-3-Clause | Verified |
| simple-dftd3 / dftd3-python 1.4 | D3(BJ) dispersion for the DFT interaction energies | LGPL-3.0-or-later | Verified (conda metadata) |
| AmberTools (tleap, sander) | Independent check in CI that Athanor's OpenMM set-up gives AMBER's energies (not shipped) | GPL-3.0 / LGPL / public domain components | Verified |
| TensorFlow 2.15 | Running the original Deep DNAshape in CI to check Athanor's port (not shipped) | Apache-2.0 | Verified |
| Biopython | Sequence utilities; independent cross-check of the melting temperatures in CI (Tm_NN) | Biopython License Agreement / BSD-3-Clause | Verified |
| 3DNA 2.4.9 (Lu & Olson, Nucleic Acids Res. 31, 5108, 2003; Nat. Protoc. 3, 1213, 2008) | DNA Structure: the standard base geometries and the B- and A-DNA nucleotide sets ship with Athanor (with 3DNA's copyright notice); the program itself is the reference the analysis and rebuild are checked against in CI (not shipped) | BSD-3-Clause | Verified |
| parasail | Pairwise alignment | BSD | Verified |
| AutoDock Vina | Docking | Apache-2.0 | To verify |

## Data sources

| Source | Use | Terms | Status |
|---|---|---|---|
| GenBank / NCBI | Mitochondrial genome sequences | Public, no general restrictions; some individual records may carry their own | To verify |
| GetOrganelleDB | Curated mitochondrial sequence and gene-label database (prototype only) | GPL-3.0 | Verified |
| MMFF94 validation suite and parameter tables (T. A. Halgren) | MMFF94 correctness test; fitted-parameter coverage check | Distributed via CCL and RDKit (BSD-3-Clause) | Verified |
| NIST CCCBDB (SRD 101), NIST Chemistry WebBook (SRD 69) | Experimental reference values for benchmarks | US government data | Verified |
| NCBI RefSeq mitogenomes NC_007596.2, NC_000934.1, NC_005129.2, NC_012920.1 | Ancient DNA Lab references (fetched at pack build; hashes in the pack manifest) | NCBI: no general restrictions; individual records may carry their own | Verified source |
| Allen Ancient DNA Resource v66.p1 (1240K) | Personal Genome reference data: ancient and present-day genotypes, dates, places and published haplogroups (summaries only; the raw data is not redistributed) | CC0 1.0 | Verified |
| PhyloTree Build 17 (RSRS), as packaged for Haplogrep 3 | mtDNA haplogroup tree and mutation weights | MIT | Verified |
| Y-chromosome tree and GRCh37 marker positions, as packaged with Yleaf | Y haplogroup calling | Yleaf is GPL-3.0; the underlying YFull tree's own terms are **to verify before any commercial release** | To verify |
| UCSC hg38ToHg19 liftOver chain | Converting files on GRCh38 to GRCh37 | Free for academic and non-profit use; commercial use of UCSC tools requires a licence — **to verify before any commercial release** | To verify |
| ADMIXTOOLS, EIGENSOFT | Independent cross-check of qpAdm in CI (not shipped) | GPL-3.0 | Verified |
| msprime | Simulated genomes with known answers for the CI gates (not shipped) | GPL-3.0 | Verified |
| ENA PRJEB103048 reads and GenBank PZ244049–PZ244060 (Hohle Fels mammoth ivory, Scientific Reports 2026) | CI reproduction test only: downloaded at test time, never redistributed | INSDC public data | Verified source |
| Ensembl (release 110) | Gene coordinates and transcript models | No restrictions on use | To verify |
| Ensembl GRCh37 REST (variation, VEP) and 1000 Genomes phase 3 frequencies | Personal Genome trait table: variant positions, alleles, the amino acid each allele codes, population frequencies (resolved in CI; only the resulting facts ship) | EMBL-EBI Terms of Use; 1000 Genomes data are open | To verify |
| GWAS Catalog (REST API v2) | Effect alleles, directions and effect sizes recorded for some trait variants (resolved in CI) | EMBL-EBI Terms of Use | To verify |
| NCBI PubMed (E-utilities) | Titles, journals and DOIs of the trait reports' references | NCBI: no general restrictions on bibliographic data | Verified source |
| ADMIXTURE 1.3.0, PLINK 1.9 | Independent cross-checks of the ancestry estimator and of runs of homozygosity in CI (not shipped) | ADMIXTURE: free for use, own licence; PLINK 1.9: GPL-3.0 | Verified |
| CPIC database (api.cpicpgx.org) | Medicines reports: allele definitions, allele functions and frequencies, phenotype tables, prescribing recommendations for level-A gene–drug pairs (downloaded in CI; the pack records the access date and database release) | CC0 1.0. CPIC asks for attribution, for its content to be described as subject to updates, and that its logo and acronym not be used in advertising without NIH permission | Verified |
| 1000 Genomes high-coverage phased panel (NYGC; Byrska-Bishop et al., Cell 2022) | Genotypes at the CPIC positions, used to measure how far a call from a chip's positions can be trusted, and for the PharmCAT comparison | 1000 Genomes data are open (Fort Lauderdale principles; no restrictions on use) | Verified source |
| Ensembl REST (GRCh38 sequence, GRCh37 variation) | Placing the CPIC indels and checking every CPIC reference allele against the genome; checking the GRCh37 liftover (in CI) | EMBL-EBI Terms of Use | To verify |
| RCSB Protein Data Bank (search and download services) | DNA Structure: crystal duplexes from which the A- and Z-DNA parameters, the B-DNA base-pair parameters and the Z-DNA backbone are measured (in CI; only the resulting medians ship) | PDB data are CC0 1.0 | Verified |
| Olson et al., PNAS 95, 11163 (1998), Table 1 | B-DNA dimer step parameters (published values, cited) | Published scientific data | Verified source |
| SantaLucia, PNAS 95, 1460 (1998); Owczarzy et al., Biochemistry 43, 3537 (2004) | Nearest-neighbour parameters and salt correction (published values, cited) | Published scientific data | Verified source |
| TfReg data files (Owczarzy et al. 2004 measured melting temperatures), from the openSUSE source mirror | CI check of the melting temperatures only: downloaded at test time, never redistributed | GPL | Verified |
| Ensembl GRCh37 REST (sequence/region) | DNA Structure: reference sequence 10 bp either side of each trait variant, checked against the reference allele (in CI; the short sequences ship) | EMBL-EBI Terms of Use | To verify |
| RCSB PDB entries 1KX5, 6HKT, 1CDW, 1CGP, 1IHF, 1AAY, 1TSR, 1EFA | DNA Structure scales: the nucleosome core particle, the six-nucleosome fibre and the protein–DNA complexes (coordinates reduced to what the views draw; cited in the app) | PDB data are CC0 1.0 | Verified |
| Dip-C single-cell 3D genomes of GM12878 (Tan et al., Science 361, 924, 2018; GEO GSE117876) | DNA Structure nucleus view: one cell coarse-grained to 1 Mb ships; six cells are used for the CI check | GEO data: NCBI places no restrictions; cited | Verified source |
| Ensembl GRCh37 REST (overlap) | Protein-coding gene density per chromosome for the nucleus check (in CI) | EMBL-EBI Terms of Use | To verify |
| Deep DNAshape processed data (Zenodo 10403307) | CI report only: agreement of the predictions with the authors' published data (not shipped) | CC BY 4.0 | Verified |
| Psi4 S22 database (Jurečka et al. 2006 geometries; Marshall et al. 2011 CCSD(T)/CBS energies) | CI benchmark of the interaction energies (not shipped) | LGPL-3.0 (Psi4); published data | Verified |
| gargammel 1.1.4 deamSim | CI check of the ancient-damage formulas (also in the Linux engines pack) | GPL-3.0-only | Verified (conda metadata) |
| PharmCAT 3.4.0 | Independent cross-check of the pharmacogenomic calls in CI (not shipped) | MPL-2.0 | Verified |

## Compute pack note

The optional GPU compute pack is a separate download that bundles NVIDIA's CUDA runtime libraries from PyPI. They are proprietary; NVIDIA's CUDA EULA lists these runtime components as redistributable with applications, but the exact terms must be checked against the current EULA before Athanor is released publicly.

## License compatibility note

Several engines are GPL-licensed (gargammel, ART, AdapterRemoval and BWA in the Linux engines pack; Bowtie2, FastTree, IQ-TREE, Primer3 planned). In the pack they are unmodified separate programs, run as subprocesses and distributed with their conda package metadata and the curated mitochondrial data comes from GetOrganelleDB (GPL-3.0). While the repository is private and nothing is distributed, this has no effect. Before any public release, Athanor's own license and the way these components are bundled (as separate programs rather than linked libraries) must be reviewed together. This note is not legal advice.

## 3DNA licence

The standard base geometries and the B- and A-DNA nucleotide sets in `core/athanor_core/data/dna/` come from 3DNA v2.4 (https://github.com/x3dna-dssr/x3dna-v2.4). The licence is reproduced here and in `core/athanor_core/data/dna/LICENSE-3DNA.txt`:

```
BSD 3-Clause License

Copyright (c) 1999-2026, Xiang-Jun Lu and Wilma K. Olson
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

## Deep DNAshape licence

The trained weights in `core/athanor_core/data/dna/dnashape/` come from Deep DNAshape (https://github.com/JinsenLi/deepDNAshape). Its licence is reproduced in `core/athanor_core/data/dna/dnashape/LICENSE-DeepDNAshape.txt`:

```
BSD 3-Clause License

Copyright (c) 2023, JinsenLi

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```
