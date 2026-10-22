# Athanor

**A desktop science workbench for Windows: computational chemistry, genomics, ancient-DNA reconstruction and local-only DNA report analysis, built on real engines and validated against published results.**

An *athanor* was the alchemist's furnace, built for slow, controlled transformation. This one runs on real math: every number comes from an established scientific engine, carries a label saying which method produced it, and ships with its known error bars.

> **Status: pre-alpha.** Phase 1 (the molecule workbench) is working: MMFF94, GFN2-xTB and DFT tiers, frequencies and IR spectra, provenance, and an in-app Validation view. The browser prototype in [`prototype/`](prototype/) is kept as a design reference only; its chemistry and gap-fill modules were illustrative and are superseded. See [Roadmap](docs/ROADMAP.md).

---

## Principles

1. **Real science, not set dressing.** Results come from established engines (RDKit, xtb, Psi4, OpenMM, samtools, minimap2 and others), never from hand-tuned approximations presented as computation.
2. **Every number carries its provenance.** The engine, version, method, parameters and random seed are recorded and exportable for any result.
3. **Uncertainty is shown, not hidden.** Confidence intervals, posteriors and method error bars are drawn in the interface alongside the values.
4. **Validated against literature.** Each module ships a benchmark suite that reproduces published reference values. If a benchmark fails, the module doesn't ship.
5. **Limits are labeled.** Where something is impossible (predicting an organism from a genome, recovering bases that were never sequenced), the app says so plainly.
6. **Your genetic data stays on your machine.** DNA reports are parsed locally and are never uploaded or sent anywhere.

## Modules

| Module | What it does | Engines | State |
|---|---|---|---|
| **Molecule** | Build or import molecules; optimize geometry; partial charges, dipoles, harmonic frequencies, IR spectra, animated normal modes | Tiered: MMFF94 (RDKit) → GFN2-xTB (xtb) → DFT/HF (Psi4); optional DFT on NVIDIA GPUs (GPU4PySCF) | Phase 1, working |
| **Ancient DNA Lab** | Reconstruct mitogenomes from your own ancient reads (FASTQ or BAM, processed locally): circular mapping, duplicate removal, iterative indel-aware consensus with per-base posteriors, authentication by terminal damage and fragment length, human contamination by competitive mapping. Simulated mode: damaged reads from a real mitogenome, every base scored against the known truth | AdapterRemoval, BWA aln, samtools, mapDamage2, minimap2, gargammel (Linux engines pack) | Phase 2–2b, working |
| **Genome** | Genome browsing: reads, variants, coverage, multiple assemblies; circular and comparative views | igv.js, minimap2, MAFFT | Phase 4 |
| **Sequence Workbench** | Editing, ORFs, restriction maps, primer design, cloning simulation, alignments, phylogenetic trees, CRISPR guide scoring | Primer3, MAFFT/MUSCLE, FastTree/IQ-TREE, BLAST+ | Phase 5 |
| **Structure** | Protein structures, linking each genome position to its residue in 3D, docking | Mol*, AutoDock Vina | Phase 5 |
| **Personal Genome** | Load an AncestryDNA, 23andMe, MyHeritage, FamilyTreeDNA or Living DNA raw file (or a VCF) locally and get about 70 evidence-graded reports: ancient ancestry (hunter-gatherer, farmer, steppe), Neanderthal and Denisovan, haplogroups and the ancient people who shared them, modern ancestry, relatives, appearance, food and drink, pharmacogenomics, and opt-in health | AADR, 1000 Genomes, HGDP, PhyloTree, YFull; CPIC, PGS Catalog, ClinVar later | Phase 3, first release working |

## Validation (Phase 1)

| Tier | Correctness check | Accuracy check against experiment |
|---|---|---|
| MMFF94 | Reproduces Halgren's MMFF94 validation suite (755 molecules to 1×10⁻⁴ kcal/mol) | Declines molecules outside MMFF94's fitted bond parameters |
| GFN2-xTB | The independent tblite library reproduces the xtb program; frequencies match xtb's own Hessian and ASE's vibrational analysis | Bond lengths, dipoles and IR fundamentals vs NIST data, compared with published GFN2-xTB statistics |
| DFT | Matches Psi4's own frequency analysis; the GPU tier's code path (PySCF) matches Psi4 | Same reference set at B3LYP/cc-pVTZ, compared with published B3LYP statistics (Janoschek 2001; CCCBDB scale factors; Hait & Head-Gordon 2018) |

Reference values and sources: [`core/benchmarks/reference/phase1.json`](core/benchmarks/reference/phase1.json). Run locally with `python -m athanor_core.benchmark`.

## Validation (Ancient DNA Lab)

Woolly mammoth mtDNA (NC_007596.2) simulated at 10× with Briggs damage (ν 0.03, λ 0.4, δd 0.01, δs 0.3), mapped to African (NC_000934.1) and Asian (NC_005129.2) elephant mitogenomes with three mapping iterations:

| Check | Result |
|---|---|
| Damage recovery (mapDamage2) | δd 0.0089–0.0109 vs 0.01; λ 0.377–0.389 vs 0.4 (δs differs by parameterisation; reported, not used) |
| Reference bias (African elephant) | mapped reads 2,795 → 3,763; divergence estimate 3.99% → 4.63% (true 4.76%) |
| Consensus vs truth | 99.951% correct against either elephant reference |
| Posterior calibration | expected calibration error 0.0004; residual overconfidence at single-read sites is reported |

These run as gated checks in CI on every pack build (`compute/ci_adna.py`).

### Real data: reproduction of published mitogenomes

The deposited reads of seven woolly mammoth ivory samples from Hohle Fels (Scientific Reports, 2026; reads ENA PRJEB103048, sequences GenBank PZ244049–PZ244060) are run through Athanor from the mammoth reference and compared with the published sequences:

| Sample | GenBank | Mean depth | Bases both call | Disagreements |
|---|---|---|---|---|
| HOLF001.A | PZ244049.1 | 19.6× | 16,538 | 0 |
| HOLF003.A | PZ244050.1 | 28.2× | 16,721 | 0 |
| HOLF009.A | PZ244054.1 | 33.7× | 16,767 | 0 |
| HOLF010.B | PZ244055.1 | 36.4× | 16,694 | 0 |
| HOLF022.A | PZ244058.1 | 124.1× | 16,767 | 0 |
| HOLF024.B | PZ244059.1 | 28.1× | 16,538 | 1 |
| HOLF025.B | PZ244060.1 | 16.2× | 16,767 | 0 |

Every sample shows the terminal C→T excess of authentic ancient DNA. The one disagreement (HOLF024.B, position 12,196) is a base the published sequence gives as T, while all six deposited reads covering it read C. The remaining differences are lengths in the control-region repeat (HOLF003.A, HOLF010.B), where the published sequences also carry uncalled bases.

Simulations with known answers check the rest. Contamination estimates were 0.00%, 1.95% and 10.40% against true values of 0%, 1.95% and 10.20%. The mammoth sequence, rebuilt starting from the African elephant genome (about 9% divergent), agreed with the truth at every one of 16,462 called bases (latest CI run; the simulations are random, so the figures vary slightly between runs). All of it runs in CI (`compute/ci_reproduce.py`).

## Validation (Personal Genome)

Everything the module reports is checked before release, and each report card shows what was checked.

| Check | Result |
|---|---|
| qpAdm against ADMIXTOOLS 8.0.2 | On a well-powered simulation, weights agree to 0.0006 and standard errors to 0.0003; on real population data, within two standard errors |
| qpAdm accuracy | Recovers simulated mixture proportions within 0.03 when the truth is known |
| Raw-file parsers | Files written in the AncestryDNA, 23andMe, MyHeritage, FamilyTreeDNA, Living DNA and VCF layouts (including a GRCh38 VCF, converted by liftOver) read back with every genotype exact |
| PCA placement | Genomes cut down to a genotyping chip's worth of markers land a median ~1% of the map's width from where the same genomes land with all 1.15 million markers |
| Neanderthal ancestry | On held-out genomes reduced to chip size, the 95% interval contains the all-marker estimate in ≥ 90% of cases |
| Y haplogroups | 586 of 600 ancient men on the same lineage as the published call, 204 identical |
| mtDNA haplogroups | 1,476 of 1,500 ancient mitogenomes on the same lineage, 1,246 identical |

These run as gated checks in CI whenever the reference data is rebuilt (`tools/pgdata/validate.py`, `tools/pgdata/crosscheck_qpadm.py`).

## Architecture

- **Shell:** [Tauri 2](https://tauri.app). A native Windows app with a web-technology interface, shipped as a standard installer with signed in-app updates.
- **Science core:** Python package [`core/athanor_core`](core/athanor_core). RDKit, the xtb program and Psi4 run natively on Windows; every result carries a provenance record. The shell talks to it over stdin/stdout (JSON-RPC); no network port is opened.
- **Genomics engines:** WebAssembly builds from [biowasm](https://github.com/biowasm/biowasm) (samtools, bcftools, minimap2, Bowtie2, MAFFT and others) running inside the app. No Linux layer required.
- **Optional Linux engines pack:** a private WSL2 distribution that the app installs and verifies (SHA-256 per part). It carries GPU4PySCF (self-tested against its CPU path before the DFT·GPU tier is offered) and the ancient-DNA tools (bioconda), plus reference genomes fetched from NCBI with recorded accession versions and hashes. No GPU is needed for the ancient-DNA tools.

## System requirements

- Windows 10 or 11, 64-bit
- NVIDIA GPU recommended for molecular dynamics (CUDA). Everything else runs on the CPU.
- 16 GB RAM recommended

## Install and update

Download the latest installer from [Releases](../../releases) and run it. The app isn't code-signed yet, so Windows SmartScreen will warn about an unknown publisher: choose **More info → Run anyway**. After that, Athanor checks for updates at launch and installs them in one click.

## Known limitations

- **DFT cost** grows roughly with the cube of system size. Large systems use MMFF94 or GFN2-xTB, and the tier is always labeled.
- **The WebAssembly genomics engines** are slower than native builds and capped at about 4 GB of memory. That's fine for mitochondrial, bacterial and targeted-region work; whole-human-genome mapping waits for the optional compute pack.
- **No simulation predicts an organism from its genome.** Splice experiments stop at molecular compatibility.
- **Ancient-DNA reconstruction** gives probabilities for missing bases, never certainty.
- **Contamination by competitive mapping** detects human DNA in non-human samples only, and is a lower bound. Human ancient samples need a haplotype-based method (planned).

## Citing

If Athanor contributes to your work, please cite it using [`CITATION.cff`](CITATION.cff), and cite the underlying engines listed in [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md).

## License

Copyright © 2026 Tyler Lee Jones. All rights reserved for now; see [`LICENSE`](LICENSE). Third-party components keep their own licenses.
