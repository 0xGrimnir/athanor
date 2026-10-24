# Athanor

**A desktop science workbench for Windows: computational chemistry, genomics, ancient-DNA reconstruction and local-only DNA report analysis, built on real engines and validated against published results.**

An *athanor* was the alchemist's furnace, built for slow, controlled transformation. This one runs on real math: every number comes from an established scientific engine, carries a label saying which method produced it, and ships with its known error bars.

> **Status: pre-alpha.** Working: the molecule workbench (MMFF94, GFN2-xTB and DFT tiers, frequencies and IR spectra), the Ancient DNA Lab (simulated and real reads), and Personal Genome (deep and modern ancestry, appearance, food and drink, fitness, medicines). The browser prototype in [`prototype/`](prototype/) is kept as a design reference only. See [Roadmap](docs/ROADMAP.md).

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
| **Molecule** | Build or import molecules; optimize geometry; partial charges, dipoles, harmonic frequencies, IR spectra, animated normal modes | Tiered: MMFF94 (RDKit) → GFN2-xTB (xtb) → DFT/HF (Psi4); optional DFT on NVIDIA GPUs (GPU4PySCF) | Working (0.1–0.2) |
| **Ancient DNA Lab** | Reconstruct mitogenomes from your own ancient reads (FASTQ or BAM, processed locally): circular mapping, duplicate removal, iterative indel-aware consensus with per-base posteriors, authentication by terminal damage and fragment length, human contamination by competitive mapping. Simulated mode: damaged reads from a real mitogenome, every base scored against the known truth | AdapterRemoval, BWA aln, samtools, mapDamage2, minimap2, gargammel (Linux engines pack) | Working (0.3–0.4) |
| **Personal Genome** | Load an AncestryDNA, 23andMe, MyHeritage, FamilyTreeDNA or Living DNA raw file (or a VCF) locally and get evidence-graded reports: ancient ancestry (hunter-gatherer, farmer, steppe), Neanderthal and Denisovan, haplogroups and the ancient people who shared them, continental and regional ancestry, chromosome painting, Indigenous American ancestry, runs of homozygosity, X-chromosome ancestry, appearance, food and drink, fitness, and medicines (CPIC pharmacogenomics for 19 genes); relatives and opt-in health to come | AADR (ancient and present-day genomes), PhyloTree, YFull; trait variants checked against Ensembl, 1000 Genomes and the GWAS Catalog; CPIC, checked against PharmCAT | Working (0.5–0.7) |
| **DNA Structure** | Any sequence as a 3D helix, atom by atom, in B, A or Z form, from parameters measured on crystal structures, with the backbone closed by OpenMM. Molecular dynamics (AMBER OL15, on the GPU when available). Base-pair, step and groove parameters for any duplex. Shape in solution predicted from sequence (Deep DNAshape). Pairing and stacking energies (GFN2-xTB, DFT). Melting. Ancient damage on a molecule. A zoom from helix to nucleosome, protein–DNA complexes, chromatin fibre and a whole nucleus. Personal Genome variants shown in the helix and the nucleus | 3DNA's methods (checked against 3DNA 2.4.9); OpenMM; Deep DNAshape; xtb, Psi4; PDB, Dip-C | Working (0.8–0.9) |
| **Deep Time** | Ancient nuclear genomes, traits from ancient genomes, ancient proteins, DNA decay, radiocarbon calibration | minimap2, bwa-mem2, samtools; IntCal20 | Planned |
| **Populations** | Evolution sandbox: drift, selection, bottlenecks and interbreeding, compared with real ancient data | msprime, SLiM | Planned |
| **RNA and proteins** | RNA folding; Personal Genome variants placed on AlphaFold and PDB protein structures | Mol* | Planned |

## Validation (Molecules)

| Tier | Correctness check | Accuracy check against experiment |
|---|---|---|
| MMFF94 | Reproduces Halgren's MMFF94 validation suite (755 molecules to 1×10⁻⁴ kcal/mol) | Declines molecules outside MMFF94's fitted bond parameters |
| GFN2-xTB | The independent tblite library reproduces the xtb program; frequencies match xtb's own Hessian and ASE's vibrational analysis | Bond lengths, dipoles and IR fundamentals vs NIST data, compared with published GFN2-xTB statistics |
| DFT | Matches Psi4's own frequency analysis; the GPU tier's code path (PySCF) matches Psi4 | Same reference set at B3LYP/cc-pVTZ, compared with published B3LYP statistics (Janoschek 2001; CCCBDB scale factors; Hait & Head-Gordon 2018) |

Reference values and sources: [`core/benchmarks/reference/phase1.json`](core/benchmarks/reference/phase1.json). Run locally with `python -m athanor_core.benchmark`.

## Validation (DNA Structure)

Every check runs in CI on each change (`tools/dnadata/`, `.github/workflows/dna*.yml`), and CI confirms that the shipped data files regenerate exactly.

**Helix geometry.** Checked against 3DNA 2.4.9 and the DNA-only crystal duplexes in the PDB (X-ray, 1.6 Å or better; entry list pinned in [`tools/dnadata/entries.txt`](tools/dnadata/entries.txt)).

| Check | Result |
|---|---|
| Base-pair, step and groove parameters vs 3DNA's `analyze` | Equal to its printed precision on 98 of 101 duplexes (3DNA pairs fewer bases in the other 3) |
| Rebuild vs 3DNA's `rebuild` from the same parameters | 0.001 Å RMSD (15 duplexes) |
| Each crystal rebuilt from its own parameters | Bases back to a median 0.03 Å |
| B-DNA built from sequence alone (Olson et al. 1998 steps) | Median 1.42 Å all-atom RMSD to 52 crystals |
| A- and Z-DNA from parameters measured on the other crystals | Median 1.56 Å (16 A-DNA) and 0.57 Å (28 Z-DNA) |
| Backbone closed by OpenMM (OL15, GBn2) | Every O3′–P link 1.55–1.67 Å in all 96 models (rigid build 0.68–2.90 Å); bases move at most 0.3 Å; models no further from the crystals |

**Physics.**

| Check | Result |
|---|---|
| OpenMM set-up vs AMBER (tleap + sander, same coordinates) | @@MDENERGY@@ |
| Dickerson dodecamer in dynamics (@@MDPS@@ ps, implicit water) vs its NMR structure (1NAJ) | @@MDRESULT@@ |
| Pairing and stacking energies on S22 vs CCSD(T)/CBS | Mean error 0.76 kcal/mol (GFN2-xTB) and 0.26 (B3LYP-D3(BJ)/def2-TZVP, counterpoise); adenine–thymine pair −15.9 / −17.0 vs −16.7, stacked −12.3 / −11.6 vs −11.7 |
| In a built helix | G–C pairs bind more strongly than A–T (−30.8 vs −14.2 kcal/mol, GFN2-xTB; DFT −32.1 and −14.8) |
| Melting temperature vs Owczarzy et al. 2004 (92 duplexes × 5 Na⁺ concentrations) | Mean error 1.66 °C; identical to Biopython's Tm_NN |
| Deep DNAshape port vs the original TensorFlow model | Same predictions to within 1e-6 for all 18 features, at 2, 4 and 7 bp of flank@@SHAPE@@ |
| Ancient damage vs gargammel's deamSim | @@DAMAGE@@ |

**Scales.**

| Check | Result |
|---|---|
| Nucleosome 1KX5 rebuilt from its own parameters | Bases back in place to 0.02 Å; left-handed superhelix of radius 41.8 Å and pitch 25.9 Å (Richmond & Davey 2003: 41.9 and 25.9) |
| DNA bends in protein complexes vs the papers | TBP 70° (80°), CAP 81° (90°), IHF 164° (more than 160°) |
| Six-nucleosome fibre 6HKT | Two-start zigzag: nucleosomes 6.1 nm from the one two along, 21.4 nm from their neighbour |
| Dip-C single-cell genomes (6 GM12878 cells) | Gene-rich chromosomes sit centrally: chromosome 19 inside 18 in 5 of 6 cells; radial position vs gene density ρ = −0.74 |

Limits the app states: dynamics use implicit water (no explicit ions or water molecules). Energies are for the bases alone in vacuum. Shape predictions are for free DNA. The nucleosome view threads your sequence onto a crystal's path; it does not predict where nucleosomes form. The nucleus is one real cell, not your own.

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
| Medicines against PharmCAT 3.4.0 | On the same 1000 Genomes high-coverage genotypes of 500 people (100 per continental group), the same diplotype and phenotype in at least 99.8% of calls for every gene compared (8,472 calls) |
| Medicines, chip input | 23andMe-style files written from 1000 Genomes genotypes (GRCh37, both strands, indels as I/D) read back exactly: 218,770 genotypes |
| Medicines, reliability figure | The agreement figure each gene card shows, estimated on half of 1000 Genomes, holds on the other half to within 0.007 on average |

These run as gated checks in CI whenever the reference data is rebuilt (`tools/pgdata/validate.py`, `tools/pgdata/crosscheck_qpadm.py`, `tools/pgdata/validate_pgx.py`).

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

Athanor is open source under the [Apache License 2.0](LICENSE). Copyright © 2026 Tyler Lee Jones.

If you redistribute Athanor or build on it, keep the [`NOTICE`](NOTICE) file with it, as the license requires. Third-party components keep their own licenses; see [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md).
