# Changelog

All notable changes to Athanor are recorded here. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow [Semantic Versioning](https://semver.org/).

## [Unreleased]

## [0.6.0] - 2026-10-23

### Added
- **Personal Genome: modern ancestry.**
  - **Continental and regional ancestry:** supervised admixture against present-day reference genomes in the AADR (1000 Genomes, HGDP and others), with 95% intervals from a block bootstrap and a measured detection limit for every region. Every region keeps its own reference panel; regions that validation could not tell apart are reported together. Each card also shows what reference people from your main region typically get.
  - **Chromosome painting:** the pair of continental ancestries along each of your 22 chromosomes.
  - **Indigenous American ancestry,** with its detection limit and plain statements of what it cannot show (it is not tribal membership; the references are mostly Central and South American).
  - **Runs of homozygosity,** shown against reference groups measured the same way and against the expected values for related parents.
  - **X-chromosome ancestry,** compared with your autosomes.
- **Personal Genome: trait reports (17).** Appearance: eye colour, red hair, freckles, the KITLG blond-hair variant, skin pigmentation variants, earwax type, hair thickness, hair form. Food and drink: lactose tolerance, alcohol flush, caffeine metabolism, bitter taste, vitamin D, vitamin B12 and secretor status, omega-3 and omega-6 conversion, celiac HLA types. Fitness: ACTN3. Each shows your genotype, how common the variant is, the evidence grade and what it cannot tell you.
- Every trait variant's effect allele is checked when the reference data are built: against the amino acid it codes (Ensembl VEP), its frequencies in 1000 Genomes, or the direction recorded in the GWAS Catalog. References are looked up in PubMed.
- The PCA maps now label the reference populations and show your closest ancient individuals.

### Fixed
- The maternal-line card could report a haplogroup supported by a single marker. AncestryDNA reports some mtDNA markers on the opposite strand; those are now read on the right strand, and a haplogroup is reported only when at least 5 markers along its lineage agree. Otherwise the card shows how many of your markers matched.
- The Denisovan card said Denisovan ancestry raises the statistic; it lowers it, and Neanderthal ancestry raises it. The summary now says which reference regions your value matches.
- Haplogroup timelines no longer repeat the same ancient genome on many rows.

### Changed
- Reference data version 4 (download it again from the Personal Genome tab; the app offers the update).
- Runs of homozygosity follow PLINK 1.9's algorithm exactly (checked against PLINK in CI).
- Polygenic reports (height, baldness, body mass, chronotype, sleep, coffee and alcohol intake) move to the imputation step: a raw chip file covers too few of a modern score's variants.

## [0.5.0] - 2026-10-22

### Added
- **Personal Genome** (new tab). Load the raw data file from AncestryDNA, 23andMe, MyHeritage, FamilyTreeDNA or Living DNA, or a VCF, and get a set of reports. The file is read on your computer and never uploaded.
  - **Your file:** testing company, genome build (GRCh38 files are converted), call rate, genetic sex, and how many markers could be used.
  - **Ancient ancestry:** published qpAdm models fitted to your genome (hunter-gatherer, early farmer and steppe herder, with a Near Eastern and a South Asian model as well), each with proportions, 95% intervals and whether the model actually fits.
  - **Your place among people past and present:** principal-component maps with present-day reference populations and ancient genomes, and the ancient individuals closest to you, on a map with their dates and places.
  - **Neanderthal ancestry** as an f4-ratio with an interval, shown against reference groups, plus a Denisovan affinity statistic and the list of Neanderthal-matching variants you carry.
  - **Maternal and paternal lines:** mtDNA and Y haplogroups from the PhyloTree and YFull trees, the ancient genomes that carry the same branch, and the oldest dated genome on each step of your lineage.
  - Every card carries an evidence grade, what it cannot tell you, the method, its sources, and how it was checked. The report saves as JSON.
- Reference data (about 150 MB) is downloaded once from the release feed and checked against its checksum. It is built in CI from public sources with recorded versions: the Allen Ancient DNA Resource (CC0), PhyloTree 17, the YFull tree as packaged with Yleaf, and the UCSC liftOver chain.
- CI gates on the Personal Genome: qpAdm is checked against ADMIXTOOLS 8.0.2 on a well-powered simulation (weights agree to 0.0006, standard errors to 0.0003) and recovers known simulated proportions; the parsers round-trip every vendor layout exactly; chip-sized genomes land where full-marker genomes land on the PCA; the Neanderthal interval covers the all-marker estimate; haplogroup calls match published ones for ancient individuals (Y: 586 of 600 on the same lineage; mtDNA: 1,476 of 1,500).

## [0.4.0] - 2026-10-21

### Added
- **Ancient DNA Lab: your own data.** A new "Your data" mode reconstructs a mitogenome from real reads: paired-end or single-end FASTQ, or BAM files such as those deposited with a publication. Files are chosen with the Windows file dialog and read where they are; nothing is uploaded.
  - Circular mapping (reference elongated by its own start, positions folded back) and per-library duplicate removal.
  - Iterative, indel-aware reconstruction from the species' or a relative's genome until the sequence stops changing; stops and reports calls that alternate between states.
  - Damage measured against the sample's own consensus; reads cut short by the sequencing read length get no 3′ damage term.
  - Authentication: terminal C→T and G→A with 95% intervals against the molecule interior, whole-molecule length distribution, mapDamage2 profiles.
  - Present-day human contamination by competitive mapping (Feuerborn et al. 2020), with an exact binomial interval and the length and damage of the human-assigned reads.
  - Consensus FASTA with an N mask (depth and posterior thresholds), a conventional majority-rule consensus, and a full JSON report with provenance (input SHA-256, settings, every command).
  - Optional comparison with a published sequence, listing every difference with its read support.
- CI reproduces seven published woolly mammoth mitogenomes (Hohle Fels; ENA PRJEB103048, GenBank PZ244049–PZ244060) from their deposited reads: six identical at every base both sequences call; the seventh differs at one base, where the published sequence has T and all six deposited reads have C. Simulations check the contamination estimate (0%, 2%, 10%) and reconstruction starting from the African elephant genome.

### Changed
- No new Linux engines pack is needed; the existing pack (version 2) has every tool.

## [0.3.1] - 2026-10-20

### Added
- Ancient DNA runs show a continuous progress bar: percentage, time left, time elapsed and the current step of the planned run. Progress updates every second inside each step and never stalls; the time estimate learns this machine's speed separately for each kind of step. Replaying CI timings, the estimate is within about 5% from the halfway point.

## [0.3.0] - 2026-10-20

### Added
- **Ancient DNA Lab** (new tab). Simulates ancient reads from a real mitogenome with gargammel (Briggs-model damage, log-normal fragment lengths, ART Illumina errors, optional present-day human contamination), merges pairs with AdapterRemoval, maps with BWA aln (ancient-DNA settings), estimates damage with mapDamage2, and calls a damage-aware consensus whose per-base posteriors are checked for calibration against the known truth. Iterative mapping against a relative's genome reduces reference bias. Every command and tool version is recorded; reference genomes come from NCBI with accession versions and SHA-256 hashes.
- Woolly mammoth vs African and Asian elephant demo, with CI gates on damage recovery, reference-bias reduction, accuracy and calibration.

### Changed
- The optional pack is now the **Linux engines pack** (version 2, about 2.4 GB): it adds the ancient-DNA tools and installs without an NVIDIA GPU. A failed GPU self-test turns off only the GPU tier.

## [0.2.5] - 2026-10-20

### Added
- Analytic Hessians and IR intensities (atomic polar tensors from the coupled-perturbed SCF equations) for the GPU tier, closed-shell molecules. CI checks the formula on CPU against finite differences (frequencies within 1.5 cm⁻¹, intensities within 2%); the GPU self-test checks the GPU path the same way, and the GPU tier uses it only after that check passes on the user's machine. Otherwise frequencies use finite differences as before.
- "GPU ready" offers to re-run the self-test.

## [0.2.4] - 2026-10-20

### Fixed
- **Stop** hung while a GPU calculation was running: the "is the GPU core running" check waited on the lock the running calculation holds.
- Geometry optimization now stops on Gaussian's actual default criteria (max/RMS force and max/RMS displacement); it previously kept iterating to about ten times tighter forces than the stated threshold.

### Added
- DFT optimizations start with a GFN2-xTB pre-optimization (on by default, recorded in the provenance). Ethanol at B3LYP/def2-SVP from an RDKit structure: 50 DFT steps before, 18 after.
- Live progress for optimizations (step, energy, max force against the target) and frequency jobs (gradient n of N).
- The PySCF/GPU tier reuses each step's converged density as the next step's initial guess.

## [0.2.3] - 2026-10-20

### Fixed
- GPU self-test failed with "NotImplementedType object is not callable": GPU4PySCF does not implement Mulliken analysis on GPU objects, so charges are now computed with PySCF's CPU routine from the GPU density on both paths (regression test added). The self-test also compares charges.
- First GPU self-test on real hardware (RTX 5080, driver 610.88) passed: caffeine B3LYP/def2-SVP/RI, GPU vs CPU ΔE 2.2×10⁻¹⁰ Eh, max Δgradient 1.8×10⁻⁶ Eh/bohr, 6.1 s vs 71.3 s.
- Retrying GPU compute no longer re-downloads an already installed compute pack, and self-test errors are reported instead of stopping the test.

## [0.2.2] - 2026-10-19

### Changed
- The installer build now determines which engine files need relocating (16 of about 35,000); first launch touches only those.

## [0.2.1] - 2026-10-19

### Changed
- First-launch engine relocation skips files that cannot change on Windows, reads the rest in parallel, and shows a live file count.

## [0.2.0] - 2026-10-19

### Added
- Optional GPU compute pack: DFT on NVIDIA GPUs through GPU4PySCF, in a private WSL2 distribution the app installs, verifies and self-tests (GPU vs CPU) before offering the DFT·GPU tier.
- PySCF DFT engine; CI checks its CPU path against Psi4.

## [0.1.0] - 2026-10-19

### Added
- Science core (`core/athanor_core`): MMFF94 (RDKit), GFN2-xTB (xtb program) and DFT (Psi4) engines behind one interface; shared harmonic vibrational analysis with IR intensities; provenance on every result; JSON-RPC server.
- Validation: MMFF94 validation suite, cross-checks against the xtb program, Psi4 and ASE, and a literature benchmark with published-statistics gates.
- Desktop UI: molecule workbench (build/import, optimize, frequencies and IR, properties, animated normal modes, error bars, provenance) and a Validation view.
- Windows shell (Tauri 2) with first-launch engine unpacking, signed in-app updates, and a release pipeline that smoke-tests the packed runtime before publishing.
- Repository setup: README, license, security policy, citation metadata, third-party notices and roadmap.
- Browser prototype preserved in `prototype/` as the interface reference.

### Known issues carried from the prototype
- The prototype's molecule force field and ancient-DNA gap-fill were illustrative, not computational. Both are being replaced with real engines in Phases 1 and 2.
