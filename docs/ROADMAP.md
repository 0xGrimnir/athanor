# Roadmap

Athanor is a workbench for the molecules of life and their history. It is organised in layers, each useful on its own and each feeding the next:

**Molecules → DNA structure → genomes → populations → deep time**

Every module ships only after its benchmark suite reproduces published reference values within each method's stated error. Athanor does not try to replace mature free tools that already do a job well (IGV, UGENE, Benchling, ChimeraX); where a module needs a small piece of that functionality, it builds only that piece.

## Order of work

1. **Now:** Personal Genome, build step 2 (modern ancestry, appearance, food and drink, fitness)
2. **Next:** DNA Structure, first release
3. **Then:** Deep Time, starting with ancient nuclear genomes and then traits from ancient genomes, interleaved with Personal Genome steps 3–5
4. **Later:** Populations, RNA and proteins, Radiocarbon dating, and the research problems

---

## Molecules (complete; 0.1–0.2.x)

The chemistry layer. Its 3D viewer and quantum engines are the foundation of DNA Structure.

Done:
- Three chemistry tiers behind one interface, with a provenance record on every result:
  - **MMFF94** (RDKit): geometry optimization; declines molecules outside its fitted parameters
  - **GFN2-xTB** (xtb program): charges, dipoles, Hessian-based frequencies and IR spectra
  - **DFT/HF** (Psi4): energies, gradients, dipoles, frequencies; optional GPU DFT (GPU4PySCF) with analytic GPU Hessians
- Correctness checks on every build: Halgren's MMFF94 validation suite; tblite vs the xtb program; Psi4's and ASE's vibrational analyses; PySCF vs Psi4
- Accuracy benchmark against NIST reference geometries, dipoles and frequencies, gated at twice each method's published error. B3LYP passes all gates: bonds 0.0097 Å, dipoles 4.2%, scaled frequencies 0.88%
- Molecule workbench UI: build/import, optimize, IR spectrum, animated normal modes, error bars, Validation view

Remaining (only what DNA Structure needs):
- Noncovalent interaction energies (base pairing and stacking), benchmarked against the S22 set (Jurečka et al., Phys. Chem. Chem. Phys. 8, 1985, 2006), which includes adenine–thymine and guanine–cytosine pairs
- Molecular orbital and electron-density views

Parked:
- MMFF94 accuracy set for organic molecules (blocked on access to experimental structure tables)
- Conformer search

## Ancient DNA Lab (working; 0.3.0–0.4.0)

Done:
- **Simulated mode:** gargammel reads from real mitogenomes (Briggs-model damage, log-normal fragment lengths, ART Illumina errors, optional human contamination), every base scored against the known truth. CI gates on damage recovery, reference-bias reduction, accuracy (≥ 99.9%) and calibration (ECE ≤ 0.002)
- **Import:** paired-end or single-end FASTQ, or BAM; read in place, never uploaded
- **Circular mapping:** BWA aln with ancient-DNA settings to the reference elongated by its own start, positions folded back; per-library duplicate removal
- **Iterative, indel-aware reconstruction** from the species' or a relative's reference until the sequence stops changing, with per-base posteriors
- **Authentication:** terminal C→T / G→A with 95% intervals, whole-molecule length distribution, coverage, mapDamage2 profiles
- **Contamination:** competitive mapping against the reconstructed sequence and the human mitogenome (Feuerborn et al. 2020), exact binomial interval
- **Validation in CI:** contamination estimates on simulations (0%, 2%, 10%); mammoth reconstructed from the African elephant genome; reproduction of the published Hohle Fels mammoth mitogenomes (ENA PRJEB103048 → GenBank PZ244049–PZ244060)

Next:
- **Read viewer:** the reads behind any position of a reconstruction, so a disputed base can be inspected directly
- **Mitogenome tree:** where a reconstructed mitogenome sits among published relatives
- More reproduction datasets: a library without damage repair; a Neanderthal mitogenome with an established contamination estimate
- Human samples: contamination needs a different method (e.g. haplotype-based), since competitive mapping cannot separate human from human

Open questions:
- mapDamage2's δs differs from gargammel's single-stranded deamination parameter on every reference (≈ 0.55 vs 0.3): likely different parameterisations
- Residual overconfidence at single-read sites; candidate fix: a site-heterogeneous divergence prior

## Personal Genome (working; 0.5.0–0.6.0)

Load a consumer raw DNA file and get a library of reports, each graded for evidence, validated against data with known answers, and computed on the user's computer. The file never leaves the machine.

### Foundations
- **Inputs:** AncestryDNA, 23andMe, MyHeritage, FamilyTreeDNA and Living DNA raw files; VCF from whole-genome sequencing. Detect the vendor and genome build; settle strand from unambiguous markers; lift over between GRCh37 and GRCh38.
- **Imputation** (optional, Linux engines pack): Beagle 5 against the 1000 Genomes panel, with accuracy reported per allele-frequency bin from masked-genotype tests.
- **Reference data pack:** the Allen Ancient DNA Resource (AADR), 1000 Genomes and HGDP at consumer-chip positions, archaic genomes, PhyloTree and Y-chromosome trees; later CPIC and PharmVar tables, ClinVar, PGS Catalog scores. Every file has a recorded version and hash.
- **Every report card** shows the result, an evidence grade (A established, B moderate, C weak or fun), the effect size, what it cannot tell you, its sources and how it was validated.
- **Gene lookup:** any gene or variant named in a report links to its Ensembl page and, where relevant, to DNA Structure or the protein view.

### Report catalogue (about 70)

**Your file (3).** File check: vendor, chip, build, call rate. Genetic sex check. Imputation coverage.

**Deep ancestry (10), done.** Ancient ancestry proportions (qpAdm on AADR). Closest ancient individuals on a map. Neanderthal and Denisovan ancestry. Neanderthal-derived variants. mtDNA and Y haplogroups. Ancient individuals who shared them. Haplogroup timeline.

**Modern ancestry (6), done.** Continental and regional composition (supervised admixture against present-day genomes in the AADR, with panel frequencies shrunk towards the pooled frequency so panels of different sizes compete fairly). Chromosome painting (continental). Indigenous American ancestry, with a measured detection limit, stating the reference-panel limits plainly and that it is not tribal membership. Runs of homozygosity. X-chromosome ancestry. (Placement among reference populations is the PCA in Deep ancestry.)

**Relatives: two files (4).** Relationship estimate from shared DNA. Shared-segment chromosome browser. X-chromosome sharing. Parent–child check.

**Appearance (8 done, 2 with imputation).** Done: eye colour (HERC2, OCA2), red hair (MC1R), freckling (IRF4, MC1R), the KITLG blond-hair variant, skin pigmentation variants (SLC24A5, SLC45A2), earwax type (ABCC11), hair thickness (EDAR), hair form (TCHH). HIrisPlex-S is not used: its model parameters are not published (its web tool re-estimates them), its terms reserve all rights, and part of it is patented. Asparagus smell, cilantro taste and the photic sneeze reflex were dropped: the GWAS Catalog has no association for the cilantro variant, its asparagus record does not say which allele impairs smelling, and the two photic sneeze studies report opposite alleles. With imputation (step 4): male-pattern baldness and height.

**Food and drink (8 done, 2 with imputation).** Done: lactose tolerance (LCT, with the African and Arabian variants), alcohol flush (ALDH2, ADH1B), caffeine metabolism (CYP1A2, AHR), bitter taste (TAS2R38), vitamin D (GC, CYP2R1), vitamin B12 and secretor status (FUT2), omega-3 and omega-6 conversion (FADS1), celiac-compatible HLA types (mainly useful to rule out). With imputation (step 4): body-mass polygenic score, coffee and alcohol intake tendencies.

**Fitness and sleep (1 done, 2 with imputation).** Done: ACTN3. With imputation (step 4): chronotype and sleep duration.

**Medicines: pharmacogenomics (11).** CPIC level-A gene–drug pairs only: CYP2C19, CYP2C9 and VKORC1, CYP2D6 (copy-number limits stated), SLCO1B1, DPYD, TPMT and NUDT15, UGT1A1, CYP3A5, G6PD, HLA-B tag variants. Framed as something to show a clinician, never a reason to change a medicine.

**Health, opt-in and separated (10).** Blood type. Factor V Leiden and prothrombin. HFE haemochromatosis. Alpha-1 antitrypsin. Carrier status, with chip limits. APOE behind an explicit opt-in. Polygenic risk scores with each score's published accuracy and ancestry transferability. ClinVar scan in which any rare pathogenic call is marked as needing clinical confirmation (about 40% were false positives in clinical follow-up: Tandy-Connor et al., Genet. Med. 2018).

**Declined, with the reason shown.** Personality and intelligence or education scores. Genotype-based diet plans (DIETFITS, Gardner et al., JAMA 2018). "Athletic potential" scores. Famous-ancestor claims. Sexual orientation.

### Validation
- **Ancestry:** reference individuals cut down to each vendor's chip positions, scored against their known populations with leave-one-out
- **Haplogroups:** published calls for AADR and 1000 Genomes individuals (done for deep ancestry)
- **Pharmacogenomics:** GeT-RM consensus star alleles
- **Polygenic scores:** must match pgsc_calc on the same input
- **Relatives:** known pairs in 1000 Genomes, plus simulated pedigrees
- **Imputation:** masked-genotype concordance
- **Rare-variant gate:** CI fails if any rare variant is shown without the confirmation warning

### Build order
1. **Done (0.5.0):** foundations, file check, deep ancestry, haplogroups, Neanderthal and Denisovan
2. **Done (0.6.0):** modern ancestry, appearance, food and drink, fitness. Every trait variant's effect allele is checked in CI against the amino acid it codes (Ensembl VEP), its population frequencies (1000 Genomes) or the direction recorded in the GWAS Catalog
3. Pharmacogenomics
4. Relatives, imputation, and the polygenic reports that need it (height, baldness, body mass, chronotype, sleep duration, coffee and alcohol intake). A raw chip file covers well under half of a modern score's variants, and the PGS Catalog's own calculator (pgsc_calc) by default refuses to score below 75% overlap, so these wait for imputation
5. Health reports (opt-in)

### Open questions
- Regional ancestry is limited by reference panels of a few dozen people per region; finer regions need larger panels
- Several world regions have too few reference genomes to report (listed on the report cards)

## DNA Structure (new)

The double helix as a real molecule, built from measured parameters rather than drawn as a generic spiral. Uses the Molecules layer's 3D viewer and quantum engines.

- **Any sequence in 3D, atom by atom,** using sequence-dependent base-pair step parameters (Olson et al., PNAS 95, 11163, 1998), so each sequence bends and twists as real DNA does; B-, A- and Z-forms
- **Groove shape:** minor-groove width, roll, propeller twist and helix twist predicted from sequence (DNAshape; Zhou et al., Nucleic Acids Res. 41, W56, 2013), shown as colour on the helix
- **Melting:** a temperature slider unzips the helix, AT-rich regions first, from nearest-neighbour thermodynamics (SantaLucia, PNAS 95, 1460, 1998)
- **Why the helix holds together:** base-pairing and stacking energies computed with the Molecules engines
- **Ancient damage in place:** cytosine deamination at the single-stranded ends of an ancient molecule, linked to the damage the Ancient DNA Lab measures
- **Your variants in 3D:** open a Personal Genome variant and see the base change in the helix
- **Zoom through the scales:** helix → nucleosome (147 bp; PDB 1KX5) → protein–DNA complexes → chromatin fibre → whole-chromosome 3D models from public Hi-C data (4D Nucleome)
- **Later:** DNA dynamics on the GPU (OpenMM with a current nucleic-acid force field)

Benchmarks:
- Built structures vs high-resolution crystal structures (e.g. the Dickerson–Drew dodecamer), by atomic RMSD
- Groove widths vs DNAshape's published values and crystal structures
- Melting temperatures vs the published measured sets
- Pairing and stacking energies vs the S22 reference values
- Dynamics (later): helical parameters vs published simulations and NMR data

## Deep Time (new)

Ancient biology beyond the mitochondrial genome, building on the Ancient DNA Lab.

1. **Ancient nuclear genomes:** native whole-genome mapping (minimap2, bwa-mem2, samtools in the Linux engines pack); genetic sex, heterozygosity, population-size history, divergence and interbreeding between species. *Benchmark:* reproduce the published values for a published ancient genome.
2. **Traits from ancient genomes:** evidence-graded trait reports for ancient individuals (e.g. pigmentation in ancient humans; coat colour and cold-adaptation variants in mammoths), reusing the Personal Genome report cards. *Benchmark:* reproduce published trait calls for published individuals.
3. **Ancient proteins:** placing species in the tree of life from proteins preserved in tooth enamel, well beyond the age limit of DNA. Scope (starting from published peptide sequences or from mass-spectrometry data) to be settled. *Benchmark:* reproduce a published phylogenetic placement.
4. **DNA decay:** expected DNA fragment lengths and survival from a sample's age and burial temperature (Allentoft et al., Proc. R. Soc. B 279, 4724, 2012), used as an authentication check against what a sample actually yields. *Benchmark:* published fragment-length data for dated samples.
5. **Radiocarbon dating:** convert a radiocarbon measurement to a calendar date range with the IntCal20 calibration curve (Reimer et al., Radiocarbon 62, 725, 2020), linked to the dates on every ancient sample. *Benchmark:* match OxCal's calibrated ranges on a set of test dates.

## Populations (new)

An evolution sandbox: simulate populations (msprime, already in the toolchain; SLiM for selection) and watch drift, selection, bottlenecks and interbreeding play out, then compare the results with real ancient data, such as a population shrinking toward extinction.

*Benchmarks:* analytic expectations (loss of heterozygosity, fixation probabilities, coalescent times) and published simulation results.

## RNA and proteins (new)

- **RNA folding:** predict how an RNA strand folds onto itself, with base-pair probabilities. *Benchmark:* published accuracy on experimentally solved structures.
- **Protein variant view:** AlphaFold and Protein Data Bank structures (Mol* viewer), with a Personal Genome variant placed on the residue it changes and its surroundings shown. *Benchmark:* residue mapping checked against UniProt and Ensembl.

## Research problems

- **Chip false-positive detection** in consumer DNA reports: flag rare-variant calls likely to be genotyping errors, validated against sequencing data. Very rare pathogenic calls on chips are mostly wrong (Weedon et al., BMJ 372, n214, 2021). Likely the first paper.
- **Mitonuclear compatibility** at the contact surfaces between mitochondrial and nuclear protein subunits.

## Cut, with the reason

These were planned earlier and dropped because good free tools already do them:

- **Full genome browser, live genome databases, Circos-style circular and comparative plots:** IGV, the UCSC and Ensembl browsers, Circos and CGView. Athanor keeps only the read viewer and gene lookup above.
- **Sequence editing, ORF finding, restriction maps, primer design, assembly simulation, codon optimisation, BLAST:** UGENE, SnapGene Viewer, Benchling, Primer3 and NCBI BLAST.
- **CRISPR guide design and molecular docking:** CRISPOR, AutoDock Vina with ChimeraX or PyMOL.

Parked: a small cloning-simulation module (the original "splice lab" idea), to be decided later.

## Platform

Done:
- Windows installer with signed one-click updates from a public release feed
- Linux engines pack (optional, WSL2): GPU DFT self-tested on the user's GPU (verified on an RTX 5080: GPU and CPU agree to 2×10⁻¹⁰ Eh on caffeine; 11.7× faster than CPU PySCF); analytic GPU Hessians; ancient-DNA tools (bioconda) and NCBI reference genomes with recorded accession versions and hashes

Planned:
- Native whole-genome mapping in the engines pack (for Deep Time)
- Smaller updates (download only changed parts instead of the full installer)

## Before any public release

- Choose Athanor's licence (open source is required for a JOSS paper and for scientific trust); review GPL/LGPL components and NVIDIA CUDA redistribution terms together
- The Y-chromosome tree ships via Yleaf (GPL-3); the underlying YFull tree's terms need checking, and an alternative tree may be needed
- The UCSC liftOver chain is free for academic use only; a commercial release needs a licence or a different chain
- Data and tool licences to check before shipping: SNPedia (CC BY-NC-SA; avoided if Athanor is sold), PharmGKB (CC BY-SA 4.0), CPIC, ClinVar, gnomAD, PGS Catalog (per score), Ensembl and GWAS Catalog data (used to build the trait table), Beagle and PLINK 2 (GPL-3); ADMIXTURE and PLINK 1.9 are used only in CI and not shipped, DNAshape, ViennaRNA, SLiM (GPL-3), IntCal20, 4D Nucleome data use, AlphaFold DB (CC BY 4.0)
- Code-signing certificate for the installer
