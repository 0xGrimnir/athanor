# Roadmap

Each phase ships only after its benchmark suite reproduces published reference values within each method's stated error.

## Phase 1: Molecule workbench (nearly complete; 0.1–0.2.x)

Done:
- Three chemistry tiers behind one interface, with a provenance record on every result:
  - **MMFF94** (RDKit): geometry optimization; declines molecules outside its fitted parameters
  - **GFN2-xTB** (xtb program): charges, dipoles, Hessian-based frequencies and IR spectra
  - **DFT/HF** (Psi4): energies, gradients, dipoles, frequencies
- Correctness checks on every build: Halgren's MMFF94 validation suite; tblite vs the xtb program; Psi4's and ASE's vibrational analyses; PySCF vs Psi4
- Accuracy benchmark against NIST reference geometries, dipoles and frequencies, gated at twice each method's published error. B3LYP passes all gates: bonds 0.0097 Å, dipoles 4.2%, scaled frequencies 0.88%
- Geometry optimization with Gaussian's default convergence criteria, GFN2-xTB pre-optimization for DFT, live progress
- Molecule workbench UI: build/import, optimize, IR spectrum, animated normal modes, error bars, Validation view
- Windows installer with signed one-click updates from a public release feed

Remaining:
- MMFF94 accuracy set for organic molecules (currently water only): blocked on access to experimental structure tables (CCCBDB per-molecule pages; Halgren 1996, J. Comput. Chem. 17, 587)
- Molecular orbitals and electron-density views
- Conformer search
- Smaller updates (download only changed parts instead of the full installer)

## Phase 1b: Molecular dynamics

- OpenMM on CUDA, with force-field assignment and trajectory playback
- **Benchmark:** published liquid densities and diffusion coefficients

## Phase 2: Ancient DNA Lab, simulated data (working; 0.3.0)

Done:
- gargammel simulation from real mitogenomes (Briggs-model damage, log-normal fragment lengths, ART Illumina errors, optional present-day human contamination)
- AdapterRemoval merging, BWA aln with ancient-DNA settings, mapDamage2 damage estimation
- Damage-aware consensus with per-base posteriors, checked for calibration against the known truth
- Iterative mapping against relatives' genomes to reduce reference bias
- Woolly mammoth vs African/Asian elephant demo; CI gates on damage recovery, reference-bias reduction, accuracy (≥ 99.9%) and calibration (ECE ≤ 0.002)

Open questions carried forward:
- mapDamage2's δs differs from gargammel's single-stranded deamination parameter on every reference (≈ 0.55 vs 0.3): likely different parameterisations; to be resolved from both programs' definitions
- Residual overconfidence at single-read sites (prior toward the reference); candidate fix: site-heterogeneous divergence prior

## Phase 2b: Ancient DNA Lab, real data (working; 0.4.0)

Done:
- **Import:** paired-end FASTQ (adapter trimming and merging), single-end FASTQ, or BAM (e.g. the files deposited with a publication; reads restored to sequencing orientation); read in place, never uploaded
- **Circular mapping:** BWA aln to the reference elongated by its own start (as in CircularMapper), positions folded back; per-library duplicate removal
- **Iterative, indel-aware reconstruction** from the species' or a relative's reference until the sequence stops changing; posteriors keep the starting reference as the prior
- **Damage measured against the sample's own consensus**; reads cut short by the read length carry no 3′ damage term
- **Authentication:** terminal C→T / G→A with 95% intervals against the molecule interior, whole-molecule length distribution, coverage, mapDamage2 profiles
- **Contamination:** competitive mapping against the reconstructed sequence and the human mitogenome (Feuerborn et al. 2020), exact binomial interval, length and damage of the human-assigned reads
- **Outputs:** consensus FASTA with an N mask, majority-rule consensus, per-base posteriors and depths, full JSON report with provenance
- **Validation in CI:** contamination estimates on simulations (0%, 2%, 10%); mammoth reconstructed starting from the African elephant genome; reproduction of the published Hohle Fels mammoth mitogenomes (ENA PRJEB103048 → GenBank PZ244049–PZ244060)

Next:
- More reproduction datasets (a non-UDG library; a Neanderthal mitogenome with an established contamination estimate)
- Human samples: contamination needs a different method (e.g. haplotype-based), since competitive mapping against human cannot separate human from human
- Scope: mitochondrial-scale data; ancient nuclear genomes later, with native whole-genome mapping

## Phase 3: Personal Genome (first release working; 0.5.0)

Load a consumer raw DNA file and get a library of reports, the way third-party upload sites offer dozens of them, but with every report graded for evidence, validated against data with known answers, and computed on the user's computer. The file never leaves the machine.

### Foundations
- **Inputs:** AncestryDNA, 23andMe, MyHeritage, FamilyTreeDNA and Living DNA raw files; VCF from whole-genome sequencing. Detect the vendor, chip version and genome build; normalise strands; lift over between GRCh37 and GRCh38.
- **Imputation** (optional, Linux engines pack): Beagle 5 against the 1000 Genomes panel fills in untyped variants for polygenic scores and haplotype methods. Accuracy is reported per allele-frequency bin from masked-genotype tests.
- **Reference data pack:** 1000 Genomes and HGDP genotypes at consumer-chip positions, the Allen Ancient DNA Resource (AADR), Neanderthal and Denisovan genomes, PhyloTree/ISOGG haplogroup trees, CPIC and PharmVar tables, ClinVar, PGS Catalog scores. Every file has a recorded version and hash.
- **Every report card** shows the result, an evidence grade (A established, B moderate, C weak or fun), the effect size, what it cannot tell you, its sources and how it was validated. Reports export as a document.

### Report catalogue (about 70)

**Your file (3).** File check: vendor, chip, build, call rate. Genetic sex check. Imputation coverage.

**Deep ancestry: Ice Age to Bronze Age (10)** — built first; it connects to the Ancient DNA Lab.
- Ancient ancestry proportions: hunter-gatherer, early farmer, steppe herder and other regional ancient sources (published qpAdm-style models on AADR)
- Your closest ancient individuals, on a map with dates
- Neanderthal ancestry, as a percentile among people of similar ancestry, not a falsely precise percentage
- Denisovan ancestry
- Your Neanderthal-derived variants and what each is known to be associated with
- mtDNA haplogroup and Y haplogroup (if male)
- Ancient individuals who carried your mtDNA or Y haplogroup (AADR), with places and dates
- Haplogroup timeline: when each branch arose, with uncertainty

**Modern ancestry (7).** Continental and regional composition (supervised admixture against 1000 Genomes and HGDP). Chromosome painting (local ancestry). Placement among reference populations (PCA). Indigenous American ancestry, which states the reference-panel limits plainly and that it is not tribal membership. Runs of homozygosity. Ancestry on the X chromosome.

**Relatives: two files (4).** Relationship estimate from shared DNA (IBD segments and kinship). Shared-segment chromosome browser. X-chromosome sharing. Parent–child check.

**Appearance (12).** Eye, hair and skin colour (HIrisPlex-S, a validated forensic model with published accuracy). Freckling. Red hair (MC1R). Earwax type and underarm odour (ABCC11). Hair thickness and curl (EDAR, TCHH). Male-pattern baldness (polygenic). Height (polygenic, with variance explained by ancestry). Unusual-trait corner: photic sneeze, asparagus smell detection, cilantro taste (grade C, labelled as such).

**Food and drink (11).** Lactose tolerance (LCT and the known African and Middle Eastern variants). Alcohol flush (ALDH2, ADH1B). Caffeine metabolism (CYP1A2, AHR). Bitter taste (TAS2R38). Vitamin D, vitamin B12 and omega-3 conversion (GC/CYP2R1, FUT2, FADS1/2: effect on blood levels, not a diet). Celiac-compatible HLA types (mainly useful to rule out). Body-mass polygenic score, showing how little it determines. Coffee and alcohol intake tendencies (GWAS).

**Fitness and sleep (3).** ACTN3 sprint/power variant (real but small effect). Chronotype: morning or evening person (polygenic). Sleep duration (polygenic).

**Medicines: pharmacogenomics (11).** CPIC level-A gene–drug pairs only. CYP2C19 (clopidogrel, SSRIs, PPIs); CYP2C9 and VKORC1 (warfarin); CYP2D6 (codeine, tamoxifen; copy-number limits of chips stated); SLCO1B1 (statins); DPYD (fluoropyrimidines); TPMT and NUDT15 (thiopurines); UGT1A1 (irinotecan; Gilbert syndrome); CYP3A5 (tacrolimus); G6PD; HLA-B*57:01 and HLA-B*58:01 tag variants. Framed as something to show a clinician, never a reason to change a medicine.

**Health, opt-in and separated (10).** Blood type (ABO, with stated accuracy; Rh limited on chips). Factor V Leiden and prothrombin G20210A. HFE haemochromatosis. Alpha-1 antitrypsin (Z and S). Carrier status for common recessive conditions, with chip limits. APOE, behind an explicit opt-in that explains what it does and doesn't mean. Polygenic risk scores (coronary disease, type 2 diabetes, breast and prostate cancer) with each score's published accuracy and how well it transfers across ancestries. ClinVar scan in which any rare pathogenic call is marked as needing clinical confirmation: about 40% of such calls in consumer raw data were false positives in a clinical follow-up (Tandy-Connor et al., Genet. Med. 2018).

**Declined, with the reason shown.** Personality and intelligence or education scores (little predictive power, high potential for misreading). "Ideal diet type" and genotype-based diet plans (the DIETFITS trial, Gardner et al., JAMA 2018, found that genotype did not predict which diet worked). "Athletic potential" scores. Royal, Viking or famous-ancestor claims. Sexual orientation.

### Validation
- **Ancestry:** 1000 Genomes and HGDP individuals cut down to each vendor's chip positions, scored against their known populations with leave-one-out. Reproduce published ancient-ancestry proportions for modern populations.
- **Haplogroups:** 1000 Genomes individuals with published Y haplogroups (Poznik et al. 2016) and mtDNA haplogroups.
- **Pharmacogenomics:** GeT-RM consensus star alleles for reference samples (CDC Genetic Testing Reference Materials).
- **Polygenic scores:** results must match the PGS Catalog's reference calculator (pgsc_calc) on the same input; accuracy is quoted from the score's published evaluation.
- **Relatives:** known relative pairs in 1000 Genomes, plus simulated pedigrees.
- **Imputation:** masked-genotype concordance.
- **Rare-variant gate:** CI fails if any rare variant is shown as a finding without the confirmation warning.

### Build order
1. **Done (0.5.0):** foundations, file check, deep ancestry, haplogroups, Neanderthal and Denisovan
2. Modern ancestry, appearance, food and drink, fitness and sleep
3. Pharmacogenomics
4. Relatives and imputation
5. Health reports (opt-in)

### Open questions carried forward
- The Y-chromosome tree ships via Yleaf (GPL-3); the underlying YFull tree's own terms need checking before any commercial release, and an alternative tree may be needed
- The UCSC liftOver chain is free for academic use only; a commercial release needs either a licence or a different chain
- Imputation (Beagle) is not in yet, so polygenic scores and haplotype methods wait for it

### Licences to check before shipping data
SNPedia (CC BY-NC-SA: non-commercial, so avoided if Athanor is sold), PharmGKB (CC BY-SA 4.0), CPIC, ClinVar (public domain), gnomAD, PGS Catalog (licence per score), AADR, HIrisPlex-S model terms, Beagle and PLINK 2 (GPL-3), ADMIXTURE.

## Phase 4: Genome module and live data

- igv.js browser: BAM/CRAM, VCF, BED, bigWig, GFF; pileups, coverage, Sashimi plots
- Live Ensembl, UCSC and NCBI access from the user's machine
- Alignment-based identity (Needleman-Wunsch, MAFFT), with the identity definition and denominator stated on every value
- Circular and comparative views: synteny ribbons, GC skew, multi-ring plots
- **Benchmark:** published pairwise identity values

## Phase 5: Sequence workbench and structure

- Editing, ORF finding, restriction maps, Primer3 primer design, Gibson and Golden Gate assembly simulation
- MAFFT/MUSCLE alignment, FastTree/IQ-TREE phylogenetics, BLAST+
- Codon optimization from real codon-usage tables
- CRISPR guide design with published on-target and off-target scores
- Mol* structures (PDB, AlphaFold), linking each genome position to its residue in 3D; AutoDock Vina docking

## Linux engines pack (optional, WSL2)

Done:
- GPU DFT (GPU4PySCF), self-tested on the user's GPU before use; verified on an RTX 5080 (GPU and CPU agree to 2×10⁻¹⁰ Eh on caffeine; 11.7× faster than CPU PySCF)
- Analytic GPU Hessians and IR intensities, validated on the user's GPU before use (0.2.5)
- Ancient-DNA tools (bioconda) and NCBI reference genomes with recorded accession versions and hashes (pack 2); no GPU required

Planned:
- Native whole-genome mapping (minimap2, bwa-mem2, samtools) for the genomics modules and ancient nuclear genomes

## Afterward: open research problems

- Mitonuclear compatibility prediction at the contact surfaces between mitochondrial and nuclear subunits
- Chip false-positive detection in consumer DNA reports

## Before any public release

- Choose Athanor's license (open source is required for a JOSS paper); review GPL/LGPL components and NVIDIA CUDA redistribution terms together
- Code-signing certificate for the installer
