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

## Phase 3: Genome module and live data

- igv.js browser: BAM/CRAM, VCF, BED, bigWig, GFF; pileups, coverage, Sashimi plots
- Live Ensembl, UCSC and NCBI access from the user's machine
- Alignment-based identity (Needleman-Wunsch, MAFFT), with the identity definition and denominator stated on every value
- Circular and comparative views: synteny ribbons, GC skew, multi-ring plots
- DNA Report cross-checked against ClinVar and gnomAD
- **Benchmark:** published pairwise identity values

## Phase 4: Sequence workbench and structure

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
