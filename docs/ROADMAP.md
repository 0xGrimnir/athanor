# Roadmap

Each phase ships only after its benchmark suite reproduces published reference values within each method's stated error.

## Phase 1: Desktop shell and chemistry engine (in progress)

Done:
- Tiered chemistry behind one interface, with a provenance record on every result:
  - **MMFF94** (RDKit): geometry optimization; declines molecules outside its fitted parameters
  - **GFN2-xTB** (xtb program): charges, dipoles, Hessian-based frequencies and IR spectra
  - **DFT/HF** (Psi4): energies, gradients, dipoles, frequencies
- Correctness checks on every build: Halgren's MMFF94 validation suite; tblite vs the xtb program; Psi4's and ASE's vibrational analyses
- Accuracy benchmark against NIST reference geometries, dipoles and frequencies, gated at twice each method's published error
- Molecule workbench UI: build/import, optimize, IR spectrum, animated normal modes, error bars, Validation view
- Tauri 2 Windows shell and installer pipeline with a packed science runtime
- Published installers (0.1.0, 0.2.0) with signed in-app updates from a public release feed
- DFT benchmark passes all gates (B3LYP: bonds 0.0097 Å, dipoles 4.2%, frequencies 0.88%)

Remaining:
- MMFF94 accuracy set for organic molecules (currently water only): blocked on access to experimental structure tables (CCCBDB per-molecule pages; Halgren 1996 Part IV)
- Smaller updates (only changed parts instead of the full installer)
- Molecular orbitals and electron-density views
- Conformer search

## Phase 1b: Molecular dynamics

- OpenMM on CUDA, with force-field assignment and trajectory playback
- **Benchmark:** published liquid densities and diffusion coefficients

## Phase 2: Ancient DNA Lab (working, 0.3.0)

Done: gargammel simulation from real mitogenomes, AdapterRemoval merging, BWA aln with ancient-DNA settings, mapDamage2 damage estimation, damage-aware consensus with calibrated posteriors, iterative mapping against relatives' genomes, scoring against the truth, CI gates. Remaining: contamination estimation, user-supplied FASTQ/BAM input (real ancient data), more organisms.

Original plan:

- Read simulation modeled on gargammel (Renaud et al. 2017): empirical fragment lengths, the Briggs et al. (2007) damage model, Illumina-style sequencing errors, optional contamination
- Mapping with Bowtie2 using ancient-DNA settings
- Damage-pattern estimation, validated against mapDamage2 (Jónsson et al. 2013)
- Consensus calling with per-base genotype likelihoods
- Per-site error posteriors from a substitution model, calibrated on regions with ancient coverage
- Reference-bias analysis across multiple relatives (Günther & Nettelblad 2019)
- **Benchmark:** recover the damage parameters that were simulated in

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

## Heavy-compute pack

- GPU DFT (GPU4PySCF in WSL2), installed from inside the app and self-tested on the user's GPU (verified on an RTX 5080 in 0.2.3: GPU and CPU agree to 2×10⁻¹⁰ Eh on caffeine; 11.7× faster than CPU PySCF)
- Analytic GPU Hessians and IR intensities, validated on the user's GPU before use (0.2.5)
- With the genomics modules: native whole-genome mapping (minimap2, bwa-mem2, samtools) in the same distribution

## Afterward: open research problems

- Mitonuclear compatibility prediction at the contact surfaces between mitochondrial and nuclear subunits
- Chip false-positive detection in consumer DNA reports
