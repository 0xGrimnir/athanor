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

Remaining:
- First published installer (runtime smoke test on Windows being fixed)
- Signed in-app updates (needs update hosting and a signing key)
- MMFF94 accuracy set for organic molecules (currently water only)
- Molecular orbitals and electron-density views
- Conformer search

## Phase 1b: Molecular dynamics

- OpenMM on CUDA, with force-field assignment and trajectory playback
- **Benchmark:** published liquid densities and diffusion coefficients

## Phase 2: Ancient DNA Lab

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

- GPU DFT (GPU4PySCF in WSL2), installed from inside the app and self-tested on the user's GPU (0.2.0)
- Next: analytic GPU Hessians and IR intensities, so frequencies of large molecules do not need 6N gradients
- With the genomics modules: native whole-genome mapping (minimap2, bwa-mem2, samtools) in the same distribution

## Afterward: open research problems

- Mitonuclear compatibility prediction at the contact surfaces between mitochondrial and nuclear subunits
- Chip false-positive detection in consumer DNA reports
