# Roadmap

Each phase ships only after its benchmark suite reproduces published reference values within each method's stated error.

## Phase 1: Desktop shell and chemistry engine (in progress)

- Tauri 2 shell, Windows installer, signed in-app updates, GitHub Actions release pipeline
- Tiered chemistry:
  - **MMFF94** (RDKit): geometry optimization, conformer search
  - **GFN2-xTB**: partial charges, dipoles, Hessian-based vibrational frequencies and IR spectra
  - **DFT/HF** (Psi4): molecular orbitals, electron density, accurate energies
- Molecular dynamics via OpenMM on CUDA
- Provenance record attached to every result
- **Benchmark:** experimental geometries, dipoles and vibrational frequencies from the NIST CCCBDB reference set

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

## Later: optional heavy-compute pack

- WSL2-based GPU DFT and whole-genome-scale mapping

## Afterward: open research problems

- Mitonuclear compatibility prediction at the contact surfaces between mitochondrial and nuclear subunits
- Chip false-positive detection in consumer DNA reports
