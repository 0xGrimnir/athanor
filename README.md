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
| **Ancient DNA Lab** | Simulate damaged ancient reads from a real mitogenome, map them to relatives' genomes, estimate damage, call a damage-aware consensus with calibrated per-base posteriors, quantify and reduce reference bias by iterative mapping, and score every base against the known truth | gargammel, AdapterRemoval, BWA aln, mapDamage2, minimap2 (Linux engines pack) | Phase 2, working |
| **Genome** | Genome browsing: reads, variants, coverage, multiple assemblies; circular and comparative views | igv.js, minimap2, MAFFT | Phase 3 |
| **Sequence Workbench** | Editing, ORFs, restriction maps, primer design, cloning simulation, alignments, phylogenetic trees, CRISPR guide scoring | Primer3, MAFFT/MUSCLE, FastTree/IQ-TREE, BLAST+ | Phase 4 |
| **Structure** | Protein structures, linking each genome position to its residue in 3D, docking | Mol*, AutoDock Vina | Phase 4 |
| **DNA Report** | Load a Promethease or raw genotype report locally; triage by magnitude and rarity; cross-check against ClinVar and gnomAD | In-app parser | Prototype; ClinVar/gnomAD in Phase 3 |

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

## Citing

If Athanor contributes to your work, please cite it using [`CITATION.cff`](CITATION.cff), and cite the underlying engines listed in [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md).

## License

Copyright © 2026 Tyler Lee Jones. All rights reserved for now; see [`LICENSE`](LICENSE). Third-party components keep their own licenses.
