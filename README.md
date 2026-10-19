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
| **Molecule** | Build or import molecules; optimize geometry; partial charges, dipoles, harmonic frequencies, IR spectra, animated normal modes | Tiered: MMFF94 (RDKit) → GFN2-xTB (xtb) → DFT/HF (Psi4) | Phase 1, working |
| **Ancient DNA Lab** | Simulate damaged ancient reads, map them, call a consensus with per-base uncertainty, quantify reference bias | Briggs et al. (2007) damage model, Bowtie2, samtools/bcftools | Phase 2 |
| **Genome** | Genome browsing: reads, variants, coverage, multiple assemblies; circular and comparative views | igv.js, minimap2, MAFFT | Phase 3 |
| **Sequence Workbench** | Editing, ORFs, restriction maps, primer design, cloning simulation, alignments, phylogenetic trees, CRISPR guide scoring | Primer3, MAFFT/MUSCLE, FastTree/IQ-TREE, BLAST+ | Phase 4 |
| **Structure** | Protein structures, linking each genome position to its residue in 3D, docking | Mol*, AutoDock Vina | Phase 4 |
| **DNA Report** | Load a Promethease or raw genotype report locally; triage by magnitude and rarity; cross-check against ClinVar and gnomAD | In-app parser | Prototype; ClinVar/gnomAD in Phase 3 |

## Validation (Phase 1)

| Tier | Correctness check | Accuracy check against experiment |
|---|---|---|
| MMFF94 | Reproduces Halgren's MMFF94 validation suite (755 molecules to 1×10⁻⁴ kcal/mol) | Declines molecules outside MMFF94's fitted bond parameters |
| GFN2-xTB | The independent tblite library reproduces the xtb program; frequencies match xtb's own Hessian and ASE's vibrational analysis | Bond lengths, dipoles and IR fundamentals vs NIST data, compared with published GFN2-xTB statistics |
| DFT | Matches Psi4's own frequency analysis | Same reference set at B3LYP/cc-pVTZ, compared with published B3LYP statistics (Janoschek 2001; CCCBDB scale factors; Hait & Head-Gordon 2018) |

Reference values and sources: [`core/benchmarks/reference/phase1.json`](core/benchmarks/reference/phase1.json). Run locally with `python -m athanor_core.benchmark`.

## Architecture

- **Shell:** [Tauri 2](https://tauri.app). A native Windows app with a web-technology interface, shipped as a standard installer with signed in-app updates.
- **Science core:** Python package [`core/athanor_core`](core/athanor_core). RDKit, the xtb program and Psi4 run natively on Windows; every result carries a provenance record. The shell talks to it over stdin/stdout (JSON-RPC); no network port is opened.
- **Genomics engines:** WebAssembly builds from [biowasm](https://github.com/biowasm/biowasm) (samtools, bcftools, minimap2, Bowtie2, MAFFT and others) running inside the app. No Linux layer required.
- **Optional heavy-compute pack (later):** WSL2-based GPU DFT and whole-genome-scale mapping.

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
