# Changelog

All notable changes to Athanor are recorded here. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow [Semantic Versioning](https://semver.org/).

## [Unreleased]

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
