# Changelog

All notable changes to Athanor are recorded here. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow [Semantic Versioning](https://semver.org/).

## [Unreleased]

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
