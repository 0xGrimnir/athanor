Athanor 0.1.0, the first pre-release: the Phase 1 molecule workbench.

- Three chemistry tiers: MMFF94 (RDKit), GFN2-xTB (xtb) and DFT (Psi4), with geometry optimization, dipoles, charges, harmonic frequencies, IR spectra and animated normal modes.
- Every result carries its provenance (engine, version, method, parameters, input hash).
- Validated on every build: Halgren's MMFF94 suite; xtb against the independent tblite code; vibrational analysis against Psi4 and ASE. DFT accuracy vs NIST data (B3LYP): bond lengths 0.0097 Å MAD, dipoles 4.2% regularized RMS, scaled frequencies 0.88% RMS, each within its published benchmark.
- Signed in-app updates.

**First launch** unpacks the science engines (about a minute). Windows SmartScreen may warn about an unknown publisher because the installer is not code-signed yet: choose **More info → Run anyway**.
