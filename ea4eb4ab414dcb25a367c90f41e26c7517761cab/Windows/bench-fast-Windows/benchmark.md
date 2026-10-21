# Athanor Phase 1 benchmark

## mmff94  (0.0 s)

- **Geometry** (MMFF94): bonds n=1 MAD=0.0110 A, mean signed=+0.0110 A; angles n=1 MAD=0.50 deg
  - not supported by this tier: H2, HF, HCl, N2, CO, F2, Cl2

## gfn2-xtb  (24.3 s)

- **Geometry** (GFN2-xTB): bonds n=8 MAD=0.0116 A, mean signed=+0.0093 A; angles n=1 MAD=2.77 deg
- **Dipoles** (GFN2-xTB at GFN2-xTB geometry): n=20 regularized RMS=55.8%
- **Frequencies, all modes** (GFN2-xTB): n=69 raw MAD=47 cm-1, fitted scale=0.9923 (RMS rel 4.65%)
- **Frequencies, gas phase harmonic modes** (GFN2-xTB): n=40 raw MAD=56 cm-1, fitted scale=0.9914 (RMS rel 3.52%)
- [PASS] bond length MAD (A): measured 0.01164 vs published 0.0117 (limit 0.0234; Kouam2025)
- [PASS] dipole regularized RMS (%): measured 55.81 vs published 32.09 (limit 64.18; CompositeIR2024)
