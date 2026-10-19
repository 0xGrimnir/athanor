# Athanor Phase 1 benchmark

## dft  (4677.6 s)

- **Geometry** (B3LYP/cc-pVTZ): bonds n=8 MAD=0.0097 A, mean signed=+0.0039 A; angles n=1 MAD=0.06 deg
- **Dipoles** (B3LYP/aug-cc-pVTZ at B3LYP/cc-pVTZ geometry): n=20 regularized RMS=4.2%
- **Frequencies, all modes** (B3LYP/cc-pVTZ): n=69 raw MAD=59 cm-1, fitted scale=0.9662 (RMS rel 1.48%); literature scale 0.967 -> MAD=14 cm-1, RMS rel 1.46%
- **Frequencies, gas phase harmonic modes** (B3LYP/cc-pVTZ): n=40 raw MAD=68 cm-1, fitted scale=0.9665 (RMS rel 0.90%); literature scale 0.967 -> MAD=13 cm-1, RMS rel 0.88%
- [PASS] bond length MAD (A): measured 0.009674 vs published 0.0125 (limit 0.025; Janoschek2001)
- [PASS] dipole regularized RMS (%): measured 4.169 vs published 6.98 (limit 13.96; HaitHeadGordon2018)
- [PASS] scaled frequency RMS relative error (gas-phase modes): measured 0.008774 vs published 0.018 (limit 0.036; CCCBDB-R22)

