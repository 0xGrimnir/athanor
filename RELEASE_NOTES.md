Athanor 0.2.5: fast GPU frequencies and IR spectra.

- **Analytic second derivatives on the GPU.** Frequencies and IR intensities on the DFT-GPU tier come from one SCF plus the coupled-perturbed equations instead of 6N separate gradient calculations, so IR spectra of larger molecules take minutes instead of hours.
- **Validated before use.** Click **GPU ready** and re-run the self-test: it now also compares the analytic Hessian and IR intensities with finite differences on your GPU. The fast path switches on only if they agree; otherwise frequencies keep using finite differences.
