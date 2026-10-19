Athanor 0.2.0: optional GPU compute.

- **GPU compute pack.** Click **Enable GPU compute** in the top bar to run DFT on your NVIDIA GPU through GPU4PySCF. Athanor enables WSL2 if needed (one Windows administrator prompt, possibly a restart), downloads the pack (about 1.8 GB, every part checked by SHA-256), and runs a self-test that compares GPU and CPU results on the same molecules. The **DFT·GPU** tier appears only if the self-test passes.
- The GPU tier's code path (PySCF) is checked against Psi4 on every build.
- Requires an NVIDIA driver with CUDA 13 support (version 580 or newer).

**First launch after updating** re-checks the science engines. Windows SmartScreen may warn about an unknown publisher because the installer is not code-signed yet: choose **More info → Run anyway**.
