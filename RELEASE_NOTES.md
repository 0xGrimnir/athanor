Athanor 0.2.3: GPU self-test fix.

- Fixes the GPU self-test error "TypeError: 'NotImplementedType' object is not callable". GPU4PySCF does not implement Mulliken charge analysis on the GPU, so Athanor now computes charges from the GPU density with PySCF's CPU routine, identically on both paths.
- **GPU: retry** re-runs the self-test without downloading the compute pack again.
