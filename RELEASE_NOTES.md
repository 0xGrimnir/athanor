Athanor 0.2.4: faster DFT optimizations, live progress, Stop fixed.

- **Stop works during GPU calculations** (it hung in 0.2.3).
- **DFT optimizations are about 3x faster from a fresh structure:** they now start with a quick GFN2-xTB pre-optimization (checkbox under the DFT options, on by default, recorded in the provenance), and stop on Gaussian's standard convergence criteria instead of over-converging.
- **Live progress:** the busy indicator shows the optimization step, energy and largest force against its target, or the Hessian gradient count for frequency jobs.
