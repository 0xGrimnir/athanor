Athanor 0.2.1: faster first launch.

- The one-time "Relocating engine paths" step after installing or updating now skips the engine files that cannot contain build paths, reads the rest in parallel, and shows a live file count.
- Includes everything in 0.2.0, including the optional GPU compute pack (**Enable GPU compute** in the top bar).

Updating from inside Athanor re-runs this one-time step with the faster code.
