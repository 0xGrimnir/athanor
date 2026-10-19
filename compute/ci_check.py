"""CPU-only checks of the compute pack image (CI has no GPU).

Imports every engine package, confirms GPU detection fails cleanly without a
GPU, and runs the GPU tier's code path on CPU (PySCF) for water.
"""
import json
import sys

import cupy  # noqa: F401  (imports without a GPU; fails only on use)
import gpu4pyscf
import pyscf

from athanor_core.engines import get_engine
from athanor_core.engines.pyscf_dft import gpu_available
from athanor_core.molecule import Molecule

ok, why = gpu_available()
assert not ok, "CI unexpectedly reports a usable GPU"
r = get_engine("dft-pyscf", basis="def2-SVP", density_fitting=True).compute(Molecule.from_smiles("O"))
assert -76.4 < r.energy < -76.3, r.energy
print(json.dumps({"pyscf": pyscf.__version__, "gpu4pyscf": gpu4pyscf.__version__,
                  "gpu_detection_without_gpu": why, "water_b3lyp_def2svp_ri": r.energy}))
sys.exit(0)
