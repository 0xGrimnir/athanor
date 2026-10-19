"""CPU-only checks of the compute pack image (CI has no GPU).

GPU4PySCF initialises CUDA when imported, so without a GPU its import must
fail with a CUDA driver error (which also proves the package and its CUDA
runtime libraries load). The GPU tier's code path then runs on CPU (PySCF)
for water, and GPU detection must report the missing GPU cleanly.
"""
import json
import sys
from importlib.metadata import version

import pyscf

from athanor_core.engines import get_engine
from athanor_core.engines.pyscf_dft import gpu_available
from athanor_core.molecule import Molecule

try:
    import gpu4pyscf  # noqa: F401
    raise SystemExit("gpu4pyscf imported without a GPU; unexpected on CI")
except Exception as e:  # noqa: BLE001
    import_error = f"{type(e).__name__}: {e}"
    assert "CUDA" in import_error or "cuda" in import_error, import_error

ok, why = gpu_available()
assert not ok, "CI unexpectedly reports a usable GPU"
r = get_engine("dft-pyscf", basis="def2-SVP", density_fitting=True).compute(Molecule.from_smiles("O"))
assert -76.4 < r.energy < -76.3, r.energy
print(json.dumps({"pyscf": pyscf.__version__, "gpu4pyscf": version("gpu4pyscf-cuda13x"),
                  "cupy": version("cupy-cuda13x"), "gpu4pyscf_import_without_gpu": import_error[:160],
                  "gpu_detection": why[:160], "water_b3lyp_def2svp_ri_cpu": r.energy}))
sys.exit(0)
