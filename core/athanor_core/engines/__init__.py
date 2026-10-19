from .base import Engine, EngineResult


def get_engine(tier: str, **kwargs) -> Engine:
    """Engines import lazily so a missing optional engine never breaks the others."""
    tier = tier.lower()
    if tier == "mmff94":
        from .mmff import MMFF94
        return MMFF94(**kwargs)
    if tier in ("gfn2-xtb", "xtb", "gfn2"):
        from .xtb import GFN2xTB
        return GFN2xTB(**kwargs)
    if tier == "dft":
        from .dft import Psi4DFT
        return Psi4DFT(**kwargs)
    raise ValueError(f"unknown tier {tier!r}")


__all__ = ["Engine", "EngineResult", "get_engine"]
