"""The ancient-DNA progress model: continuous, monotone, and self-correcting."""
from pathlib import Path

from athanor_core.adna.pipeline import Progress, Step


def test_progress_is_monotone_and_eta_adapts(tmp_path: Path):
    events = []
    # plan expects 1 s per stage; the machine is ~2x slower
    plan = [("a", 1.0), ("b", 1.0), ("c", 1.0)]
    prog = Progress(plan, events.append)
    step = Step([], prog)
    for name in "abc":
        step.run(name, ["sleep", "2.1"], tmp_path)
    fr = [e["fraction"] for e in events]
    assert fr == sorted(fr), fr                      # never goes backwards
    assert all(0 <= f < 1 for f in fr)
    assert events[-1]["step"] == 3 and events[-1]["steps"] == 3
    # three 2.1 s stages; starting a process adds up to ~0.5 s each on slow Windows runners
    assert prog.idx == 3 and 6.2 < prog.act_done < 8.0, prog.act_done
    # after two slow stages the ETA for the last stage reflects the slowdown (> the 1 s expected)
    third = [e for e in events if e["stage"] == "c"]
    assert third[0]["eta_seconds"] >= 1.5, third[0]
