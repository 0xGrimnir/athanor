"""Real sequencing reads in, one FASTQ of molecules out.

Accepted inputs (processed locally; nothing is uploaded):

  * paired-end FASTQ: adapters trimmed and overlapping pairs merged with
    AdapterRemoval v3, as for the simulated data. Only merged pairs are kept:
    each is a whole molecule, so both of its ends carry the deamination
    signal (standard practice for ancient DNA; Schubert et al. 2016).
  * single-end FASTQ: adapters trimmed with AdapterRemoval v3 unless the
    reads are marked as already trimmed.
  * BAM (for example the files deposited with a publication): every primary
    read is turned back into its sequenced orientation (reads on the reverse
    strand are reverse-complemented). Paired reads are merged as above.

A single-end read that ends before the sequencing read length had its
adapter removed, so its 3' end is the molecule's end. A read that runs to
the full read length may be a longer molecule cut off by the sequencer; its
3' end carries no damage information. Such reads are marked with the ``:p``
name suffix (see pileup.py) so the consensus model and the damage estimate
treat their 3' ends accordingly.
"""
from __future__ import annotations

import gzip
import hashlib
from pathlib import Path

import numpy as np

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _fastq_records(path: Path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            s, _, q = fh.readline().strip(), fh.readline(), fh.readline().strip()
            yield h[1:].split()[0], s, q


def sha256(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def bam_to_fastq(bam: Path, out_dir: Path) -> dict:
    """Primary reads of a BAM in sequencing orientation. Returns paths of the
    single-end reads and (if any) the two mates, plus whether the BAM held
    only mapped reads (then it was filtered by mapping before deposit)."""
    import pysam
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {"single": out_dir / "single.fq", "r1": out_dir / "r1.fq", "r2": out_dir / "r2.fq"}
    fh = {k: open(v, "w") for k, v in paths.items()}
    n = {"single": 0, "r1": 0, "r2": 0, "unmapped": 0, "total": 0}
    mates: dict = {}
    with pysam.AlignmentFile(str(bam), "rb", check_sq=False) as b:
        for r in b.fetch(until_eof=True):
            if r.is_secondary or r.is_supplementary or r.query_sequence is None:
                continue
            n["total"] += 1
            n["unmapped"] += r.is_unmapped
            seq = r.query_sequence
            qual = r.query_qualities
            q = "".join(chr(x + 33) for x in qual) if qual is not None else "I" * len(seq)
            if r.is_reverse:
                seq, q = seq.translate(COMP)[::-1], q[::-1]
            if r.is_paired:
                # keep a pair only once both mates are seen
                key = r.query_name
                mate = 1 if r.is_read1 else 2
                mates.setdefault(key, {})[mate] = (seq, q)
                if len(mates[key]) == 2:
                    for m in (1, 2):
                        s, qq = mates[key][m]
                        fh[f"r{m}"].write(f"@{key}\n{s}\n+\n{qq}\n")
                        n[f"r{m}"] += 1
                    del mates[key]
            else:
                fh["single"].write(f"@{r.query_name}\n{seq}\n+\n{q}\n")
                n["single"] += 1
    for f in fh.values():
        f.close()
    return {"paths": {k: v for k, v in paths.items() if n[k] > 0}, "counts": n,
            "all_reads_mapped": n["total"] > 0 and n["unmapped"] == 0}


def prepare_reads(step, work: Path, inputs: list[dict], min_length: int = 30, threads: int = 4) -> dict:
    """Turn each input into molecules and write them all to one FASTQ with
    names ``<input>.<n>:f`` (whole molecule) or ``<input>.<n>:p`` (3' end not
    observed). Returns the FASTQ path and per-input statistics."""
    from .pipeline import BIO

    rdir = work / "reads"
    rdir.mkdir(parents=True, exist_ok=True)
    out_path = rdir / "molecules.fq.gz"
    stats, lengths_full, lengths_partial = [], [], []
    with gzip.open(out_path, "wt", compresslevel=4) as out:
        for i, inp in enumerate(inputs):
            kind, path = inp["kind"], Path(inp["path"]).resolve()
            idir = rdir / f"in{i}"
            idir.mkdir(exist_ok=True)
            st = {"index": i, "kind": kind, "file": path.name, "sha256": sha256(path)}
            singles: list[tuple[Path, bool]] = []    # (fastq, needs adapter trimming)
            pairs = None
            if kind == "bam":
                conv = step.python(f"read BAM input {i + 1}", lambda: bam_to_fastq(path, idir))
                st["bam_reads"] = conv["counts"]
                st["all_reads_mapped"] = conv["all_reads_mapped"]
                if "single" in conv["paths"]:
                    singles.append((conv["paths"]["single"], not inp.get("trimmed", True)))
                if "r1" in conv["paths"]:
                    pairs = (conv["paths"]["r1"], conv["paths"]["r2"])
            elif kind == "fastq":
                if inp.get("path2"):
                    pairs = (path, Path(inp["path2"]).resolve())
                    st["file2"] = Path(inp["path2"]).name
                    st["sha256_2"] = sha256(inp["path2"])
                else:
                    singles.append((path, not inp.get("trimmed", False)))
            else:
                raise ValueError(f"unknown input kind {kind!r}")

            sources: list[tuple[Path, bool]] = []   # (fastq, whole molecules?)
            if pairs:
                step.run(f"trim and merge input {i + 1} (AdapterRemoval)", [
                    BIO / "bin/adapterremoval3", "--in-file1", pairs[0], "--in-file2", pairs[1],
                    "--out-prefix", idir / "ar", "--merge", "--threads", threads], idir)
                sources.append((idir / "ar.merged.fastq.gz", True))
            for fq, trim in singles:
                if trim:
                    step.run(f"trim adapters input {i + 1} (AdapterRemoval)", [
                        BIO / "bin/adapterremoval3", "--in-file1", fq, "--out-file1", idir / "se.trimmed.fastq.gz",
                        "--threads", threads], idir)
                    fq = idir / "se.trimmed.fastq.gz"
                sources.append((fq, False))

            n_in = n_kept = n_full = n_partial = 0
            for fq, whole in sources:
                recs = list(_fastq_records(fq))
                read_len = int(inp.get("read_length") or max((len(s) for _, s, _ in recs), default=0))
                st.setdefault("max_read_length", read_len)
                st["max_read_length"] = max(st["max_read_length"], read_len)
                for _, s, q in recs:
                    n_in += 1
                    if len(s) < min_length:
                        continue
                    partial = (not whole) and len(s) >= read_len
                    n_kept += 1
                    if partial:
                        n_partial += 1
                        lengths_partial.append(len(s))
                    else:
                        n_full += 1
                        lengths_full.append(len(s))
                    out.write(f"@{i}.{n_kept}{':p' if partial else ':f'}\n{s}\n+\n{q}\n")
            st.update({"reads_in": n_in, "reads_kept": n_kept, "whole_molecules": n_full,
                       "cut_by_read_length": n_partial, "min_length": min_length})
            stats.append(st)
    lf = np.array(lengths_full)
    hist = np.bincount(lf, minlength=1).tolist() if lf.size else []
    return {"fastq": out_path, "inputs": stats,
            "fragment_lengths": {
                "whole_molecules": int(lf.size),
                "median": float(np.median(lf)) if lf.size else None,
                "mean": float(lf.mean()) if lf.size else None,
                "histogram": hist,
                "cut_by_read_length": len(lengths_partial),
                "note": ("Lengths of whole molecules only. Reads that ran to the full sequencing length are "
                         "longer molecules of unknown length and are counted separately.")}}
