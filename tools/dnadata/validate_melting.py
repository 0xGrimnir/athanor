"""Melting temperatures against measurements and an independent implementation.

1. Owczarzy et al., Biochemistry 43, 3537 (2004): 92 DNA duplexes (10-30 bp) melted at five sodium
   concentrations (69, 119, 220, 621 and 1020 mM; 2 uM total strand). The measured values are read
   from the data files of TfReg (Weber, UFMG; GPL-3.0), downloaded here and not redistributed.
2. Biopython's Bio.SeqUtils.MeltingTemp.Tm_NN with the same parameters and salt corrections must agree
   to 0.01 C on every duplex (an implementation check).

Usage: python validate_melting.py OUT.json
"""
from __future__ import annotations

import io
import json
import sys
import tarfile
import urllib.request
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "core"))
from athanor_core.dna import melting as M  # noqa: E402

SRC = "https://download.opensuse.org/repositories/home:/drgweber/Debian_12/tfreg_8.1.orig.tar.gz"
OUT = Path(sys.argv[1])
results = []


def gate(name, ok, detail):
    results.append({"name": name, "pass": bool(ok), "detail": detail})
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


data = urllib.request.urlopen(urllib.request.Request(SRC, headers={"User-Agent": "athanor-ci"}), timeout=300).read()
rows = []
with tarfile.open(fileobj=io.BytesIO(data)) as t:
    m = next(x for x in t.getmembers() if x.name.endswith("data/owczarzy04.dat"))
    for line in t.extractfile(m).read().decode().splitlines()[1:]:
        f = line.split()
        if len(f) >= 5:
            rows.append((f[0], f[1], float(f[2]), float(f[3]) / 1000, float(f[4]) * 1e-6))
seqs = sorted({r[0] for r in rows})
print(f"{len(rows)} measurements of {len(seqs)} duplexes at Na+ {sorted({r[3] for r in rows})} M")
bad = [r for r in rows if M.revcomp(r[0]) != r[1][::-1]]
print(f"{len(bad)} rows whose second strand is not the complement")

table = {}
for salt in ("owczarzy2004", "santalucia1998"):
    err = {}
    for s, _, t_obs, na, ct in rows:
        e = M.tm(s, ct, na, salt)["tm"] - t_obs
        err.setdefault(na, []).append(e)
    table[salt] = {f"{na * 1000:.0f} mM": {"mean_abs": float(np.mean(np.abs(v))), "mean": float(np.mean(v)),
                                         "max_abs": float(np.max(np.abs(v))), "n": len(v)} for na, v in sorted(err.items())}
    allv = np.concatenate(list(err.values()))
    table[salt]["all"] = {"mean_abs": float(np.mean(np.abs(allv))), "mean": float(np.mean(allv)), "n": len(allv)}
    print(salt, json.dumps(table[salt], indent=1))
o = table["owczarzy2004"]
gate("melting temperatures vs 92 measured duplexes at 1 M Na+ (Owczarzy 2004): mean error <= 2.0 C",
     o["1020 mM"]["mean_abs"] <= 2.0, f"mean |error| {o['1020 mM']['mean_abs']:.2f} C, mean {o['1020 mM']['mean']:+.2f} C, "
     f"worst {o['1020 mM']['max_abs']:.1f} C")
gate("melting temperatures vs the same duplexes at 69-621 mM Na+ (Owczarzy 2004 salt correction): mean error <= 2.0 C",
     o["all"]["mean_abs"] <= 2.0, "; ".join(f"{k} {v['mean_abs']:.2f}" for k, v in o.items() if k != "all"))

from Bio.SeqUtils import MeltingTemp as mt  # noqa: E402
diff = []
for s in seqs:
    for na, salt, corr in ((1.0, "santalucia1998", 5), (0.1, "santalucia1998", 5), (0.1, "owczarzy2004", 6)):
        a = M.tm(s, 2e-6, na, salt)["tm"]
        b = mt.Tm_NN(s, nn_table=mt.DNA_NN3, saltcorr=corr, Na=na * 1000, K=0, Tris=0, Mg=0, dNTPs=0,
                     dnac1=1000, dnac2=1000)
        diff.append(abs(a - b))
gate("same melting temperatures as Biopython's Tm_NN (Allawi & SantaLucia parameters, both salt corrections)",
     max(diff) <= 0.01, f"largest difference {max(diff):.4f} C over {len(diff)} cases")
OUT.write_text(json.dumps({"results": results, "errors": table}, indent=1))
sys.exit(0 if all(r["pass"] for r in results) else 1)
