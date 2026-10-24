"""Medicines: pharmacogenomics from a raw DNA file.

What is read, and how
---------------------
Allele definitions, allele functions, the function-to-phenotype tables and the
prescribing recommendations all come from CPIC (the Clinical Pharmacogenetics
Implementation Consortium; CC0), as CI downloads them into the reference pack
(tools/pgdata/build_pgx.py). Only gene-drug pairs CPIC rates level A are used.

Star alleles are called the way CPIC's tables are meant to be read and the
way PharmCAT's Named Allele Matcher reads them, so that CI can compare the two
on the same genomes (tools/pgdata/validate_pgx.py):

* a blank cell in a definition means the reference base;
* genotypes are unphased, so every pair of alleles that could produce the
  genotype at each tested position is a candidate diplotype;
* a candidate's score is the number of tested positions its two alleles
  define; the best-scoring candidates are reported, and when they disagree on
  the phenotype the result says so instead of picking one;
* DPYD and RYR1 alleles are single variants, so (as CPIC specifies) the two
  lowest-function variants found make the diplotype.

How far a call can be trusted
-----------------------------
A consumer chip reads only some of the positions in a definition. The pack
carries the 1000 Genomes high-coverage phased panel at every defined position;
for each gene the app calls all 2,504 unrelated 1000 Genomes people twice,
once with every position and once with only the positions your file has, and
reports how often the reduced call gives the same phenotype as the full one
among people who got the same reduced call as you.

Everything here runs on the user's computer.
"""
from __future__ import annotations

from collections import Counter

import numpy as np

from .rawfile import BASES, MISSING, RawGenotypes

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
NOVEL = 15            # an allele index no definition allows
UNTESTED = 255
# genes where an unlisted base at a defined position is read as the reference
# (PharmCAT does the same: these are toxicity genes, where a "not called"
# would hide an actionable variant found elsewhere in the gene)
TO_REFERENCE = {"CACNA1S", "DPYD", "G6PD", "NUDT15", "RYR1", "TPMT"}
VARIANT_GENES = {"DPYD", "RYR1"}
# unphased NAT2 genotypes often fit several diplotypes; as PharmCAT does (with UK Biobank frequencies), the
# commonest is chosen - here by star-allele frequencies on the phased 1000 Genomes haplotypes
BY_FREQUENCY = {"NAT2"}
SEVERITY = {"Malignant Hyperthermia associated": 0, "No function": 0, "Uncertain function": 1,
            "Decreased function": 1, "Normal function": 2}


def _num(v):
    try:
        return f"{float(v):g}"
    except (TypeError, ValueError):
        return None if v is None else str(v).strip().lower()


def _key(c: Counter) -> tuple:
    return tuple(sorted((_num(k), n) for k, n in c.items()))


class Gene:
    """One gene's CPIC definition table as matrices."""

    def __init__(self, symbol: str, d: dict):
        self.symbol = symbol
        self.d = d
        self.pos = d["positions"]
        self.alleles = [a for a in d["alleles"] if not a.get("structural")]
        self.names = [a["name"] for a in self.alleles]
        nA, nP = len(self.alleles), len(self.pos)
        self.allowed = np.zeros((nA, nP), dtype=np.uint16)
        self.defined = np.zeros((nA, nP), dtype=bool)
        # the reference row's cell; a set where CPIC's reference row has an IUPAC code (CYP2D6 *1)
        self.refmask = np.array([sum(1 << j for j in p.get("refset", [0])) for p in self.pos], dtype=np.uint16)
        for i, a in enumerate(self.alleles):
            for k, s in enumerate(a["def"]):
                if s is None:
                    self.allowed[i, k] = self.refmask[k]      # blank cell: the reference row's base
                else:
                    self.defined[i, k] = True
                    for j in s:
                        self.allowed[i, k] |= 1 << j
        self.ref = next((i for i, a in enumerate(self.alleles) if a.get("ref")), 0)
        # "wobble" cells: a defined cell allowing more than one base (an IUPAC code in CPIC's table). As in
        # PharmCAT, such a cell scores only when the base it is matched to is not the reference
        pop = np.vectorize(lambda x: bin(int(x)).count("1"))(self.allowed) if self.allowed.size else self.allowed
        self.wobble = self.defined & (pop > 1)
        self.wobble[self.ref] = False
        self.lookups = {}
        for r in d["results"]:
            self.lookups.setdefault(_key(Counter(r["key"])), r)
        self.by_name = {a["name"]: a for a in self.alleles}
        self.mode = "variants" if symbol in VARIANT_GENES else "diplotype"
        self.genotype_only = not d["results"]          # CPIC gives no phenotype (CYP4F2, VKORC1)
        self.freq = (d.get("hap_freq") or {}).get("all") or {}

    # ------------------------------------------------------------------ matching
    def _bits(self, geno: np.ndarray, T: np.ndarray, col: int) -> np.ndarray:
        x = geno[T, col].astype(np.int64)
        return np.left_shift(1, np.minimum(x, NOVEL))

    def consistent(self, geno: np.ndarray, haploid: bool = False):
        """(C, score, T): C[i, j] when alleles i and j (i <= j) could give the
        genotype at every tested position; score = defined tested positions."""
        T = np.flatnonzero(geno[:, 0] != UNTESTED)
        AM = self.allowed[:, T]
        score = self.defined[:, T].sum(1)
        M1 = (AM & self._bits(geno, T, 0)) != 0
        if haploid:
            return M1.all(1), score, T
        M2 = (AM & self._bits(geno, T, 1)) != 0
        P1, P2 = np.packbits(M1, axis=1), np.packbits(M2, axis=1)
        full = np.packbits(np.ones(len(T), dtype=bool))
        C = (((P1[:, None, :] & P2[None, :, :]) | (P2[:, None, :] & P1[None, :, :])) == full).all(2)
        C &= np.triu(np.ones(C.shape, dtype=bool))
        return C, score, T

    def match(self, geno: np.ndarray, haploid: bool = False):
        """geno: (positions, 2) allele indices, UNTESTED where not read.
        Returns (list of diplotypes as tuples of names, score) - all the
        best-scoring candidates - or ([], 0) when nothing fits."""
        if not (geno[:, 0] != UNTESTED).any():
            return [], 0
        C, score, _ = self.consistent(geno, haploid)
        if not C.any():
            return [], 0
        if haploid:
            best = score[C].max()
            return [(self.names[i],) for i in np.flatnonzero(C & (score == best))], int(best)
        S = score[:, None] + score[None, :] - self._wobble_penalty(geno)
        best = S[C].max()
        ij = np.argwhere(C & (S == best))
        return [tuple(sorted((self.names[i], self.names[j]), key=_order)) for i, j in ij], int(best)

    def _wobble_penalty(self, geno: np.ndarray) -> np.ndarray:
        """Per pair, the fewest wobble cells that must take a reference base,
        choosing at each heterozygous position which allele gets which base."""
        A = len(self.alleles)
        pen = np.zeros((A, A), dtype=np.int64)
        for k in np.flatnonzero(self.wobble.any(0) & (geno[:, 0] != UNTESTED)):
            x, y = (int(v) for v in np.minimum(geno[k], NOVEL))
            ax = (self.allowed[:, k] >> x) & 1 == 1
            ay = (self.allowed[:, k] >> y) & 1 == 1
            rx = bool((self.refmask[k] >> x) & 1) if x != NOVEL else False
            ry = bool((self.refmask[k] >> y) & 1) if y != NOVEL else False
            w = self.wobble[:, k].astype(np.int64)
            o1 = ax[:, None] & ay[None, :]                  # first allele takes x
            o2 = ay[:, None] & ax[None, :]
            p1 = w[:, None] * rx + w[None, :] * ry
            p2 = w[:, None] * ry + w[None, :] * rx
            big = 10 ** 6
            pen += np.minimum(np.where(o1, p1, big), np.where(o2, p2, big)) % big
        return pen

    def alternatives(self, geno: np.ndarray, haploid: bool, called: set, phenotype: str | None, limit=4):
        """Lower-scoring diplotypes that also fit the genotypes and give a
        different phenotype, built only from alleles with evidence at a tested
        position (e.g. TPMT *3B/*3C behind a *1/*3A call)."""
        C, score, T = self.consistent(geno, haploid)
        tested = np.zeros(len(self.pos), dtype=bool)
        tested[T] = True
        # only alleles every one of whose variant cells was read: the others are the "alleles your file
        # cannot see", reported separately
        nonref = self.defined & ~(self.allowed == self.refmask[None, :])
        evid = np.array([i == self.ref or (bool(nonref[i].any()) and bool(tested[nonref[i]].all()))
                         for i in range(len(self.alleles))])
        out = []
        pairs = [(i,) for i in np.flatnonzero(C)] if haploid else [tuple(x) for x in np.argwhere(C)]
        for pr in pairs:
            if not all(evid[i] for i in pr):
                continue
            names = tuple(sorted((self.names[i] for i in pr), key=_order))
            if "/".join(names) in called:
                continue
            p = self.phenotype(names)
            r = (p or {}).get("result") or "Indeterminate"
            if r != phenotype:
                out.append({"diplotype": "/".join(names), "phenotype": r})
        out.sort(key=lambda x: x["diplotype"])
        return out[:limit]

    def by_frequency(self, geno: np.ndarray):
        """The diplotype with the highest product of haplotype frequencies
        among all that fit, and how much of the total it carries."""
        C, _, _ = self.consistent(geno, False)
        ij = np.argwhere(C)
        f = np.array([max(self.freq.get(n, 0.0), 1e-5) for n in self.names])
        w = f[ij[:, 0]] * f[ij[:, 1]] * np.where(ij[:, 0] == ij[:, 1], 1.0, 2.0)
        best = int(np.argmax(w))
        names = tuple(sorted((self.names[ij[best, 0]], self.names[ij[best, 1]]), key=_order))
        by_ph = {}
        for (i, j), x in zip(ij, w):
            r = (self.phenotype((self.names[i], self.names[j])) or {}).get("result") or "Indeterminate"
            by_ph[r] = by_ph.get(r, 0.0) + x
        tot = float(w.sum())
        return [names], {"chosen_by": "haplotype frequency", "diplotype_probability": float(w[best]) / tot,
                         "phenotype_probabilities": {k: v / tot for k, v in sorted(by_ph.items(), key=lambda x: -x[1])}}

    def _is_ref_cell(self, i):
        """Positions where allele i's definition is just the reference base."""
        return self.allowed[i] == self.refmask

    def variants_found(self, geno: np.ndarray, haploid: bool = False):
        """DPYD/RYR1: each non-reference allele is a variant; copies of each
        found (0-2), with a note when a two-position allele (DPYD HapB3) was
        read from one position only."""
        found, notes = [], []
        tested = geno[:, 0] != UNTESTED
        g = np.minimum(geno.astype(np.int64), NOVEL)
        cols = g[:, :1] if haploid else g
        nonref = tested & (((self.refmask[:, None] >> cols) & 1) == 0).any(1)
        cand = np.flatnonzero((self.defined & nonref[None, :]).any(1))
        for i in cand:
            a = self.alleles[i]
            if i == self.ref:
                continue
            ks = np.flatnonzero(self.defined[i])
            ks = [k for k in ks if geno[k, 0] != UNTESTED]
            if not ks:
                continue
            prefer = [k for k in ks if self.pos[k].get("primary")]
            use = prefer or ks
            if len(np.flatnonzero(self.defined[i])) > len(use) and not prefer:
                notes.append(f"{a['name']} was read from {', '.join(self.pos[k]['rsid'] or self.pos[k]['hgvs'] for k in use)} only")
            counts = []
            for k in use:
                alleles = geno[k, :1] if haploid else geno[k]
                counts.append(int(sum(1 for x in alleles if x != NOVEL and (self.allowed[i, k] >> int(x)) & 1
                                      and not (self.refmask[k] >> int(x)) & 1)))
            c = min(counts)
            if c:
                found.append((i, c))
        # a variant that is part of a longer allele found as often (DPYD
        # c.1129-5923C>G inside HapB3) is that allele, as in PharmCAT
        keep = []
        for i, c in found:
            di = set(np.flatnonzero(self.defined[i]))
            if any(j != i and c2 >= c and di < set(np.flatnonzero(self.defined[j])) for j, c2 in found):
                continue
            keep.append((self.alleles[i]["name"], c))
        return keep, notes

    # ---------------------------------------------------------------- phenotype
    def value(self, name: str):
        a = self.by_name.get(name)
        if a is None:
            return None
        return a.get("activity") if self.d["lookup"] == "ACTIVITY_SCORE" else a.get("function")

    def phenotype(self, names) -> dict | None:
        vals = [self.value(n) for n in names]
        if any(v is None for v in vals):
            return None
        return self.lookups.get(_key(Counter(vals)))

    def variant_diplotype(self, found, haploid=False) -> tuple:
        """The two lowest-function variants (homozygous ones count twice)."""
        n = 1 if haploid else 2
        refname = self.alleles[self.ref]["name"]
        pool = []
        for name, c in found:
            pool += [name] * c
        if self.d["lookup"] == "ACTIVITY_SCORE":
            pool.sort(key=lambda x: float(self.value(x) or 1))
        else:
            pool.sort(key=lambda x: SEVERITY.get(self.value(x), 1))
        pool = pool[:n]
        return tuple(sorted(pool + [refname] * (n - len(pool)), key=_order))

    def call(self, geno: np.ndarray, haploid: bool = False, alternatives: bool = False) -> dict:
        """The full call: diplotypes, phenotype, and what fits equally well."""
        tested = int((geno[:, 0] != UNTESTED).sum())
        novel = [k for k in range(len(self.pos)) if geno[k, 0] != UNTESTED and (NOVEL in geno[k, :1 if haploid else 2])]
        notes = []
        if novel and self.symbol in TO_REFERENCE:
            geno = geno.copy()
            for k in novel:
                geno[k][geno[k] == NOVEL] = 0
            notes.append("a base no CPIC allele lists was read as the reference at "
                         + ", ".join(self.pos[k]["rsid"] or self.pos[k]["hgvs"] for k in novel))
        if tested == 0:
            return {"tested": 0, "diplotypes": [], "phenotype": None, "status": "not tested", "notes": notes}
        if self.mode == "variants":
            found, n2 = self.variants_found(geno, haploid)
            dips = [self.variant_diplotype(found, haploid)]
            notes += n2
            extra = {"variants": [{"name": n, "copies": c, "function": self.by_name[n].get("function"),
                                   "activity": self.by_name[n].get("activity")} for n, c in found]}
        else:
            dips, _ = self.match(geno, haploid)
            extra = {}
            if dips and self.symbol in BY_FREQUENCY and self.freq and not haploid:
                dips, extra = self.by_frequency(geno)
            if not dips:
                why = (" (a base no CPIC allele lists, at " + ", ".join(self.pos[k]["rsid"] or self.pos[k]["hgvs"] for k in novel) + ")") if novel else ""
                return {"tested": tested, "diplotypes": [], "phenotype": None, "status": "no call" + why, "notes": notes}
        out = {"tested": tested, "diplotypes": ["/".join(dp) for dp in dips], "notes": notes, **extra}
        if self.genotype_only:
            out.update(status="called" if len(dips) == 1 else "ambiguous", phenotype=None, genotype_only=True)
            return out
        phen = [self.phenotype(dp) for dp in dips]
        results = sorted({(p or {}).get("result") or "Indeterminate" for p in phen})
        if len(results) == 1:
            p = phen[0] or {}
            out.update(status="called", phenotype=results[0], score=p.get("score"),
                       description=p.get("description"), priority=p.get("priority"))
            if alternatives and self.mode == "diplotype":
                out["alternatives"] = self.alternatives(geno, haploid, set(out["diplotypes"]), results[0])
        else:
            out.update(status="ambiguous", phenotype=None, possible=results)
        return out


def _order(name: str):
    """Sort star alleles numerically (*2 before *10), reference first."""
    import re
    m = re.match(r"\*(\d+)", name or "")
    return (0, int(m.group(1)), name) if m else (1, 0, name or "")


# ---------------------------------------------------------------------- input
class PgxGenotyper:
    """Allele indices for CPIC positions from a raw file (GRCh37 positions)."""

    def __init__(self, g: RawGenotypes, pos37: np.ndarray, strand: str):
        self.g = g
        self.flip = {"the reference's + strand": False, "the opposite strand": True}.get(strand)
        called = g.a1 != MISSING
        self.by_pos, self.by_rs = {}, {}
        for i in np.flatnonzero(called & (pos37 > 0)):
            self.by_pos.setdefault((int(g.chrom[i]), int(pos37[i])), []).append(int(i))
        for i in np.flatnonzero(called & (g.rsid > 0)):
            self.by_rs.setdefault(int(g.rsid[i]), int(i))
        self.chips = self.ambiguous = 0

    def _at(self, chrom, pos37, indel: bool):
        """The row at a position whose call is of the right kind (a base, or I/D):
        chips list an indel and a SNV at the same coordinate."""
        for i in self.by_pos.get((chrom, pos37), []) if pos37 else []:
            if (self.g.a1[i] > 3) == indel:
                return i
        return None

    def row(self, chrom: int, p: dict):
        rs = int(p["rsid"][2:]) if (p.get("rsid") or "").startswith("rs") else None
        indel = p["kind"] != "snv"
        by_rs = self.by_rs.get(rs) if rs else None
        if by_rs is not None and (self.g.a1[by_rs] > 3) != indel:
            by_rs = None
        if not indel:
            i = self._at(chrom, p.get("pos37"), False)
            return i if i is not None else by_rs
        return by_rs if by_rs is not None else self._at(chrom, p.get("pos37"), True)

    def genotype(self, gene: Gene, chrom: int) -> np.ndarray:
        geno = np.full((len(gene.pos), 2), UNTESTED, dtype=np.uint8)
        for k, p in enumerate(gene.pos):
            i = self.row(chrom, p)
            if i is None:
                continue
            b1, b2 = int(self.g.a1[i]), int(self.g.a2[i])
            if p["kind"] == "snv":
                if b1 > 3 or b2 > 3:
                    continue
                x = [BASES[b1], BASES[b2]]
                allowed = set(p["alleles"])
                if self.flip is True:
                    x = [COMP[c] for c in x]
                elif self.flip is None:
                    if len(allowed) == 2 and COMP[p["alleles"][0]] == p["alleles"][1]:
                        self.ambiguous += 1
                        continue
                    if not set(x) <= allowed and {COMP[c] for c in x} <= allowed:
                        x = [COMP[c] for c in x]
                geno[k] = [p["alleles"].index(c) if c in allowed else NOVEL for c in x]
            else:
                chip = p.get("chip") or {}
                if b1 < 4 or b2 < 4 or not chip:
                    continue           # an indel needs an I/D call
                geno[k] = [chip.get(BASES[b], NOVEL) for b in (b1, b2)]
        return geno


# ---------------------------------------------------------------- reliability
class Panel:
    """1000 Genomes high-coverage phased genotypes at every CPIC position."""

    def __init__(self, npz, genes: dict):
        self.z = npz
        self.pop = npz["superpop"].astype(str)
        self.male = npz["male"].astype(bool)
        self.use = ~npz["related"].astype(bool)
        self.genes = genes

    def genotypes(self, gene: str) -> np.ndarray | None:
        k = f"hap_{gene}"
        if k not in self.z:
            return None
        h = self.z[k]                       # (2 * samples, positions)
        return np.stack([h[0::2], h[1::2]], axis=2)    # (samples, positions, 2)

    def agreement(self, gene: Gene, chrom: int, tested: np.ndarray, user_call: dict) -> dict | None:
        G = self.genotypes(gene.symbol)
        if G is None:
            return None
        known = (G[self.use] != UNTESTED).all(axis=(0, 2))      # positions the panel has
        cols_full = np.flatnonzero(known)
        cols_mask = np.flatnonzero(known & tested)
        if len(cols_full) == 0:
            return None
        X = chrom == 23
        idx = np.flatnonzero(self.use)
        cache_f, cache_m = {}, {}
        rows = []
        for s in idx:
            hap = X and self.male[s]
            g = G[s]
            f_key = (hap, g[cols_full].tobytes())
            m_key = (hap, g[cols_mask].tobytes())
            if f_key not in cache_f:
                gg = np.full_like(g, UNTESTED)
                gg[cols_full] = g[cols_full]
                cache_f[f_key] = _summary(gene.call(gg, haploid=hap))
            if m_key not in cache_m:
                gg = np.full_like(g, UNTESTED)
                gg[cols_mask] = g[cols_mask]
                cache_m[m_key] = _summary(gene.call(gg, haploid=hap))
            rows.append((self.pop[s], cache_f[f_key], cache_m[m_key]))
        mine = _summary(user_call)
        same = [r for r in rows if r[2] == mine]
        allp = [r for r in rows if r[2][0] is not None]
        res = {"panel": len(rows), "positions_in_panel": int(len(cols_full)), "positions_used": int(len(cols_mask)),
               "same_call": len(same),
               "phenotype_agrees": (sum(r[1][0] == r[2][0] for r in same) / len(same)) if same else None,
               "diplotype_agrees": (sum(r[1][1] == r[2][1] for r in same) / len(same)) if same else None,
               "overall": (sum(r[1][0] == r[2][0] for r in allp) / len(allp)) if allp else None,
               "by_population": {}}
        for p in sorted({r[0] for r in rows}):
            sp = [r for r in same if r[0] == p]
            ap = [r for r in rows if r[0] == p]
            res["by_population"][p] = {"n": len(ap), "same_call": len(sp),
                                       "phenotype_agrees": (sum(r[1][0] == r[2][0] for r in sp) / len(sp)) if sp else None,
                                       "carrier_full": None}
        if same:
            others = Counter(r[1][0] for r in same if r[1][0] != mine[0])
            res["full_instead"] = [{"phenotype": k, "share": v / len(same)} for k, v in others.most_common(4)]
        return res


def _summary(c: dict) -> tuple:
    dips = tuple(sorted(c.get("diplotypes") or []))
    ph = c.get("phenotype") or ("uncertain: " + " or ".join(c["possible"]) if c.get("possible") else None)
    if c.get("genotype_only"):
        ph = "diplotype " + " or ".join(dips) if dips else None
    return (ph, dips)


# ------------------------------------------------------------------ the cards
DISCLAIMER = ("Show this to a doctor or pharmacist; never start, stop or change a medicine because of it. "
              "A consumer chip reads only some of the positions a clinical pharmacogenomic test reads, and chip calls "
              "of rare variants are often wrong (Weedon et al. 2021).")


def _pct(x):
    return "n/a" if x is None else f"{100 * x:.0f}%"


def pgx_cards(pgx: dict, panel: Panel | None, g: RawGenotypes, pos37: np.ndarray, strand: str, sex: str,
              make_card, step=None) -> list[dict]:
    gt = PgxGenotyper(g, pos37, strand)
    refs = pgx.get("references", {})
    genes = {s: Gene(s, d) for s, d in pgx["genes"].items()}
    status, calls, cards = {}, {}, []
    drugs_by_gene = {}
    for d in pgx["drugs"]:
        for gs in d["genes"]:
            drugs_by_gene.setdefault(gs, []).append(d["name"])
    for n, (sym, gene) in enumerate(sorted(genes.items())):
        if step:
            step(f"medicines: {sym}", 0.9 + 0.08 * n / max(1, len(genes)))
        d = pgx["genes"][sym]
        chrom = d["chrom"]
        haploid = chrom == 26 or (chrom == 23 and sex == "XY")
        geno = gt.genotype(gene, chrom)
        c = gene.call(geno, haploid, alternatives=True)
        tested = geno[:, 0] != UNTESTED
        rel = (panel.agreement(gene, chrom, tested, c)
               if (panel is not None and sym != "CYP2D6" and c["status"] in ("called", "ambiguous")) else None)
        cyp2d6 = sym == "CYP2D6"
        if c["status"] == "called" and not cyp2d6 and not c.get("genotype_only"):
            why = _unreliable(gene, geno, rel, [tuple(x.split("/")) for x in c["diplotypes"]], c)
            if why:
                c = {**c, "status": "unreliable", "presumptive": c["phenotype"], "phenotype": None, "why": why}
        calls[sym] = c
        if c["status"] == "called" and not cyp2d6 and c.get("phenotype"):
            status[sym] = c["score"] if d["lookup"] == "ACTIVITY_SCORE" else c["phenotype"]
        cards.append(_gene_card(sym, d, gene, c, rel, geno, drugs_by_gene.get(sym, []), refs, make_card, haploid))
    cards.insert(0, _summary_card(pgx, genes, calls, status, refs, make_card))
    return cards


MIN_AGREEMENT = 0.90


def _unreliable(gene: Gene, geno, rel, c_names=None, call=None) -> str | None:
    """Why a call from these positions should not be used, or None. With the 1000 Genomes panel: a full
    reading gives the same phenotype in fewer than 90% of people with this call. Without it (no panel
    genotypes, e.g. MT-RNR1): an allele seen in at least 5% of some population, which would change the
    phenotype in place of either called allele, was not read."""
    if rel and rel.get("phenotype_agrees") is not None:
        if rel["phenotype_agrees"] < MIN_AGREEMENT:
            alt = (rel.get("full_instead") or [{}])[0]
            return (f"a full reading gives the same phenotype in only {_pct(rel['phenotype_agrees'])} of the "
                    f"{rel['same_call']:,} 1000 Genomes people with this call"
                    + (f" (most often instead: {alt['phenotype']}, {_pct(alt['share'])})" if alt.get("phenotype") else ""))
        return None
    pp = (call or {}).get("phenotype_probabilities")
    if pp and max(pp.values()) < MIN_AGREEMENT:
        return (f"several diplotypes fit your genotypes and the likeliest phenotype carries only "
                f"{_pct(max(pp.values()))} of their combined 1000 Genomes frequency")
    dip = (c_names or [None])[0]
    if dip is None:
        return None
    cur = (gene.phenotype(dip) or {}).get("result")
    for u in _alleles_untested(gene, geno):
        if u["max_frequency"] < 0.05:
            continue
        if gene.mode == "variants":
            found = [(v["name"], v["copies"]) for v in (call or {}).get("variants", [])] + [(u["allele"], 1)]
            swaps = [gene.variant_diplotype(found, len(dip) == 1)]
        else:
            swaps = [(u["allele"],)] if len(dip) == 1 else [(u["allele"], dip[1]), (dip[0], u["allele"])]
        if any((gene.phenotype(sw) or {}).get("result") != cur for sw in swaps):
            return (f"{u['allele']} ({u['function'] or 'function not assigned'}, up to {100 * u['max_frequency']:.0f}% "
                    f"in {u['where']}) would change the result and is at positions your file does not have")
    return None


def _alleles_untested(gene: Gene, geno) -> list[dict]:
    """Alleles none of whose defining positions your file read, most frequent first."""
    tested = geno[:, 0] != UNTESTED
    out = []
    for i, a in enumerate(gene.alleles):
        if i == gene.ref or not gene.defined[i].any():
            continue
        if not (gene.defined[i] & tested).any():
            f = a.get("freq") or {}
            vals = [v for v in f.values() if isinstance(v, (int, float))]
            top = max(vals) if vals else None
            if top is not None and top >= 0.01:
                where = max(f, key=lambda k: f[k] if isinstance(f[k], (int, float)) else -1)
                out.append({"allele": a["name"], "function": a.get("function"), "max_frequency": top, "where": where})
    out.sort(key=lambda r: -r["max_frequency"])
    return out


def _rare_in_call(gene: Gene, c: dict) -> list[str]:
    names = set()
    for dp in c.get("diplotypes") or []:
        names.update(dp.split("/"))
    rare = []
    for nme in sorted(names):
        a = gene.by_name.get(nme)
        if not a or a.get("ref"):
            continue
        f = [v for v in (a.get("freq") or {}).values() if isinstance(v, (int, float))]
        if f and max(f) < 0.01:
            rare.append(nme)
    return rare


def _gene_card(sym, d, gene, c, rel, geno, drug_names, refs, make_card, haploid):
    ntot = len(gene.pos)
    tested = int((geno[:, 0] != UNTESTED).sum())
    untested = _alleles_untested(gene, geno)
    rare = _rare_in_call(gene, c)
    cyp2d6 = sym == "CYP2D6"
    drugs = ", ".join(sorted(set(drug_names))) or "none"
    ok = c["status"] in ("called", "ambiguous", "unreliable")
    if c["status"] == "not tested":
        summary = f"Your file has none of the {ntot} positions CPIC uses to define {sym} alleles."
    elif c["status"].startswith("no call"):
        summary = f"No combination of {sym} alleles fits your genotypes{c['status'][7:]}; a clinical test is needed."
    else:
        dips = c["diplotypes"]
        dtxt = dips[0] if len(dips) == 1 else " or ".join(dips[:4]) + (" (and others)" if len(dips) > 4 else "")
        if cyp2d6:
            summary = (f"Alleles {dtxt} from single-letter variants only. CYP2D6 is often deleted or duplicated, which "
                       "a chip cannot see, so no phenotype is given; a clinical CYP2D6 test (with copy number) is needed "
                       "before it could guide any medicine.")
            if c.get("phenotype"):
                summary += f" If you have the usual two copies, these alleles would mean: {c['phenotype']} (presumptive)."
        elif c["status"] == "unreliable":
            summary = (f"Not determinable from your file. From the positions it has, the call would be {dtxt} "
                       f"({c['presumptive']}), but {c['why']}. Treat it as unknown; a clinical test is needed.")
        elif c.get("genotype_only"):
            summary = (f"{sym} {dtxt}. CPIC assigns this gene no phenotype; its guideline uses the genotype directly "
                       f"(for {drugs}).")
        elif c["status"] == "called":
            summary = f"{c['phenotype']} ({sym} {dtxt})."
            if c.get("chosen_by"):
                pp = (c.get("phenotype_probabilities") or {}).get(c["phenotype"])
                summary += (f" Several diplotypes fit an unphased genotype here; this is the commonest by 1000 Genomes "
                            f"haplotype frequencies, and diplotypes giving this phenotype carry {_pct(pp)} of the total.")
        else:
            summary = f"{sym} {dtxt}: the positions your file has cannot separate {' and '.join(c['possible'])}."
        if rel and rel.get("phenotype_agrees") is not None and not cyp2d6 and c["status"] != "unreliable":
            summary += (f" Among {rel['same_call']:,} 1000 Genomes people given the same result from the same positions, "
                        f"a full reading gave the same phenotype in {_pct(rel['phenotype_agrees'])}.")
        if rare:
            summary += (f" {', '.join(rare)} is rare (under 1% in every population CPIC lists): chip calls of rare variants "
                        "are often false, so it needs confirming with a clinical test.")
    agree = (rel or {}).get("phenotype_agrees")
    grade = ("C" if cyp2d6 or not ok or c["status"] == "unreliable"
             else "A" if (agree is not None and agree >= 0.95 and not rare) else "B")
    limits = [DISCLAIMER,
              f"Your file has {tested} of the {ntot} positions CPIC uses for {sym}; an allele defined only at missing "
              "positions reads as the reference allele (*1 or 'Reference'), as in PharmCAT."]
    if untested:
        limits.append("Common alleles your file cannot see: " + "; ".join(
            f"{u['allele']} ({u['function'] or 'function not assigned'}, up to {100 * u['max_frequency']:.0f}% in {u['where']})"
            for u in untested[:5]) + ".")
    if haploid and d["chrom"] == 23:
        limits.append("Read as one X chromosome (genetic sex XY).")
    if d["chrom"] == 23 and not haploid:
        limits.append("With two X chromosomes, G6PD activity in a carrier of one deficient allele varies from cell to cell (X inactivation); CPIC calls this 'Variable'.")
    if d["chrom"] == 26:
        limits.append("Mitochondrial variants can be present in only some of your mitochondria (heteroplasmy), which a chip does not measure.")
    if sym == "RYR1":
        limits.append("Only a fraction of the RYR1 variants linked to malignant hyperthermia are on a chip; finding none does not rule out susceptibility, which is why CPIC's result without a listed variant is 'Uncertain Susceptibility'.")
    method = ("Your genotype at each CPIC-defined position (GRCh37 positions from the UCSC liftover, checked against "
              "Ensembl in CI) matched to every pair of CPIC alleles; the best-scoring pairs are reported, as in PharmCAT's "
              "Named Allele Matcher. Allele functions and the phenotype for each combination are CPIC's. ")
    if sym in VARIANT_GENES:
        method += "As CPIC specifies for this gene, the two lowest-function variants found make the diplotype. "
    method += ("Reliability: the unrelated people of the 1000 Genomes high-coverage panel, called once with every position "
               "and once with only yours. " f"CPIC level-A medicines: {drugs}.")
    res = {"gene": sym, "call": c, "reliability": rel, "untested_alleles": untested, "rare_in_call": rare,
           "positions": [{"rsid": p.get("rsid"), "hgvs": p.get("hgvs"), "name": p.get("name"),
                          "genotype": _gt_text(p, geno[k], haploid)} for k, p in enumerate(gene.pos)
                         if geno[k, 0] != UNTESTED],
           "drugs": sorted(set(drug_names))}
    srcs = [{"title": "CPIC allele definition, function and phenotype tables for " + sym,
             "url": f"https://www.clinpgx.org/page/{d.get('ref_page') or 'pgxGeneRef'}"}]
    if refs.get("pharmcat"):
        srcs.append(_cite(refs["pharmcat"]))
    if rare and refs.get("weedon2021"):
        srcs.append(_cite(refs["weedon2021"]))
    if refs.get("1000g30x"):
        srcs.append(_cite(refs["1000g30x"]))
    return make_card(f"pgx_{sym.lower().replace('-', '_')}", sym, "Medicines", grade, summary, res, limits,
                     method, srcs, ok)


def _gt_text(p, g2, haploid):
    al = p["alleles"]
    def nm(x):
        return "?" if x == NOVEL else al[x] if x < len(al) else "?"
    return nm(g2[0]) if haploid else f"{nm(g2[0])}/{nm(g2[1])}"


def _cite(r):
    return {"title": f"{r['first_author']} et al., {r['title']}, {r['journal']} ({r['year']})", "url": r["url"]}


def recommendations(pgx: dict, status: dict, calls: dict) -> list[dict]:
    """Every CPIC level-A drug, with the recommendation rows that match the
    genes we could call, and what is still unknown."""
    def same(a, b):
        return _num(a) == _num(b) or str(a).strip().lower() == str(b).strip().lower()

    out = []
    for d in pgx["drugs"]:
        unknown = [gs for gs in d["genes"] if gs not in status]
        fits = [r for r in d["recs"] if all(same(v, status[gs]) for gs, v in r["key"].items() if gs in status)]
        # a gene we could not call: CPIC's own "No Result" row when it has one
        rows = [r for r in fits if all(str(v).strip() in ("No Result", "n/a") for gs, v in r["key"].items()
                                       if gs not in status)]
        depends = bool(unknown) and not rows
        out.append({"drug": d["name"], "genes": d["genes"], "unknown": unknown, "depends": depends,
                    "guideline": d["guideline"],
                    "rows": [{"population": r.get("population"), "recommendation": r.get("recommendation"),
                              "classification": r.get("classification"), "implications": r.get("implications"),
                              "comments": r.get("comments"), "key": r["key"], "priority": r.get("priority")}
                             for r in rows]})
    return out


def _summary_card(pgx, genes, calls, status, refs, make_card):
    pr = {}
    for sym, c in calls.items():
        if sym in status:
            pr[sym] = c.get("priority") or (pgx["genes"][sym].get("phenotypes", {}).get(c.get("phenotype")) or {}).get("ehrpriority")
    # CPIC's EHR priority: "Abnormal/Priority/High Risk" (or "Priority/High Risk") marks a result that changes
    # a recommendation; the normal labels vary ("Normal/Routine/Low Risk", "Normal Risk", "none")
    flagged = sorted(s for s, p in pr.items() if p and "high risk" in p.lower())
    recs = recommendations(pgx, status, calls)
    actionable = [r for r in recs if r["rows"] and any(g in flagged for g in r["genes"])]
    needs_test = sorted({r["drug"] for r in recs if "CYP2D6" in r["genes"]})
    called = sorted(status)
    if flagged and actionable:
        s = (f"Results that change a CPIC recommendation in {', '.join(flagged)}, for "
             f"{len(actionable)} medicine{'s' if len(actionable) != 1 else ''}: "
             + ", ".join(sorted({r['drug'] for r in actionable})[:12]) + ("…" if len(actionable) > 12 else "") + ".")
    elif flagged:
        s = f"Results CPIC marks as high priority in {', '.join(flagged)}."
    else:
        s = "No result that changes a CPIC recommendation in the genes your file could read."
    unrel = sorted(k for k, c in calls.items() if c["status"] == "unreliable")
    if unrel:
        s += f" Not determinable from your file: {', '.join(unrel)}."
    s += f" Called {len(called)} of {len(genes)} genes."
    if needs_test:
        s += f" {len(needs_test)} medicines depend on CYP2D6, which needs a clinical test."
    limits = [DISCLAIMER,
              "Only CPIC level-A gene-drug pairs are included; other guidelines (DPWG, FDA labels) are not.",
              "HLA-A and HLA-B (abacavir, allopurinol, carbamazepine, oxcarbazepine, phenytoin, fosphenytoin) need HLA typing, which a chip cannot do reliably, so they are not reported.",
              "CPIC guidelines are updated; the pack records the date its data was downloaded. Check cpicpgx.org for the current recommendation."]
    res = {"flagged": flagged, "called": called, "drugs": recs, "cpic": pgx.get("source"),
           "genes": {s: {"phenotype": c.get("phenotype"), "diplotypes": c.get("diplotypes"), "status": c["status"],
                         "priority": pr.get(s)} for s, c in calls.items()}}
    srcs = [{"title": "CPIC: Clinical Pharmacogenetics Implementation Consortium guidelines (CC0)",
             "url": "https://cpicpgx.org/guidelines/"}]
    if refs.get("relling2011"):
        srcs.append(_cite(refs["relling2011"]))
    if refs.get("weedon2021"):
        srcs.append(_cite(refs["weedon2021"]))
    return make_card("medicines", "Medicines and your genes", "Medicines", "B", s, res, limits,
                     "Each gene is called from your file (see its card); the phenotype then selects the CPIC "
                     "recommendation for each medicine. Medicines whose guideline also needs a gene your file cannot "
                     "settle show that gene as unknown.", srcs, bool(called))
