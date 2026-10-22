"""Haplogroups from genotypes: mtDNA (PhyloTree 17) and Y chromosome (YFull).

Trees:
  mtDNA  PhyloTree Build 17 rooted at the Reconstructed Sapiens Reference
         Sequence (van Oven & Kayser, Hum. Mutat. 30, E386, 2009; Behar et al.,
         Am. J. Hum. Genet. 90, 675, 2012), as distributed for Haplogrep 3
         (Schönherr et al., Nucleic Acids Res. 51, W263, 2023; MIT licence),
         with Haplogrep's mutation weights (rare, stable mutations weigh more
         than recurrent hotspots).
  Y      the YFull tree with GRCh37 positions as distributed with Yleaf
         (Ralf et al., Mol. Biol. Evol. 35, 1291, 2018; GPL-3).

Calling: every branch is defined by mutations (position, ancestral base,
derived base). For each branch, the genotyped positions give weighted counts
of derived (D) and ancestral (A) calls. A branch is a candidate when it has
more derived than ancestral support; the call is the candidate whose path
from the root has the highest score (sum of D minus twice the sum of A along
the path), so a single stray derived call on another branch cannot outweigh
a consistent lineage. The report lists the path and the counts behind each
step so the call can be checked.
"""
from __future__ import annotations

import re
import xml.etree.ElementTree as ET

import numpy as np

MUT = re.compile(r"^(\d+)([ACGT])(!*)$")


# ------------------------------------------------------------------ building (data pack)
def build_mt_tree(tree_xml: str, root_sequence: str, weights_text: str) -> dict:
    """PhyloTree XML (RSRS-rooted) -> nodes with (pos, ancestral, derived) mutations.
    Ancestral bases come from walking the tree from the RSRS sequence."""
    weights = {}
    for line in weights_text.splitlines():
        f = line.split("\t")
        if len(f) >= 2:
            try:
                weights[f[0]] = float(f[1])
            except ValueError:
                pass
    root = ET.fromstring(tree_xml)
    nodes: list[dict] = []

    def walk(el, parent: int, state: dict):
        name = el.get("name")
        muts = []
        st = dict(state)
        for p in el.findall("./details/poly"):
            txt = (p.text or "").strip()
            m = MUT.match(txt)
            if not m:
                continue                     # insertions, deletions: not typed on genotyping chips
            pos, der = int(m.group(1)), m.group(2)
            anc = st.get(pos, root_sequence[pos - 1] if pos <= len(root_sequence) else "N")
            if anc != der:
                key = f"{pos}{der}"
                muts.append([pos, anc, der, round(weights.get(key, weights.get(txt.rstrip("!"), 5.0)), 2)])
            st[pos] = der
        idx = len(nodes)
        nodes.append({"name": name, "parent": parent, "muts": muts})
        for ch in el.findall("./haplogroup"):
            walk(ch, idx, st)

    for top in root.findall("./haplogroup"):
        walk(top, -1, {})
    return {"tree": "PhyloTree 17 (RSRS)", "nodes": nodes}


def build_y_tree(tree_json: dict, positions_text: str, root_name: str | None = None) -> dict:
    """Yleaf's YFull tree (parent -> children) and SNP table (GRCh37):
    columns chr, SNP name, haplogroup, position, mutation, ancestral, derived."""
    children = tree_json
    all_children = {c for v in children.values() for c in v}
    roots = [k for k in children if k not in all_children]
    root = root_name or roots[0]
    nodes: list[dict] = []
    index: dict = {}

    stack = [(root, -1)]
    while stack:
        name, parent = stack.pop()
        if name in index:
            continue
        index[name] = len(nodes)
        nodes.append({"name": name, "parent": parent, "muts": []})
        for c in reversed(children.get(name, [])):
            stack.append((c, index[name]))
    seen = set()
    for line in positions_text.splitlines():
        f = line.split("\t")
        if len(f) < 7:
            continue
        name, hg, pos, anc, der = f[1], f[2], f[3], f[5].upper(), f[6].upper()
        if hg not in index or not pos.isdigit() or len(anc) != 1 or len(der) != 1 or anc == der:
            continue
        key = (int(pos), der)
        if key in seen:
            continue            # the same mutation listed under several names
        seen.add(key)
        nodes[index[hg]]["muts"].append([int(pos), anc, der, 1.0, name])
    return {"tree": "YFull (via Yleaf)", "nodes": nodes}


# ------------------------------------------------------------------ calling
class Tree:
    def __init__(self, data: dict):
        self.nodes = data["nodes"]
        self.name = data.get("tree", "")
        self.parent = np.array([n["parent"] for n in self.nodes], dtype=np.int64)
        self.index = {n["name"]: i for i, n in enumerate(self.nodes)}
        depth = np.zeros(len(self.nodes), dtype=np.int64)
        for i, p in enumerate(self.parent):        # parents precede children in the list
            depth[i] = depth[p] + 1 if p >= 0 else 0
        self.depth = depth
        self.positions = sorted({m[0] for n in self.nodes for m in n["muts"]})

    def path(self, i: int) -> list[int]:
        out = []
        while i >= 0:
            out.append(i)
            i = int(self.parent[i])
        return out[::-1]

    def is_ancestor(self, a: int, b: int) -> bool:
        """a is b or an ancestor of b."""
        while b >= 0:
            if a == b:
                return True
            b = int(self.parent[b])
        return False

    def call(self, observed: dict[int, str]) -> dict:
        """observed: position -> base (on the + strand of the tree's reference)."""
        n = len(self.nodes)
        D, A = np.zeros(n), np.zeros(n)
        nd, na = np.zeros(n, dtype=np.int64), np.zeros(n, dtype=np.int64)
        typed = 0
        for i, node in enumerate(self.nodes):
            for m in node["muts"]:
                b = observed.get(m[0])
                if b is None:
                    continue
                typed += 1
                w = m[3]
                if b == m[2]:
                    D[i] += w; nd[i] += 1
                elif b == m[1]:
                    A[i] += w; na[i] += 1
        PD, PA = np.zeros(n), np.zeros(n)
        for i in range(n):
            p = self.parent[i]
            PD[i] = D[i] + (PD[p] if p >= 0 else 0)
            PA[i] = A[i] + (PA[p] if p >= 0 else 0)
        cand = np.flatnonzero((D > 0) & (D > A))
        if typed == 0:
            return {"call": None, "reason": "none of the tree's positions are genotyped in this file", "typed_mutations": 0}
        if len(cand) == 0:
            return {"call": self.nodes[0]["name"], "reason": "no branch has more derived than ancestral calls",
                    "typed_mutations": typed, "path": []}
        score = PD[cand] - 2 * PA[cand]
        order = np.lexsort((-self.depth[cand], -score))
        best = int(cand[order[0]])
        runner = [int(c) for c in cand[order[1:6]] if not self.is_ancestor(int(c), best)]
        path = self.path(best)
        return {"call": self.nodes[best]["name"], "node": best,
                "typed_mutations": typed,
                "path": [{"name": self.nodes[i]["name"], "derived": int(nd[i]), "ancestral": int(na[i]),
                          "defining_mutations": len(self.nodes[i]["muts"])} for i in path],
                "path_derived": int(nd[path].sum()), "path_ancestral": int(na[path].sum()),
                "alternatives": [{"name": self.nodes[c]["name"], "derived": int(nd[c]), "ancestral": int(na[c]),
                                  "score": float(PD[c] - 2 * PA[c])} for c in runner[:3]],
                "score": float(PD[best] - 2 * PA[best])}

    def lookup(self, name: str, trim: bool = True) -> int | None:
        """Find a node by name. With ``trim`` (mtDNA names such as
        'H1c3a+16189' or 'U5b1a1*'), extra notation and then trailing
        characters are removed until a node matches."""
        if not name:
            return None
        name = name.strip()
        if name in self.index:
            return self.index[name]
        if not trim:
            return None
        base = re.split(r"[+\s(*]", name)[0]
        while base:
            if base in self.index:
                return self.index[base]
            base = base[:-1]
        return None
