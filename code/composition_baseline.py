# -*- coding: utf-8 -*-
"""Can composition stand in for connectivity?

A structural descriptor has to earn its place against the cheapest thing that
could replace it. For connectivity the cheapest thing is stoichiometry: A2MX4
is commonly layered and A4MX6 commonly holds isolated octahedra, so a lookup
from the anion-to-metal ratio might do the work without any graph analysis.

This measures how far that goes, over the whole indexed database. The lookup is
fitted and scored on the same structures, which is the most favourable number
such a rule can reach; anything it cannot do here it cannot do at all.

    python composition_baseline.py [cod_dimensionality.jsonl.gz]

Writes composition_baseline.json and prints the table quoted in the paper.
"""
import collections
import gzip
import io
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    HERE, "deposit_npj", "data", "cod_dimensionality.jsonl.gz")
OUT = os.path.join(HERE, "composition_baseline.json")

# Light elements an organic cation is built from, plus the bridging anions
# themselves: what is left is the framework metal.
NONMETAL = {"H", "D", "C", "N", "O", "S", "Se", "Te", "P", "F",
            "Cl", "Br", "I", "B", "Si", "As"}
ANION = {"oxide": {"O"}, "halide": {"F", "Cl", "Br", "I"},
         "chalcogen": {"S", "Se", "Te"}}
BIN = 0.5           # X/M is binned to a half, so neighbouring stoichiometries pool


def ratio(formula, cls):
    """Bridging anions per framework metal, straight from the formula."""
    anion = metal = 0
    for sym, num in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        if not sym:
            continue
        n = int(num) if num else 1
        if sym in ANION[cls]:
            anion += n
        elif sym not in NONMETAL:
            metal += n
    return (anion / metal) if metal else None


def main():
    rows = []
    # The baseline can only be scored where a ratio can be read off the
    # deposited formula, which is a slightly smaller set than the indexed one.
    # Both totals are recorded so the paper can say which it is quoting.
    n_indexed = n_no_ratio = 0
    with gzip.open(SRC, "rt", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            r = json.loads(line)
            if r.get("d_top") is None:          # not indexed
                continue
            n_indexed += 1
            q = ratio(r["formula"], r["anion_class"])
            if q is None:
                n_no_ratio += 1
                continue
            rows.append((r["anion_class"], round(q / BIN) * BIN, r["d_top"]))

    tab = collections.defaultdict(collections.Counter)
    for cls, q, d in rows:
        tab[(cls, q)][d] += 1
    majority = {k: v.most_common(1)[0][0] for k, v in tab.items()}

    hit = sum(1 for cls, q, d in rows if majority[(cls, q)] == d)
    commonest = collections.Counter(d for _, _, d in rows).most_common(1)[0]

    per_rank = {}
    for k in (0, 1, 2, 3):
        g = [x for x in rows if x[2] == k]
        per_rank[k] = (sum(1 for cls, q, d in g if majority[(cls, q)] == d), len(g))
    extended = [x for x in rows if x[2] != 0]
    ext_hit = sum(1 for cls, q, d in extended if majority[(cls, q)] == d)

    res = {
        "n": len(rows),
        "indexed_total": n_indexed,
        "no_ratio": n_no_ratio,
        "lookup_pct": round(100.0 * hit / len(rows), 1),
        "commonest_rank": commonest[0],
        "commonest_pct": round(100.0 * commonest[1] / len(rows), 1),
        "per_rank": {str(k): {"recalled": v[0], "n": v[1],
                              "pct": round(100.0 * v[0] / v[1], 1)}
                     for k, v in per_rank.items()},
        "extended_recalled": ext_hit,
        "extended_n": len(extended),
        "extended_pct": round(100.0 * ext_hit / len(extended), 1),
        "ranks_ever_predicted": sorted(set(majority.values())),
        "by_class": {},
    }
    for cls in sorted(ANION):
        g = [x for x in rows if x[0] == cls]
        h = sum(1 for c, q, d in g if majority[(c, q)] == d)
        res["by_class"][cls] = {"n": len(g), "pct": round(100.0 * h / len(g), 1)}

    json.dump(res, io.open(OUT, "w", encoding="utf-8", newline="\n"), indent=1)

    print("structures with a usable X/M ratio: %d" % res["n"])
    print("composition lookup, best case     : %.1f per cent" % res["lookup_pct"])
    print("always predicting %dD              : %.1f per cent"
          % (res["commonest_rank"], res["commonest_pct"]))
    print("\n rank        n    recalled      per cent")
    for k in ("0", "1", "2", "3"):
        v = res["per_rank"][k]
        print("  %sD  %9d %10d %11.1f" % (k, v["n"], v["recalled"], v["pct"]))
    print("\nany extended sublattice: %d of %d = %.1f per cent"
          % (res["extended_recalled"], res["extended_n"], res["extended_pct"]))
    print("ranks the rule ever predicts:", res["ranks_ever_predicted"])
    print("wrote " + os.path.basename(OUT))


if __name__ == "__main__":
    main()
