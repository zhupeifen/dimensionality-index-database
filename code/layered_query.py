# -*- coding: utf-8 -*-
"""A worked query: layered compositions the database is confident about.

The paper claims that an indexed database turns "find me a layered compound"
from a literature search into a query. This is that query, run.

Three conditions, each doing a different job. The sublattice must be
two-dimensional. Its fractional term must be below 0.05, so that no near-miss
contact would promote it under a slightly more generous bridging criterion. And
the composition must have been deposited at least twice, with every deposition
of it indexing two-dimensional - independent determinations agreeing, which is
the closest thing a database offers to a replicate.

    python layered_query.py [cod_dimensionality.jsonl.gz]

Writes layered_query.json.
"""
import collections
import gzip
import io
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

def data(name):
    """Find a data file whether run from the working tree or from the deposit."""
    for p in (os.path.join(HERE, "deposit_npj", "data", name),
              os.path.join(HERE, os.pardir, "data", name),
              os.path.join(HERE, name)):
        if os.path.exists(p):
            return p
    return os.path.join(HERE, name)

SRC = sys.argv[1] if len(sys.argv) > 1 else data("cod_dimensionality.jsonl.gz")
OUT = os.path.join(HERE, "layered_query.json")

FRAC_MAX = 0.05
MIN_DEPOSITIONS = 2
N_EXAMPLES = 25


def main():
    by_formula = collections.defaultdict(list)
    n_2d = n_firm = 0
    with gzip.open(SRC, "rt", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            r = json.loads(line)
            if r.get("d_top") is None:
                continue
            if r["d_top"] == 2:
                n_2d += 1
                if (r.get("frac") or 0.0) < FRAC_MAX:
                    n_firm += 1
            by_formula[r["formula"]].append(r)

    # The query as stated: firm two-dimensional records first, then the
    # composition must carry at least two of them, then no deposition of that
    # composition may index anything other than two dimensions.
    repeat, unanimous, by_class, rows = 0, 0, collections.Counter(), []
    for formula, recs in by_formula.items():
        firm = [r for r in recs
                if r["d_top"] == 2 and (r.get("frac") or 0.0) < FRAC_MAX]
        if len(firm) < MIN_DEPOSITIONS:
            continue
        repeat += 1
        if not all(r["d_top"] == 2 for r in recs):
            continue
        unanimous += 1
        cls = firm[0]["anion_class"]
        by_class[cls] += 1
        rows.append({"formula": formula, "depositions": len(recs), "cls": cls})

    rows.sort(key=lambda r: (-r["depositions"], r["formula"]))
    res = {"n_2d": n_2d, "n_firm": n_firm, "n_repeat": repeat,
           "n_unanimous": unanimous, "by_class": dict(by_class),
           "examples": rows[:N_EXAMPLES],
           # the whole list, for the supplementary table; it is the paper's
           # most directly usable output and a count alone cannot be checked
           "all": rows}
    json.dump(res, io.open(OUT, "w", encoding="utf-8", newline="\n"), indent=1)

    print("two-dimensional sublattices            : %d" % n_2d)
    print("  of which the term is below %.2f       : %d" % (FRAC_MAX, n_firm))
    print("compositions with %d+ firm 2D depositions : %d" % (MIN_DEPOSITIONS, repeat))
    print("  and no deposition indexing otherwise  : %d" % unanimous)
    print("     " + ", ".join("%s %d" % (k, v) for k, v in by_class.most_common()))
    print("\nmost corroborated:")
    for r in rows[:8]:
        print("   %-24s %3d depositions  %s" % (r["formula"], r["depositions"], r["cls"]))
    print("\nwrote " + os.path.basename(OUT))


if __name__ == "__main__":
    main()
