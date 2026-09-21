# -*- coding: utf-8 -*-
"""The ceiling on any model that predicts connectivity from composition.

A composition model returns one rank per formula. Where a formula's own
depositions carry different ranks, such a model is wrong about all but the
commonest of them, however good it becomes. That is arithmetic rather than a
benchmark, and it is the one limit a better estimator cannot move.

Adds a "ceiling" block to composition_ml.json in place.

    python composition_ceiling.py [cod_dimensionality.jsonl.gz]
"""
import collections
import gzip
import io
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def data(name):
    for p in (os.path.join(HERE, "deposit_npj", "data", name),
              os.path.join(HERE, os.pardir, "data", name),
              os.path.join(HERE, name)):
        if os.path.exists(p):
            return p
    return os.path.join(HERE, name)


SRC = sys.argv[1] if len(sys.argv) > 1 else data("cod_dimensionality.jsonl.gz")
OUT = data("composition_ml.json")


def main():
    by = collections.defaultdict(collections.Counter)
    with gzip.open(SRC, "rt", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            r = json.loads(line)
            if r.get("d_top") is None:
                continue
            by[r["formula"]][r["d_top"]] += 1

    total = sum(sum(c.values()) for c in by.values())
    multi = {f: c for f, c in by.items() if len(c) > 1}
    deps = sum(sum(c.values()) for c in multi.values())
    forced = sum(sum(c.values()) - max(c.values()) for c in multi.values())

    J = json.load(io.open(OUT, encoding="utf-8"))
    J["ceiling"] = {
        "indexed": total,
        "compositions_disagreeing": len(multi),
        "depositions_affected": deps,
        "necessarily_wrong": forced,
        "ceiling_overall_pct": round(100.0 * (total - forced) / total, 1),
        "ceiling_on_affected_pct": round(100.0 * (deps - forced) / deps, 1),
        "note": ("a one-label-per-formula model cannot exceed these, whatever "
                 "features or estimator it uses"),
    }
    json.dump(J, io.open(OUT, "w", encoding="utf-8", newline="\n"), indent=1)

    c = J["ceiling"]
    print("compositions whose depositions disagree : %d" % c["compositions_disagreeing"])
    print("depositions belonging to them           : %d" % c["depositions_affected"])
    print("necessarily mislabelled by any such model: %d" % c["necessarily_wrong"])
    print("ceiling over the whole archive          : %.1f per cent" % c["ceiling_overall_pct"])
    print("ceiling on the affected depositions     : %.1f per cent" % c["ceiling_on_affected_pct"])
    print("\nupdated " + os.path.basename(OUT))


if __name__ == "__main__":
    main()
