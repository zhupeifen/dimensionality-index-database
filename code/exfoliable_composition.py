# -*- coding: utf-8 -*-
"""Overlap of the two-dimensional sublattices with the Mounet exfoliable set, by composition.

exfoliable_overlap.py matches by COD entry number, which undercounts: the
Mounet table gives one source entry per compound, so other depositions of the
same compound go unmatched. Matching by reduced composition overcounts, since
it also catches polymorphs of a listed compound. The two together bound the
overlap from both sides.

    python exfoliable_composition.py <cod_dimensionality.jsonl[.gz]> <EE_and_PE_structures.txt> [out.json]

EE_and_PE_structures.txt is distributed by Mounet et al., Materials Cloud,
doi:10.24435/materialscloud:2017.0008/v3, and is not redistributed here.
"""
import gzip
import json
import sys
import warnings

from pymatgen.core import Composition

warnings.filterwarnings("ignore")


def reduced(formula):
    try:
        return Composition(formula).reduced_formula
    except Exception:
        return None


def main():
    db, table = sys.argv[1], sys.argv[2]
    out = sys.argv[3] if len(sys.argv) > 3 else "exfoliable_composition.json"

    mounet = set()
    with open(table, encoding="utf-8") as f:
        next(f)                                    # header
        for line in f:
            parts = line.split()
            if len(parts) > 3:                     # FORMULA_3D is the fourth column
                r = reduced(parts[3])
                if r:
                    mounet.add(r)

    opener = gzip.open if db.endswith(".gz") else open
    twod = matched = 0
    comps, matched_comps = set(), set()
    with opener(db, "rt", encoding="utf-8") as f:
        for line in f:
            rec = json.loads(line)
            if rec.get("d_top") != 2:
                continue
            twod += 1
            r = reduced(rec.get("formula", ""))
            comps.add(r)
            if r in mounet:
                matched += 1
                matched_comps.add(r)

    res = {"mounet_compositions": len(mounet), "twod_structures": twod,
           "matched_structures": matched,
           "pct_of_twod": round(100.0 * matched / twod, 1),
           "twod_compositions": len(comps), "matched_compositions": len(matched_comps)}
    json.dump(res, open(out, "w", encoding="utf-8"), indent=1)
    print(res)


if __name__ == "__main__":
    main()
