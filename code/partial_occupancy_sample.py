# -*- coding: utf-8 -*-
"""How many indexed structures carry partially occupied sites?

The indexer reads every site as its majority species (cod_wide_index.symbols),
without enumerating orderings. This measures how often that reading is needed:
a seeded random sample of indexed structures, each parsed with pymatgen and
counted as disordered if any site is not fully occupied by one species.

    python partial_occupancy_sample.py <cod_dimensionality.jsonl.gz> <cif_dir> [n] [out.json]

<cif_dir> holds <cod>.cif for the sampled entries (or a COD tree).
"""
import gzip
import json
import math
import os
import random
import sys
import warnings

warnings.filterwarnings("ignore")
SEED = 20261006


def cif(root, cod):
    for p in (os.path.join(root, cod + ".cif"),
              os.path.join(root, cod[0], cod[1:3], cod[3:5], cod + ".cif")):
        if os.path.exists(p):
            return p
    return None


def main():
    from pymatgen.core import Structure
    db, root = sys.argv[1], sys.argv[2]
    n = int(sys.argv[3]) if len(sys.argv) > 3 else 1000
    out = sys.argv[4] if len(sys.argv) > 4 else "partial_occupancy_sample.json"
    ids = []
    with gzip.open(db, "rt", encoding="utf-8") as f:
        for line in f:
            r = json.loads(line)
            if "skip" not in r:
                ids.append(r["cod"])
    random.seed(SEED)
    sample = random.sample(ids, n)
    disordered = ordered = failed = 0
    for cod in sample:
        p = cif(root, cod)
        try:
            st = Structure.from_file(p)
        except Exception:                                # noqa: BLE001
            failed += 1
            continue
        if st.is_ordered:
            ordered += 1
        else:
            disordered += 1
    m = disordered + ordered
    f = disordered / m
    half = 1.96 * math.sqrt(f * (1 - f) / m)
    res = {"sampled": n, "parsed": m, "failed": failed, "disordered": disordered,
           "pct_disordered": round(100 * f, 1),
           "ci95_pct": [round(100 * (f - half), 1), round(100 * (f + half), 1)]}
    json.dump(res, open(out, "w", encoding="utf-8"), indent=1)
    print(res)


if __name__ == "__main__":
    main()
