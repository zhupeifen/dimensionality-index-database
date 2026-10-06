# -*- coding: utf-8 -*-
"""Does the majority-species reading of partially occupied sites move the rank?

Uses the same seeded sample of indexed structures as partial_occupancy_sample.py.
Two checks:
  1. the rank distribution of the ordered and the disordered structures in the
     sample, side by side;
  2. for every disordered structure, the integer rank recomputed after removing
     each site whose total occupancy is below 0.5, the minor partners of split
     sites and sparsely occupied positions that the majority-species reading
     counts as full atoms. Sites at 0.5 and above are kept, so a site shared
     equally between two positions survives.
The rank is computed with cod_wide_index.index_structure, as in the database.

    python partial_occupancy_test.py <cod_dimensionality.jsonl.gz> <cif_dir> [n] [out.json]
"""
import gzip
import importlib.util
import json
import os
import random
import sys
import warnings
from collections import Counter

warnings.filterwarnings("ignore")
SEED = 20261006
HERE = os.path.dirname(os.path.abspath(__file__))


def indexer():
    for p in (os.path.join(HERE, "cod_wide_index.py"),
              os.path.join(HERE, "deposit_npj", "code", "cod_wide_index.py")):
        if os.path.exists(p):
            spec = importlib.util.spec_from_file_location("cwi", p)
            m = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(m)
            return m
    raise SystemExit("cod_wide_index.py not found")


def cif(root, cod):
    for p in (os.path.join(root, cod + ".cif"),
              os.path.join(root, cod[0], cod[1:3], cod[3:5], cod + ".cif")):
        if os.path.exists(p):
            return p
    return None


def main():
    from pymatgen.core import Structure
    cwi = indexer()
    db, root = sys.argv[1], sys.argv[2]
    n = int(sys.argv[3]) if len(sys.argv) > 3 else 1000
    out = sys.argv[4] if len(sys.argv) > 4 else "partial_occupancy_test.json"
    rank = {}
    with gzip.open(db, "rt", encoding="utf-8") as f:
        for line in f:
            r = json.loads(line)
            if "skip" not in r:
                rank[r["cod"]] = r["d_top"]
    random.seed(SEED)
    sample = random.sample(list(rank), n)

    dist = {"ordered": Counter(), "disordered": Counter()}
    moved, tested, untestable, examples, reindexed = 0, 0, 0, [], 0
    for cod in sample:
        st = Structure.from_file(cif(root, cod))
        grp = "ordered" if st.is_ordered else "disordered"
        dist[grp][rank[cod]] += 1
        if grp == "ordered":
            continue
        keep = [i for i, s in enumerate(st.sites) if s.species.num_atoms >= 0.5]
        if len(keep) == len(st):
            tested += 1                     # nothing below 0.5: same reading
            continue
        try:
            res, _ = cwi.index_structure(Structure.from_sites([st[i] for i in keep]))
        except Exception:                   # noqa: BLE001
            res = None
        if not res:
            untestable += 1                 # e.g. the anion class itself was minor
            continue
        tested += 1
        reindexed += 1
        if res["d_top"] != rank[cod]:
            moved += 1
            if len(examples) < 10:
                examples.append([cod, rank[cod], res["d_top"]])

    def pct(c):
        n = sum(c.values())
        return {"n": n, **{str(k): round(100 * c[k] / n, 1) for k in range(4)}}
    res = {"sampled": n, "distribution_pct": {k: pct(v) for k, v in dist.items()},
           "disordered_tested": tested, "disordered_reindexed_without_minor_sites": reindexed,
           "disordered_untestable": untestable,
           "rank_moved": moved,
           "rank_moved_pct_of_tested": round(100 * moved / tested, 1) if tested else None,
           "examples": examples}
    json.dump(res, open(out, "w", encoding="utf-8"), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k != "examples"}, indent=1))
    print("examples", examples)


if __name__ == "__main__":
    main()
