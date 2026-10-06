# -*- coding: utf-8 -*-
"""Why a third of the anion-choice sample was never evaluated.

anion_choice.py draws 1,500 structures, with a fixed seed, from those set aside
as "anion class not defined", and evaluates the ones that carry two or more
anion classes. The "anion class not defined" skip, however, fires whenever the
number of anion classes present is not exactly one, so it also holds structures
with a framework metal and no O, S, Se, Te or halogen at all. This script redraws
the identical sample and records, for each structure, how many anion classes are
present, or why it could not be read.

    python anion_choice_dropouts.py <cif_root_or_cache> [out.json]

CIFs are read from <cif_root_or_cache>/<cod>.cif, or from the COD tree layout.
"""
import collections
import gzip
import json
import os
import random
import sys
import warnings

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
SEED, N = 20260914, 1500                 # as in anion_choice.py, run with N = 1500
ANION_CUT = {"halide": {"Cl", "Br", "I", "F"}, "chalcogen": {"S", "Se", "Te"}, "oxide": {"O"}}


def db_path():
    for p in (os.path.join(HERE, "data", "cod_dimensionality.jsonl.gz"),
              os.path.join(HERE, "..", "data", "cod_dimensionality.jsonl.gz"),
              os.path.join(HERE, "deposit_npj", "data", "cod_dimensionality.jsonl.gz")):
        if os.path.exists(p):
            return p
    raise SystemExit("cod_dimensionality.jsonl.gz not found")


def cif(root, cod):
    for p in (os.path.join(root, cod + ".cif"),
              os.path.join(root, cod[0], cod[1:3], cod[3:5], cod + ".cif")):
        if os.path.exists(p):
            return p
    return None


def main():
    from pymatgen.core import Structure
    root = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else "anion_choice_dropouts.json"
    pool = []
    with gzip.open(db_path(), "rt", encoding="utf-8") as fh:
        for line in fh:
            r = json.loads(line)
            if r.get("skip") == "anion class not defined":
                pool.append(r["cod"])
    random.seed(SEED)
    random.shuffle(pool)
    sample = pool[:N]

    why = collections.Counter()
    per = {}
    for cod in sample:
        p = cif(root, cod)
        if p is None:
            why["cif not available"] += 1
            per[cod] = "cif not available"
            continue
        try:
            st = Structure.from_file(p)
            els = set()
            for site in st.sites:            # as cod_wide_index.symbols(): majority species
                try:
                    els.add(site.specie.symbol)
                except Exception:
                    els.add(max(site.species, key=site.species.get).symbol)
        except Exception:
            why["unreadable"] += 1
            per[cod] = "unreadable"
            continue
        n = sum(1 for cuts in ANION_CUT.values() if els & cuts)
        key = {0: "no anion class", 1: "one anion class"}.get(n, "two or more classes")
        why[key] += 1
        per[cod] = key
    res = {"sampled": len(sample), "reasons": dict(why), "per_structure": per}
    json.dump(res, open(out, "w", encoding="utf-8"), indent=1)
    print(json.dumps(res["reasons"], indent=1))


if __name__ == "__main__":
    main()
