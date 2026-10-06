# -*- coding: utf-8 -*-
"""Complete the per-deposition records behind the Figure 5 classification.

npj_disagreements.json carries, for each composition deposited at more than one
rank, a record per deposition: [cod, rank, space group number, volume per atom,
sites]. The run that first built it left out 20 depositions, among them every
deposition of five compositions, so 612 of the 617 compositions could be
classified. This adds the missing records from the CIFs, computed the way the
existing ones were (checked on seven of them): the structure as parsed by
pymatgen, its space group from SpacegroupAnalyzer at the default tolerance, and
its cell volume per site.

    python complete_disagreement_records.py <cod_dimensionality.jsonl.gz> <cif_dir> <npj_disagreements.json>

<cif_dir> holds <cod>.cif for the missing depositions (or a COD tree). The file
is rewritten in place with every composition and deposition present; run
make_summaries.py afterwards.
"""
import gzip
import json
import os
import sys
import warnings
from collections import defaultdict

warnings.filterwarnings("ignore")


def cif(root, cod):
    for p in (os.path.join(root, cod + ".cif"),
              os.path.join(root, cod[0], cod[1:3], cod[3:5], cod + ".cif")):
        if os.path.exists(p):
            return p
    return None


def space_group(st, path):
    """Space group number: pymatgen at the default tolerance, as for the existing
    records; where that cannot settle the symmetry, a looser tolerance (0.1 A);
    failing that, the number the CIF itself states."""
    import re
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    for prec in (0.01, 0.1):
        try:
            return SpacegroupAnalyzer(st, symprec=prec).get_space_group_number()
        except Exception:                                # noqa: BLE001
            pass
    text = open(path, encoding="utf-8", errors="ignore").read()
    m = re.search(r"^_(?:space_group_IT_number|symmetry_Int_Tables_number)\s+(\d+)", text, re.M)
    if m:
        return int(m.group(1))
    raise ValueError("no space group")


def main():
    from pymatgen.core import Structure
    db, root, path = sys.argv[1], sys.argv[2], sys.argv[3]
    byf = defaultdict(list)
    with gzip.open(db, "rt", encoding="utf-8") as f:
        for line in f:
            r = json.loads(line)
            if "skip" not in r:
                byf[r["formula"]].append(r)
    disagree = {f: v for f, v in byf.items() if len(v) > 1 and len({r["d_top"] for r in v}) > 1}

    doc = json.load(open(path, encoding="utf-8"))
    detail = {x["formula"]: x for x in doc["detail"]}
    added, failed = 0, []
    for f, recs in sorted(disagree.items()):
        entry = detail.setdefault(f, {"formula": f, "cls": recs[0]["anion_class"], "recs": []})
        have = {r[0] for r in entry["recs"]}
        for r in recs:
            if r["cod"] in have:
                continue
            p = cif(root, r["cod"])
            try:
                st = Structure.from_file(p)
                sg = space_group(st, p)
                entry["recs"].append([r["cod"], r["d_top"], sg, round(st.volume / len(st), 3),
                                      len(st)])
                added += 1
            except Exception as e:                      # noqa: BLE001
                failed.append((r["cod"], type(e).__name__))
    for x in detail.values():
        x.pop("category", None)                          # make_summaries assigns it
    doc["detail"] = [detail[f] for f in sorted(detail)]
    doc["n"] = len(doc["detail"])
    json.dump(doc, open(path, "w", encoding="utf-8"), indent=1)
    print("compositions", len(doc["detail"]), "of", len(disagree), "| records added", added,
          "| failed", failed)


if __name__ == "__main__":
    main()
