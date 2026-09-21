# -*- coding: utf-8 -*-
"""How much of the two-dimensional set is what an exfoliation screen returns?

The paper argues that a two-dimensional metal-anion sublattice is not the same
thing as an exfoliable material: exfoliability needs a weakly bound gap between
the layers, a property of the complete structure, while a two-dimensional
sublattice needs the electronically active network confined to a plane, whatever
holds the planes together. That argument can be replaced by a number.

Mounet and co-workers screened the experimental databases for easily and
potentially exfoliable compounds and published the list. Their table names the
source database and entry number for each, so the COD entries can be matched
against the index built here directly, without any structure matching.

Input is the table EE_and_PE_structures.txt distributed with that work
(Materials Cloud, doi:10.24435/materialscloud:2017.0008/v3).

    python exfoliable_overlap.py [EE_and_PE_structures.txt]

Writes exfoliable_overlap.json.
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

TABLE = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    HERE, "EE_and_PE_structures.txt")
INDEXED = data("cod_dimensionality.jsonl.gz")
OUT = os.path.join(HERE, "exfoliable_overlap.json")


def read_table(path):
    """SOURCE_DB, DB_ID and EE/PE off the fixed-width table, by column name."""
    rows = []
    with io.open(path, encoding="utf-8") as fh:
        head = fh.readline().split()
        i_db, i_id, i_ee = (head.index("SOURCE_DB"), head.index("DB_ID"),
                            head.index("EE/PE"))
        for line in fh:
            f = line.split()
            if len(f) <= max(i_db, i_id, i_ee):
                continue
            rows.append((f[i_db], f[i_id], f[i_ee]))
    return rows


def main():
    rows = read_table(TABLE)
    cod = {db_id for db, db_id, _ in rows if db.upper() == "COD"}
    cod_ee = {db_id for db, db_id, ee in rows
              if db.upper() == "COD" and ee.upper() == "EE"}
    print("Mounet et al.: %d exfoliable entries, %d unique COD ids, %d easily "
          "exfoliable" % (len(rows), len(cod), len(cod_ee)))

    rank, twod = {}, 0
    n_indexed = 0
    with gzip.open(INDEXED, "rt", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            r = json.loads(line)
            if r.get("d_top") is None:
                continue
            n_indexed += 1
            if r["d_top"] == 2:
                twod += 1
            if r["cod"] in cod:
                rank[r["cod"]] = r["d_top"]

    dist = collections.Counter(rank.values())
    overlap = dist[2]
    res = {"mounet_total": len(rows),
           "mounet_cod": len(cod),
           "mounet_cod_ee": len(cod_ee),
           "indexed_here": n_indexed,
           "twod_here": twod,
           "mounet_cod_indexed_here": len(rank),
           "rank_of_mounet_here": {str(k): dist[k] for k in sorted(dist)},
           "overlap": overlap,
           "overlap_pct_of_twod": round(100.0 * overlap / twod, 1),
           "overlap_pct_of_mounet": round(100.0 * overlap / len(cod), 1)}
    json.dump(res, io.open(OUT, "w", encoding="utf-8", newline="\n"), indent=1)

    print("  of their %d COD entries, %d are indexed here" % (len(cod), len(rank)))
    print("  sublattice rank of those: " +
          ", ".join("%dD %d" % (k, dist[k]) for k in sorted(dist)))
    print("  overlap: %d of %d two-dimensional sublattices = %.1f per cent"
          % (overlap, twod, res["overlap_pct_of_twod"]))
    print("\nwrote " + os.path.basename(OUT))


if __name__ == "__main__":
    main()
