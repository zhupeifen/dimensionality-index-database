# -*- coding: utf-8 -*-
"""What does excluding the two-anion structures cost?

The paper indexes a structure only where one bridging anion class is present,
and sets aside the 148,440 depositions carrying two — an oxide and a halide,
say — on the grounds that the bridging ligand is not defined without a further
choice. That is the largest single judgment in the paper: it removes more
structures than it keeps. It has not been tested.

This samples those structures, indexes each on every anion class it contains,
and asks how often the choice changes the integer rank. If the answer is rarely,
the exclusion is conservative and could be relaxed. If it is often, the
exclusion is necessary and the paper can say by how much.

    python anion_choice.py [n_sample] [cif_root]

Writes anion_choice.json.
"""
import collections
import gzip
import importlib.util
import io
import json
import os
import random
import sys
import warnings

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))

def data(name):
    """Find a data file whether run from the working tree or from the deposit."""
    for p in (os.path.join(HERE, "deposit_npj", "data", name),
              os.path.join(HERE, os.pardir, "data", name),
              os.path.join(HERE, name)):
        if os.path.exists(p):
            return p
    return os.path.join(HERE, name)


def code(name):
    """Find a code file whether run from the working tree or from the deposit."""
    for p in (os.path.join(HERE, "deposit_npj", "code", name),
              os.path.join(HERE, name)):
        if os.path.exists(p):
            return p
    return os.path.join(HERE, name)

N = int(sys.argv[1]) if len(sys.argv) > 1 else 1200
ROOT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
    "C:\\", "Users", "pzgft", "CODdata", "cif", "cif")
INDEXED = data("cod_dimensionality.jsonl.gz")
OUT = os.path.join(HERE, "anion_choice.json")
SEED = 20260914


def load_indexer():
    path = code("cod_wide_index.py")
    spec = importlib.util.spec_from_file_location("cwi", path)
    mod = importlib.util.module_from_spec(spec)
    saved = sys.argv
    sys.argv = ["cod_wide_index.py"]
    try:
        spec.loader.exec_module(mod)
    finally:
        sys.argv = saved
    return mod


def cif_path(cod):
    p = os.path.join(ROOT, cod[0], cod[1:3], cod[3:5], cod + ".cif")
    return p if os.path.exists(p) else None


def main():
    from pymatgen.core import Structure
    cwi = load_indexer()

    pool = []
    with gzip.open(INDEXED, "rt", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            r = json.loads(line)
            if r.get("skip") == "anion class not defined":
                pool.append(r["cod"])
    print("structures set aside for a competing anion class: %d" % len(pool))
    random.seed(SEED)
    random.shuffle(pool)
    sample = pool[:N]

    # index_structure picks one class per structure. Reach past it: for each
    # class present, restrict ANION_CUT to that class alone and index.
    base = {c: dict(v) for c, v in cwi.ANION_CUT.items()}
    rows, done = [], 0
    for cod in sample:
        p = cif_path(cod)
        if p is None:
            continue
        try:
            st = Structure.from_file(p)
        except Exception:
            continue
        els = {s.specie.symbol for s in st.sites
               if hasattr(s, "specie")} if st else set()
        present = [c for c, cuts in base.items() if els & set(cuts)]
        if len(present) < 2:
            continue
        got, bridges = {}, {}
        for c in present:
            cwi.ANION_CUT.clear()
            cwi.ANION_CUT[c] = dict(base[c])
            try:
                res, _ = cwi.index_structure(st)
            except Exception:
                res = None
            if res:
                got[c] = res["d_top"]
                bridges[c] = res["nbridge"]
        cwi.ANION_CUT.clear()
        cwi.ANION_CUT.update(base)
        if len(got) >= 2:
            rows.append({"cod": cod, "ranks": got, "nbridge": bridges})
        done += 1
        if done % 200 == 0:
            print("   %d/%d" % (done, len(sample)), flush=True)

    agree_rows = [r for r in rows if len(set(r["ranks"].values())) == 1]
    agree = len(agree_rows)
    spread = collections.Counter(
        max(r["ranks"].values()) - min(r["ranks"].values()) for r in rows)
    pairs = collections.Counter(
        tuple(sorted(r["ranks"])) for r in rows)

    # Agreement is cheap where nothing bridges on any candidate class: the metal
    # sublattice is isolated whichever anion is chosen, and both answers are
    # rank zero. Separating that case from a real agreement is the point of the
    # block below - the exclusion only costs something where a bridge exists.
    at_zero = sum(1 for r in agree_rows
                  if next(iter(r["ranks"].values())) == 0)
    bridging = [r for r in rows if max(r["nbridge"].values()) > 0]
    choice_dep = [r for r in bridging if len(set(r["ranks"].values())) > 1]

    def pct(a, b):
        return round(100.0 * a / b, 1) if b else None

    res = {"pool": len(pool), "sampled": len(sample), "evaluated": len(rows),
           "agree": agree,
           "agree_pct": pct(agree, len(rows)),
           "differ_pct": pct(len(rows) - agree, len(rows)),
           "spread": {str(k): v for k, v in sorted(spread.items())},
           "class_pairs": {"+".join(k): v for k, v in pairs.most_common()},
           "bridging": {
               "evaluated": len(rows),
               "agreements_at_rank_zero": at_zero,
               "agreements_at_rank_zero_pct": pct(at_zero, agree),
               "bridge_on_some_class": len(bridging),
               "bridge_pct": pct(len(bridging), len(rows)),
               "rank_choice_dependent": len(choice_dep),
               "rank_choice_dependent_pct": pct(len(choice_dep), len(bridging)),
               "note": ("every figure in this block shares the one denominator "
                        "above: the structures indexed on two or more classes"),
           }}
    json.dump(res, io.open(OUT, "w", encoding="utf-8", newline="\n"), indent=1)

    print("\nevaluated on two or more classes: %d" % len(rows))
    print("  same integer rank whichever class is chosen : %d = %s per cent"
          % (agree, res["agree_pct"]))
    print("  rank depends on the choice                  : %d = %s per cent"
          % (len(rows) - agree, res["differ_pct"]))
    print("\n  what the agreement is worth:")
    print("     agreements that are agreements at rank 0 : %d of %d = %s per cent"
          % (at_zero, agree, res["bridging"]["agreements_at_rank_zero_pct"]))
    print("     bridge on at least one candidate class   : %d of %d = %s per cent"
          % (len(bridging), len(rows), res["bridging"]["bridge_pct"]))
    print("     of those, rank depends on the choice     : %d of %d = %s per cent"
          % (len(choice_dep), len(bridging),
             res["bridging"]["rank_choice_dependent_pct"]))
    print("\n  spread in rank across the classes present:")
    for k, v in sorted(spread.items()):
        print("     %d rank(s)  %5d  %5.1f%%" % (k, v, 100.0 * v / len(rows)))
    print("\n  commonest class combinations:")
    for k, v in pairs.most_common(5):
        print("     %-22s %5d" % ("+".join(k), v))
    print("\nwrote " + os.path.basename(OUT))


if __name__ == "__main__":
    main()
