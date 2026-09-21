# -*- coding: utf-8 -*-
"""Is the database-scale distribution an artefact of the bridging cutoff?

The paper reports a distribution of integer periodic rank over the indexed
database. That distribution rests on one bridging criterion, taken from the gap
between the first coordination shell and the next group of contacts. The
companion article shows the integer topology is stable for a handful of
structures and for 703 halides; nothing shows it for the database this paper
reports on, and the sensitivity of a distribution is not the sensitivity of a
single assignment.

This re-indexes a stratified sample at the criterion used and at scaled
alternatives, and counts how many structures change integer rank.

    python cutoff_sensitivity.py [n_per_stratum] [cif_root]

Writes cutoff_sensitivity.json.
"""
import collections
import gzip
import importlib.util
import json
import io
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

N_PER = int(sys.argv[1]) if len(sys.argv) > 1 else 400
ROOT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
    "C:\\", "Users", "pzgft", "CODdata", "cif", "cif")
INDEXED = data("cod_dimensionality.jsonl.gz")
OUT = os.path.join(HERE, "cutoff_sensitivity.json")

# Scale every anion's cutoff by the same factor, so the comparison is of one
# criterion loosened or tightened rather than of twelve unrelated numbers.
SCALES = [0.95, 1.00, 1.05]


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
    base_cut = {c: dict(v) for c, v in cwi.ANION_CUT.items()}

    # a stratified sample: N_PER from each anion class at each integer rank, so
    # the rare ranks are represented as well as the dominant one
    pools = collections.defaultdict(list)
    with gzip.open(INDEXED, "rt", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            r = json.loads(line)
            if r.get("d_top") is None:
                continue
            pools[(r["anion_class"], r["d_top"])].append(r["cod"])
    random.seed(20260913)
    sample = []
    pop = {k: len(v) for k, v in pools.items()}          # true stratum sizes
    for key, codes in sorted(pools.items()):
        random.shuffle(codes)
        sample.extend((c, key) for c in codes[:N_PER])
    print("strata %d, sampled %d structures of %d indexed"
          % (len(pools), len(sample), sum(pop.values())))

    # Parse once and index at every scale while the structure is in hand. Holding
    # thousands of parsed structures to make three passes over them costs far more
    # memory than it saves time, and parsing is the expensive step either way.
    ranks = {s: {} for s in SCALES}
    parsed = []
    scaled = {s: {cls: {el: d * s for el, d in cuts.items()}
                  for cls, cuts in base_cut.items()} for s in SCALES}
    done = 0
    for cod, key in sample:
        p = cif_path(cod)
        if p is None:
            continue
        try:
            st = Structure.from_file(p)
        except Exception:
            continue
        parsed.append((cod, key))
        for s in SCALES:
            cwi.ANION_CUT.update(scaled[s])
            try:
                res, _ = cwi.index_structure(st)
            except Exception:
                continue
            if res:
                ranks[s][cod] = res["d_top"]
        done += 1
        if done % 200 == 0:
            print("   %d/%d" % (done, len(sample)), flush=True)
    print("parsed and indexed %d of %d" % (len(parsed), len(sample)))

    base = ranks[1.00]
    res = {"n_sampled": len(parsed), "n_per_stratum": N_PER, "scales": SCALES,
           "base_cutoffs": base_cut, "changed": {}, "distribution": {}}
    for s in SCALES:
        common = [c for c in base if c in ranks[s]]
        moved = sum(1 for c in common if ranks[s][c] != base[c])
        res["changed"]["%.2f" % s] = {"n": len(common), "moved": moved,
                                      "pct": round(100.0 * moved / len(common), 2)}
        # The sample deliberately over-represents the rare ranks, so a raw count
        # of it is not the database. Weight each structure by the size of the
        # stratum it was drawn from, divided by how many of that stratum were
        # actually indexed here.
        drawn = collections.Counter(key for cod, key in parsed if cod in ranks[s])
        wsum = collections.Counter()
        for cod, key in parsed:
            if cod not in ranks[s] or cod not in base:
                continue
            wsum[ranks[s][cod]] += pop[key] / drawn[key]
        tot = sum(wsum.values())
        res["distribution"]["%.2f" % s] = {str(k): round(100.0 * wsum[k] / tot, 1)
                                           for k in (0, 1, 2, 3)}
    json.dump(res, io.open(OUT, "w", encoding="utf-8", newline="\n"), indent=1)

    print("\n scale   n      rank changed        distribution 0/1/2/3 (%)")
    for s in SCALES:
        k = "%.2f" % s
        c, d = res["changed"][k], res["distribution"][k]
        print("  %.2f  %5d   %5d = %5.2f%%     %s" %
              (s, c["n"], c["moved"], c["pct"],
               "  ".join(d[x] and "%4.1f" % d[x] or " 0.0" for x in ("0", "1", "2", "3"))))
    print("\nwrote " + os.path.basename(OUT))


if __name__ == "__main__":
    main()
