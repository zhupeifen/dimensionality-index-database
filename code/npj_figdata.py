# -*- coding: utf-8 -*-
"""Panel data for Fig. 2: the rank distribution and the fractional term.

Both panels are drawn from the whole indexed database, so everything here is a
count. Panel a is the rank composition of each bridging-anion class; panel b is
the distribution of the fractional term over the structures for which it is
non-zero.

This replaces a version produced inline, whose fractional-term histogram held
16,561 structures where the database has 16,644 with a non-zero term. The 83
missing entries all sat in the first bin, so the caption and the panel disagreed
about how many structures the figure showed. Nothing else in the file changed.

    python npj_figdata.py [cod_dimensionality.jsonl.gz]

Writes npj_figdata.json.
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
OUT = os.path.join(HERE, "npj_figdata.json")

# Panel a: the index runs 0 to 3, with headroom so the top bin is not clipped.
BINS, LO, HI = 120, 0.0, 3.2
# Panel b: the fractional term lies in [0, 1); twenty bins of 0.05, which is the
# binning fig_npj_distribution.m draws with (edges = 0:0.05:1.0).
FBINS = 20
# The three bridging classes, and the display names the figure uses for them.
CLASSES = [("chalcogen", "chalcogenide"), ("halide", "halide"), ("oxide", "oxide")]


def quantile(xs, q):
    """Linear-interpolated quantile of a sorted list."""
    if not xs:
        return None
    i = q * (len(xs) - 1)
    lo = int(i)
    hi = min(lo + 1, len(xs) - 1)
    return xs[lo] + (i - lo) * (xs[hi] - xs[lo])


def main():
    hist = {c: [0] * BINS for c, _lab in CLASSES}
    dtop = {c: collections.Counter() for c, _lab in CLASSES}
    fracs = []
    n_by_class = collections.Counter()
    n_indexed = 0

    with gzip.open(SRC, "rt", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            r = json.loads(line)
            if r.get("d_top") is None:
                continue
            n_indexed += 1
            cls = r["anion_class"]
            n_by_class[cls] += 1
            dtop[cls][r["d_top"]] += 1

            d = r.get("D")
            if d is not None:
                k = int((d - LO) / (HI - LO) * BINS)
                hist[cls][min(max(k, 0), BINS - 1)] += 1

            f = r.get("frac")
            if f:                      # non-zero; None and 0.0 both excluded
                fracs.append(f)

    fracs.sort()
    fhist = [0] * FBINS
    for f in fracs:
        k = int(f * FBINS)
        fhist[min(max(k, 0), FBINS - 1)] += 1

    res = {
        "bins": BINS, "lo": LO, "hi": HI,
        "hist": hist,
        "frac_quartiles": [round(quantile(fracs, q), 4) for q in (0.25, 0.5, 0.75)],
        "frac_n": len(fracs),
        "frac_gt_half": sum(1 for f in fracs if f > 0.5),
        "frac_hist": fhist,
        "classes": [lab for _c, lab in CLASSES] + ["all"],
        "n": [n_by_class[c] for c, _lab in CLASSES] + [n_indexed],
        "counts": [[dtop[c][i] for i in range(4)] for c, _lab in CLASSES]
                  + [[sum(dtop[c][i] for c, _l in CLASSES) for i in range(4)]],
    }
    json.dump(res, io.open(OUT, "w", encoding="utf-8", newline="\n"), indent=1)

    print("indexed %d" % n_indexed)
    for (c, lab), n in zip(CLASSES, res["n"]):
        print("  %-13s n = %-7d  %s" % (lab, n, dtop[c][0:4] if False else
                                        [dtop[c][i] for i in range(4)]))
    print("  %-13s n = %-7d  %s" % ("all", res["n"][-1], res["counts"][-1]))
    print("\nfractional term non-zero for %d structures" % res["frac_n"])
    print("  quartiles %s, above 0.5: %d"
          % (res["frac_quartiles"], res["frac_gt_half"]))
    print("  histogram sums to %d" % sum(fhist))
    print("  first bin (below 0.05): %d" % fhist[0])
    print("\nwrote " + os.path.basename(OUT))


if __name__ == "__main__":
    main()
