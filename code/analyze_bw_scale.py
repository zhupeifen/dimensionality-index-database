# -*- coding: utf-8 -*-
"""Set the sublattice dimensionality index against computed electronic structure.

For each structure in the stratified sample, read the single-point EIGENVAL and
extract three quantities: the width of the topmost valence band, the width of the
lowest conduction band, and the gap between them. Bandwidth is taken per band as
max_k E - min_k E, which needs no high-symmetry path and so is comparable across
structures of different space groups.

The question is whether a structural index computed from connectivity alone
tracks how far the band-edge states disperse. Bands are reported per anion class,
because the classes sit at different absolute bandwidths and pooling them would
confound a chemical offset with a structural trend.

Usage
    python analyze_bw_scale.py <rundir> [-o bw_scale_results.json]
"""
import argparse
import collections
import json
import math
import os


def read_eigenval(path):
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh]
    nspin = int(lines[0].split()[-1])
    _nelect, nkpt, nband = (int(float(x)) for x in lines[5].split())
    eig = [[[0.0] * nband for _ in range(nkpt)] for _ in range(nspin)]
    occ = [[[0.0] * nband for _ in range(nkpt)] for _ in range(nspin)]
    i = 6
    for k in range(nkpt):
        while i < len(lines) and not lines[i].strip():
            i += 1
        i += 1                                   # k-point coordinate line
        for b in range(nband):
            f = lines[i].split()
            for s in range(nspin):
                eig[s][k][b] = float(f[1 + s])
                occ[s][k][b] = float(f[1 + nspin + s])
            i += 1
    return nkpt, nband, nspin, eig, occ


def widths(path):
    """Valence and conduction band edges, handled per spin channel.

    Taking the highest occupied band index and then the next one is wrong for a
    spin-polarised calculation: the two channels can have different fillings, so
    band n+1 in one channel may lie far below the valence edge in the other and
    the gap comes out negative. The edges are therefore found by energy, not by
    index, and the widths are taken from the bands that actually carry them.
    """
    nkpt, nband, nspin, eig, occ = read_eigenval(path)
    vbm, cbm = -1e9, 1e9
    vb_sb, cb_sb = None, None
    for s_ in range(nspin):
        for k in range(nkpt):
            for b in range(nband):
                e, o = eig[s_][k][b], occ[s_][k][b]
                if o > 0.5:
                    if e > vbm:
                        vbm, vb_sb = e, (s_, b)
                else:
                    if e < cbm:
                        cbm, cb_sb = e, (s_, b)
    if vb_sb is None or cb_sb is None:
        return None

    def span(sb):
        s_, b = sb
        v = [eig[s_][k][b] for k in range(nkpt)]
        return max(v) - min(v)

    gap = cbm - vbm
    return {"vb_width": round(span(vb_sb), 4), "cb_width": round(span(cb_sb), 4),
            "gap": round(gap, 4), "nkpt": nkpt, "nspin": nspin,
            "metallic": gap <= 0.05}


def pearson(xs, ys):
    n = len(xs)
    if n < 3:
        return None, None
    mx = sum(xs) / n
    my = sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx <= 0 or syy <= 0:
        return None, None
    return sxy / sxx, sxy / math.sqrt(sxx * syy)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rundir")
    ap.add_argument("-o", "--out", default="bw_scale_results.json")
    a = ap.parse_args()

    man = {m["cod"]: m for m in
           json.load(open(os.path.join(a.rundir, "manifest.json"), encoding="utf-8"))}
    rows, missing = [], 0
    for cod, m in man.items():
        p = os.path.join(a.rundir, cod, "EIGENVAL")
        if not os.path.exists(p) or os.path.getsize(p) < 1000:
            missing += 1
            continue
        try:
            w = widths(p)
        except Exception:
            w = None
        if not w:
            missing += 1
            continue
        w.update(m)
        rows.append(w)

    print("structures with usable EIGENVAL: {} of {}  ({} pending or failed)".format(
        len(rows), len(man), missing))
    if not rows:
        return

    print("\n{:10s} {:>5} {:>5} {:>9} {:>9} {:>8}".format(
        "class", "rank", "n", "VB width", "CB width", "gap"))
    for cls in ("oxide", "halide", "chalcogen"):
        for dt in range(4):
            g = [r for r in rows if r["cls"] == cls and r["d_top"] == dt]
            if not g:
                continue
            print("{:10s} {:>5} {:>5} {:>9.3f} {:>9.3f} {:>8.3f}".format(
                cls, "{}D".format(dt), len(g),
                sum(r["vb_width"] for r in g) / len(g),
                sum(r["cb_width"] for r in g) / len(g),
                sum(r["gap"] for r in g) / len(g)))

    print("\nvalence bandwidth against the index, within each anion class:")
    corr = {}
    for cls in ("oxide", "halide", "chalcogen"):
        g = [r for r in rows if r["cls"] == cls]
        if len(g) < 6:
            continue
        sl, r_ = pearson([x["D"] for x in g], [x["vb_width"] for x in g])
        if sl is None:
            continue
        corr[cls] = {"n": len(g), "slope": round(sl, 4), "r": round(r_, 3)}
        print("   {:10s} n={:>4}  slope {:+.3f} eV per unit D   r = {:+.3f}".format(
            cls, len(g), sl, r_))
    g = rows
    sl, r_ = pearson([x["D"] for x in g], [x["vb_width"] for x in g])
    if sl is not None:
        corr["all"] = {"n": len(g), "slope": round(sl, 4), "r": round(r_, 3)}
        print("   {:10s} n={:>4}  slope {:+.3f} eV per unit D   r = {:+.3f}".format(
            "pooled", len(g), sl, r_))

    json.dump({"n": len(rows), "correlations": corr, "rows": rows},
              open(a.out, "w"), indent=1)
    print("\nwrote " + a.out)


if __name__ == "__main__":
    main()
