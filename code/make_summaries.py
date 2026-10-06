# -*- coding: utf-8 -*-
"""Rebuild the summary files behind the paper from the deposited per-structure data.

The scripts that first wrote these summaries were lost with the machine they ran
on. This reconstruction was checked against the files the paper was built from:
npj_stats.json, npj_bymetal.json, npj_wholenet.json and the fit block of
npj_bw.json reproduce exactly (only the stored p-values, which the paper does not
use, differ in the last digits). It also writes the permissive whole-network
reading, which the paper quotes but no file held.

The classification of disagreeing compositions (npj_disagreements.json,
npj_fig2data.json; Figure 5) uses the rule stated here, which is the rule the
paper describes. For each composition, the most similar pair of depositions at
different ranks is found (smallest relative difference in cell volume per atom,
ties broken towards a shared space group):
  more than 3 per cent apart              -> different density
  within LO and the same space group      -> index instability
  within LO and different space groups    -> polymorphism
  between LO and 3 per cent               -> other (straddles the tolerances)
A composition with no pair of depositions at different ranks in its records is
"other" (none remains once complete_disagreement_records.py has run). Because
each composition is classed once, by its closest pair, a composition can hold a
same-phase pair at different ranks yet be classed otherwise; those are listed as
"also_same_phase", so the index-instability count is a floor. LO is set by --lo
(default 0.01).

    python make_summaries.py <data_dir> <out_dir> [--lo 0.01]

<data_dir> holds cod_dimensionality.jsonl.gz, wholenet_scale.json,
bw_scale_results.json and npj_disagreements.json (whose "detail" records give,
per deposition, [cod, rank, space group number, volume per atom, sites]).
"""
import argparse
import gzip
import itertools
import json
import os
import re
import statistics
from collections import Counter, defaultdict

import numpy as np
from scipy import stats

SPECTATOR = {"Li", "Na", "K", "Rb", "Cs", "Fr", "Be", "Mg", "Ca", "Sr", "Ba", "Ra",
             "H", "C", "N", "P", "B", "Si"}
ANIONS = {"Cl", "Br", "I", "F", "S", "Se", "Te", "O"}
CLASSES = ("oxide", "halide", "chalcogen")


def dump(obj, path):
    json.dump(obj, open(path, "w", encoding="utf-8"), indent=1)
    print("wrote", path)


def load_db(path):
    recs, skips = [], Counter()
    with gzip.open(path, "rt", encoding="utf-8") as f:
        for line in f:
            r = json.loads(line)
            if "skip" in r:
                skips[r["skip"]] += 1
            else:
                recs.append(r)
    return recs, skips


def ranks(rr):
    return {str(k): v for k, v in sorted(Counter(r["d_top"] for r in rr).items())}


def stats_block(recs, skips):
    byf = defaultdict(list)
    for r in recs:
        byf[r["formula"]].append(r)
    rep = [f for f, v in byf.items() if len(v) > 1]
    dis = [f for f in rep if len({r["d_top"] for r in byf[f]}) > 1]

    def same_nsite(f):
        v = byf[f]
        return any(a["nsite"] == b["nsite"] and a["d_top"] != b["d_top"]
                   for a, b in itertools.combinations(v, 2))
    return {"total_cifs": len(recs) + sum(skips.values()), "indexed": len(recs),
            "skips": dict(skips), "dtop": ranks(recs),
            "by_class": {c: {"n": sum(1 for r in recs if r["anion_class"] == c),
                             "dtop": ranks([r for r in recs if r["anion_class"] == c])}
                         for c in CLASSES},
            "frac_nonzero": sum(1 for r in recs if r["frac"] > 0),
            "nmetal_median": statistics.median(r["nmetal"] for r in recs),
            "repeat_compositions": len(rep), "repeat_disagree": len(dis),
            "disagree_same_nsite": sum(map(same_nsite, dis)),
            "frac_undefined": sum(1 for r in recs if r["nbridge"] == 0),
            "frac_defined": sum(1 for r in recs if r["nbridge"] > 0)}


def bymetal_block(recs, threshold=1900):
    bym = defaultdict(Counter)
    for r in recs:
        m = set(re.findall(r"[A-Z][a-z]?", r["formula"])) - SPECTATOR - ANIONS
        if len(m) == 1:
            bym[m.pop()][r["d_top"]] += 1
    out = {}
    for m, c in sorted(bym.items(), key=lambda kv: -sum(kv[1].values())):
        n = sum(c.values())
        if n > threshold:
            out[m] = {"n": n, "pct": [round(100 * c[k] / n, 1) for k in range(4)],
                      "mean": round(sum(k * c[k] for k in range(4)) / n, 2)}
    return out


def wholenet_block(w):
    rows = w["rows"]
    ev = [r for r in rows if r.get("rda")]

    def strict(r):
        return r["rda"] != "{}D".format(r["d_top"])

    def permissive(r):
        return str(r["d_top"]) not in r["rda"][:-1]

    def pct(rr, f):
        return round(100 * sum(map(f, rr)) / len(rr), 1)
    o = {"evaluated": len(ev), "sampled": len(rows), "differ_pct": pct(ev, strict),
         "by_class": {c: {"n": len(x), "differ_pct": pct(x, strict)}
                      for c in CLASSES for x in [[r for r in ev if r["cls"] == c]]},
         "by_rank": {str(k): pct([r for r in ev if r["d_top"] == k], strict) for k in range(4)},
         "mixed_pct": round(100 * sum(1 for r in ev if len(r["rda"]) > 2) / len(ev), 1),
         "skipped_large": len(rows) - len(ev)}
    o["halide_random_pct"] = o["by_class"]["halide"]["differ_pct"]
    o["permissive"] = {"differ_pct": pct(ev, permissive),
                       "by_class": {c: pct([r for r in ev if r["cls"] == c], permissive)
                                    for c in CLASSES},
                       "by_rank": {str(k): pct([r for r in ev if r["d_top"] == k], permissive)
                                   for k in range(4)}}
    return o


def bw_block(rows):
    rows = [r for r in rows if r.get("vb_width") is not None]

    def fit(rr):
        x = np.array([r["D"] for r in rr])
        y = np.array([r["vb_width"] for r in rr])
        L = stats.linregress(x, y)
        return {"n": len(rr), "slope": round(L.slope, 3), "r": round(L.rvalue, 3),
                "p": L.pvalue, "intercept": round(L.intercept, 3)}

    def r2(cols, y):
        X = np.column_stack([np.ones(len(y))] + cols)
        b, *_ = np.linalg.lstsq(X, y, rcond=None)
        res = y - X @ b
        return 1 - res @ res / ((y - y.mean()) @ (y - y.mean()))

    def base(rr):
        y = np.array([r["vb_width"] for r in rr])
        D = np.array([r["D"] for r in rr])
        de = np.array([r["density"] for r in rr])
        a, b, c = r2([de], y), r2([D], y), r2([D, de], y)
        df = len(y) - 3
        return {"density_r2": round(a, 3), "D_r2": round(b, 3), "both_r2": round(c, 3),
                "dR2": round(c - a, 3), "F": round((c - a) / ((1 - c) / df), 1), "df": df}
    o = {"n_total": len(rows), "n_metallic": sum(1 for r in rows if r["metallic"]),
         "pooled": fit(rows), "insulating": fit([r for r in rows if not r["metallic"]]),
         "by_class": {c: fit([r for r in rows if r["cls"] == c]) for c in CLASSES},
         "means": {c: {str(k): round(float(np.mean([r["vb_width"] for r in rows
                                                     if r["cls"] == c and r["d_top"] == k])), 3)
                       for k in range(4)} for c in CLASSES},
         "counts": {c: {str(k): sum(1 for r in rows if r["cls"] == c and r["d_top"] == k)
                        for k in range(4)} for c in CLASSES},
         "points": [{"D": r["D"], "w": round(r["vb_width"], 4), "cls": r["cls"],
                     "met": bool(r["metallic"])} for r in rows]}
    o["baseline"] = {"pooled": base(rows), **{c: base([r for r in rows if r["cls"] == c])
                                              for c in ("chalcogen", "halide", "oxide")}}
    return o


def classify(recs, lo):
    def dv(a, b):
        return abs(a[3] - b[3]) / min(a[3], b[3])
    pairs = [(a, b) for a, b in itertools.combinations(recs, 2) if a[1] != b[1]]
    if not pairs:
        return "other"
    a, b = min(pairs, key=lambda p: (dv(*p), p[0][2] != p[1][2]))
    d = dv(a, b)
    if d > 0.03:
        return "different density"
    if d <= lo:
        return "index instability" if a[2] == b[2] else "polymorphism"
    return "other"


def dv(a, b):
    return abs(a[3] - b[3]) / min(a[3], b[3])


def disagreement_blocks(detail, recs, lo):
    frac = {r["cod"]: r["frac"] for r in recs}
    order = ["different density", "polymorphism", "index instability", "other"]
    cat = {}
    for x in detail:
        cat[x["formula"]] = classify(x["recs"], lo)
    counts = Counter(cat.values())
    n_struct, above = Counter(), Counter()
    for x in detail:
        k = cat[x["formula"]]
        n_struct[k] += len(x["recs"])
        above[k] += sum(1 for r in x["recs"] if frac.get(r[0], 0) > 0.1)
    also = sorted(x["formula"] for x in detail if cat[x["formula"]] != "index instability"
                  and any(a[1] != b[1] and a[2] == b[2] and dv(a, b) <= lo
                          for a, b in itertools.combinations(x["recs"], 2)))
    dis = {"n": len(detail), "rule": {"density_tolerance": 0.03, "matched_density": lo},
           "categories": {k: counts[k] for k in order}, "also_same_phase": also,
           "detail": [dict(x, category=cat[x["formula"]]) for x in detail]}
    fig = {"order": order, "counts": [counts[k] for k in order],
           "frac_above": [round(100 * above[k] / n_struct[k], 1) if n_struct[k] else 0.0
                          for k in order],
           "n_struct": [n_struct[k] for k in order],
           "baseline_above": round(100 * sum(1 for v in frac.values() if v > 0.1) / len(frac), 1),
           "total": len(detail), "also_same_phase": len(also)}
    return dis, fig


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("data")
    ap.add_argument("out")
    ap.add_argument("--lo", type=float, default=0.01)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    recs, skips = load_db(os.path.join(a.data, "cod_dimensionality.jsonl.gz"))
    dump(stats_block(recs, skips), os.path.join(a.out, "npj_stats.json"))
    dump(bymetal_block(recs), os.path.join(a.out, "npj_bymetal.json"))
    dump(wholenet_block(json.load(open(os.path.join(a.data, "wholenet_scale.json"),
                                       encoding="utf-8"))),
         os.path.join(a.out, "npj_wholenet.json"))
    bw = json.load(open(os.path.join(a.data, "bw_scale_results.json"), encoding="utf-8"))
    dump(bw_block(bw if isinstance(bw, list) else list(bw.values())),
         os.path.join(a.out, "npj_bw_fits.json"))
    detail = json.load(open(os.path.join(a.data, "npj_disagreements.json"),
                            encoding="utf-8"))["detail"]
    dis, fig = disagreement_blocks(detail, recs, a.lo)
    dump(dis, os.path.join(a.out, "npj_disagreements.json"))
    dump(fig, os.path.join(a.out, "npj_fig2data.json"))


if __name__ == "__main__":
    main()
