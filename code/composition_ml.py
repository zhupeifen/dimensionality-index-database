# -*- coding: utf-8 -*-
"""Can a real model recover connectivity from composition, where a lookup cannot?

The paper's composition baseline is a majority label per (anion class, X/M ratio)
bin. It never returns two dimensions for any composition, and the paper reads
that as evidence that composition cannot see layers. That reading has a hole in
it: two-dimensional sublattices are 4.2 per cent of the database, so a majority
rule cannot select them almost by construction, and the zero recall may be an
artefact of majority labelling rather than a fact about chemistry.

This closes the hole. A gradient-boosted classifier is fitted on element
fractions and the anion-to-metal ratio, with balanced class weights so that the
rare ranks are not simply ignored, and scored out of sample by stratified
five-fold cross-validation. If two-dimensional recall stays near zero under a
model that is both stronger and explicitly told to care about rare classes, then
the claim is about chemistry and not about the estimator.

Balanced accuracy is reported alongside plain accuracy, because a rule that
always answers "zero-dimensional" scores 78 per cent on the second and 25 on the
first.

    python composition_ml.py [cod_dimensionality.jsonl.gz] [--fast]

Writes composition_ml.json.
"""
import collections
import gzip
import io
import json
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def data(name):
    for p in (os.path.join(HERE, "deposit_npj", "data", name),
              os.path.join(HERE, os.pardir, "data", name),
              os.path.join(HERE, name)):
        if os.path.exists(p):
            return p
    return os.path.join(HERE, name)


ARGS = [a for a in sys.argv[1:] if not a.startswith("--")]
FAST = "--fast" in sys.argv
SRC = ARGS[0] if ARGS else data("cod_dimensionality.jsonl.gz")
OUT = os.path.join(HERE, "composition_ml.json")
SEED = 20260914
FOLDS = 5

NONMETAL = {"H", "D", "C", "N", "O", "S", "Se", "Te", "P", "F",
            "Cl", "Br", "I", "B", "Si", "As"}
ANION = {"oxide": {"O"}, "halide": {"F", "Cl", "Br", "I"},
         "chalcogen": {"S", "Se", "Te"}}
CLASSES = ["oxide", "halide", "chalcogen"]


def parse(formula):
    """Element -> count, straight from the deposited formula string."""
    out = collections.Counter()
    for sym, num in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        if sym:
            out[sym] += int(num) if num else 1
    return out


def ratio(counts, cls):
    anion = sum(n for s, n in counts.items() if s in ANION[cls])
    metal = sum(n for s, n in counts.items()
                if s not in NONMETAL and s not in ANION[cls])
    return (anion / metal) if metal else None


def main():
    rows, y, groups = [], [], []
    with gzip.open(SRC, "rt", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            r = json.loads(line)
            if r.get("d_top") is None:
                continue
            c = parse(r["formula"])
            q = ratio(c, r["anion_class"])
            if q is None:
                continue                      # same exclusion as the lookup
            rows.append((c, q, r["anion_class"]))
            y.append(r["d_top"])
            groups.append(r["formula"])
    y = np.asarray(y)
    groups = np.asarray(groups)
    print("structures with a usable X/M ratio: %d" % len(y))
    print("rank distribution: %s" % dict(collections.Counter(y.tolist())))

    # Features: the fraction of each element present, the anion-to-metal ratio,
    # and the anion class. Element fractions are the standard composition-only
    # representation; nothing here knows anything about geometry.
    elements = sorted({s for c, _q, _cl in rows for s in c})
    idx = {s: i for i, s in enumerate(elements)}
    X = np.zeros((len(rows), len(elements) + 1 + len(CLASSES)), dtype=np.float32)
    for i, (c, q, cl) in enumerate(rows):
        tot = float(sum(c.values())) or 1.0
        for s, n in c.items():
            X[i, idx[s]] = n / tot
        X[i, len(elements)] = min(q, 12.0)
        X[i, len(elements) + 1 + CLASSES.index(cl)] = 1.0
    print("features: %d element fractions + ratio + %d class flags = %d"
          % (len(elements), len(CLASSES), X.shape[1]))

    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.model_selection import StratifiedKFold, GroupKFold
    from sklearn.metrics import balanced_accuracy_score
    from sklearn.utils.class_weight import compute_sample_weight

    print("distinct compositions: %d among %d structures"
          % (len(set(groups.tolist())), len(y)))

    def run(splitter, name, **kw):
        pr = np.zeros_like(y)
        for k, (tr, te) in enumerate(splitter.split(X, y, **kw), 1):
            clf = HistGradientBoostingClassifier(
                max_iter=120 if FAST else 300, learning_rate=0.1,
                max_leaf_nodes=31, early_stopping=False, random_state=SEED)
            # balanced weights: the model is told the rare ranks matter
            w = compute_sample_weight("balanced", y[tr])
            clf.fit(X[tr], y[tr], sample_weight=w)
            pr[te] = clf.predict(X[te])
            print("   %-7s fold %d/%d" % (name, k, FOLDS), flush=True)
        return pr

    # A composition is deposited many times over, so a random split puts the
    # same formula in train and test and lets the model memorise it rather than
    # generalise. The grouped split holds every deposition of a composition out
    # together, and is the number that means anything.
    pred_rand = run(StratifiedKFold(n_splits=FOLDS, shuffle=True,
                                    random_state=SEED), "random")
    pred = run(GroupKFold(n_splits=FOLDS), "grouped", groups=groups)

    acc = float((pred == y).mean())
    bal = float(balanced_accuracy_score(y, pred))
    per = {}
    for k in (0, 1, 2, 3):
        m = y == k
        per[str(k)] = {"n": int(m.sum()), "recalled": int((pred[m] == k).sum()),
                       "recall_pct": round(100.0 * float((pred[m] == k).mean()), 1)}
    ext = y != 0
    commonest = float((y == collections.Counter(y.tolist()).most_common(1)[0][0]).mean())

    per_rand = {str(k): round(100.0 * float((pred_rand[y == k] == k).mean()), 1)
                for k in (0, 1, 2, 3)}

    LU = json.load(open(data("composition_baseline.json"), encoding="utf-8"))
    res = {
        "n": int(len(y)), "folds": FOLDS, "seed": SEED,
        "n_compositions": int(len(set(groups.tolist()))),
        "n_features": int(X.shape[1]), "n_elements": len(elements),
        "model": "HistGradientBoostingClassifier, balanced sample weights",
        "scoring": "%d-fold cross-validation grouped by composition, out of sample" % FOLDS,
        "random_split_accuracy_pct": round(100.0 * float((pred_rand == y).mean()), 1),
        "random_split_balanced_pct": round(100.0 * balanced_accuracy_score(y, pred_rand), 1),
        "random_split_per_rank_recall": per_rand,
        "accuracy_pct": round(100.0 * acc, 1),
        "balanced_accuracy_pct": round(100.0 * bal, 1),
        "commonest_rank_pct": round(100.0 * commonest, 1),
        "per_rank": per,
        "extended_recall_pct": round(100.0 * float((pred[ext] == y[ext]).mean()), 1),
        "ranks_ever_predicted": sorted(int(v) for v in set(pred.tolist())),
        "lookup_accuracy_pct": LU["lookup_pct"],
        "lookup_per_rank": {k: v["pct"] for k, v in LU["per_rank"].items()},
    }
    json.dump(res, io.open(OUT, "w", encoding="utf-8", newline="\n"), indent=1)

    print("\n  RANDOM split, the same composition allowed in train and test:")
    print("      accuracy %.1f   balanced %.1f   2D recall %.1f   <- leaky, for contrast"
          % (res["random_split_accuracy_pct"], res["random_split_balanced_pct"],
             per_rand["2"]))
    print("\n  GROUPED split, every deposition of a composition held out together:")
    print("\n                      accuracy   balanced   2D recall")
    print("  always commonest    %6.1f     %6.1f     %6.1f"
          % (res["commonest_rank_pct"], 25.0, 0.0))
    print("  ratio lookup        %6.1f          -     %6.1f"
          % (LU["lookup_pct"], LU["per_rank"]["2"]["pct"]))
    print("  gradient boosting   %6.1f     %6.1f     %6.1f"
          % (res["accuracy_pct"], res["balanced_accuracy_pct"],
             per["2"]["recall_pct"]))
    print("\n  per-rank recall, out of sample:")
    for k in ("0", "1", "2", "3"):
        print("     %sD  %7d of %7d = %5.1f per cent  (lookup %5.1f)"
              % (k, per[k]["recalled"], per[k]["n"], per[k]["recall_pct"],
                 LU["per_rank"][k]["pct"]))
    print("\n  ranks the model ever predicts:", res["ranks_ever_predicted"])
    print("\nwrote " + os.path.basename(OUT))


if __name__ == "__main__":
    main()
