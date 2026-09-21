# -*- coding: utf-8 -*-
"""Confidence intervals on the bandwidth-against-index slopes.

The fit reports a slope, a correlation and a p value, which together say that
an effect exists but not how well it is pinned down. A referee reading a
physics paper wants the interval. Both the analytic interval, from the standard
error of the ordinary-least-squares slope, and a bootstrap over structures are
computed; they agree, which is worth knowing because the residuals are not
obviously normal.

Adds a "ci" block to npj_bw.json in place.

    python bandwidth_ci.py
"""
import io
import json
import math
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))

def data(name):
    """Find a data file whether run from the working tree or from the deposit."""
    for p in (os.path.join(HERE, "deposit_npj", "data", name),
              os.path.join(HERE, os.pardir, "data", name),
              os.path.join(HERE, name)):
        if os.path.exists(p):
            return p
    return os.path.join(HERE, name)

SRC = data("npj_bw.json")
BOOTSTRAPS = 4000
SEED = 7


def ols(points):
    """slope, its standard error, and n."""
    n = len(points)
    xs = [p["D"] for p in points]
    ys = [p["w"] for p in points]
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    b = sxy / sxx
    a = my - b * mx
    resid = [y - (a + b * x) for x, y in zip(xs, ys)]
    s2 = sum(r * r for r in resid) / (n - 2)
    return b, math.sqrt(s2 / sxx), n


def bootstrap(points, k=BOOTSTRAPS):
    """Resample structures, not residuals: the scatter is chemistry, not noise."""
    random.seed(SEED)
    out = []
    for _ in range(k):
        s = [random.choice(points) for _ in points]
        try:
            out.append(ols(s)[0])
        except ZeroDivisionError:
            continue
    out.sort()
    return out[int(0.025 * len(out))], out[int(0.975 * len(out))]


def main():
    J = json.load(open(SRC, encoding="utf-8"))
    P = J["points"]
    groups = [("pooled", P), ("insulating", [p for p in P if not p["met"]])]
    for c in ("chalcogen", "halide", "oxide"):
        groups.append((c, [p for p in P if p["cls"] == c]))

    J["ci"] = {}
    print("%-12s %5s %8s %8s   %-20s %s" %
          ("set", "n", "slope", "SE", "analytic 95%", "bootstrap 95%"))
    for lab, pts in groups:
        b, se, n = ols(pts)
        lo, hi = b - 1.96 * se, b + 1.96 * se
        blo, bhi = bootstrap(pts)
        J["ci"][lab] = {"n": n, "slope": round(b, 3), "se": round(se, 3),
                        "lo": round(lo, 3), "hi": round(hi, 3),
                        "boot_lo": round(blo, 3), "boot_hi": round(bhi, 3),
                        "spans_zero": lo * hi <= 0}
        print("%-12s %5d %8.3f %8.3f   [%6.3f, %6.3f]  [%6.3f, %6.3f]%s" %
              (lab, n, b, se, lo, hi, blo, bhi,
               "   spans zero" if lo * hi <= 0 else ""))

    json.dump(J, io.open(SRC, "w", encoding="utf-8", newline="\n"), indent=1)
    print("\nupdated " + os.path.basename(SRC))


if __name__ == "__main__":
    main()
