# -*- coding: utf-8 -*-
"""Single-point electronic-structure inputs for a stratified sample of the
indexed database, so sublattice dimensionality can be set against a computed
property rather than described on its own.

No relaxation. The index is computed on the deposited coordinates, so the band
structure is computed on the same coordinates; relaxing first would mean the two
quantities describe different geometries, and would cost two orders of magnitude
more time for a comparison that does not need it.

Settings come from pymatgen's MPStaticSet, so they are community-standard rather
than ours, and POTCAR selection follows the recommended PAW potentials for each
element. POTCARs are not written per structure: one copy of each is shipped and
the job script concatenates them on the cluster, which keeps the transfer to tens
of megabytes instead of a gigabyte.

About a third of the database carries partial occupancies and cannot be written
as a POSCAR without an ordering approximation, so each stratum is filled by
walking a shuffled pool and keeping ordered structures until the quota is met.

Usage
    python setup_bandwidth_scale.py <cod_dimensionality.jsonl> <outdir>
                                    [--per-cell 30] [--seed N]
"""
import argparse
import collections
import json
import os
import random
import re
import shutil
import warnings

warnings.filterwarnings("ignore")

ELEM = re.compile(r"([A-Z][a-z]?)")
CIF_ROOT = os.path.join("C:\\", "Users", "pzgft", "CODdata", "cif", "cif")
POT_ROOT = os.path.join("E:\\", "VASP", "potentials", "potpaw_PBE.54")
CLASSES = ("oxide", "halide", "chalcogen")


def cif_path(cod, root):
    return os.path.join(root, cod[0], cod[1:3], cod[3:5], cod + ".cif")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("jsonl")
    ap.add_argument("outdir")
    ap.add_argument("--per-cell", type=int, default=30)
    ap.add_argument("--seed", type=int, default=20260911)
    ap.add_argument("--max-sites", type=int, default=40)
    ap.add_argument("--max-metals", type=int, default=16)
    ap.add_argument("--cif-root", default=CIF_ROOT)
    ap.add_argument("--pot-root", default=POT_ROOT)
    a = ap.parse_args()

    from pymatgen.core import Structure
    from pymatgen.io.vasp.sets import MPStaticSet

    cand = collections.defaultdict(list)
    for line in open(a.jsonl, encoding="utf-8"):
        r = json.loads(line)
        if r.get("D") is None:
            continue
        els = set(ELEM.findall(r["formula"]))
        if "C" in els or "H" in els:
            continue                      # organic cations add states near the gap
        if r["nsite"] > a.max_sites or r["nmetal"] > a.max_metals:
            continue
        cand[(r["anion_class"], r["d_top"])].append(r)

    rng = random.Random(a.seed)
    pools = {}
    for cls in CLASSES:
        for dt in range(4):
            pool = list(cand.get((cls, dt), []))
            rng.shuffle(pool)
            pools[(cls, dt)] = pool

    print("candidate pool per stratum:")
    print("{:10s} {:>8} {:>8} {:>8} {:>8}".format("class", "0D", "1D", "2D", "3D"))
    for cls in CLASSES:
        print("{:10s} ".format(cls) + " ".join(
            "{:>8,}".format(len(pools[(cls, dt)])) for dt in range(4)))
    print()

    os.makedirs(a.outdir, exist_ok=True)
    manifest, elements = [], set()
    skipped = collections.Counter()

    for cls in CLASSES:
        for dt in range(4):
            kept = 0
            for r in pools[(cls, dt)]:
                if kept >= a.per_cell:
                    break
                d = os.path.join(a.outdir, r["cod"])
                try:
                    st = Structure.from_file(cif_path(r["cod"], a.cif_root))
                    if not st.is_ordered:
                        skipped["disordered"] += 1
                        continue
                    st = st.get_sorted_structure()
                    vs = MPStaticSet(st, user_incar_settings={
                        "NCORE": 4, "LWAVE": False, "LCHARG": False,
                        "LORBIT": 11, "NEDOS": 3001, "ALGO": "Normal"})
                    os.makedirs(d, exist_ok=True)
                    vs.incar.write_file(os.path.join(d, "INCAR"))
                    vs.kpoints.write_file(os.path.join(d, "KPOINTS"))
                    st.to(filename=os.path.join(d, "POSCAR"), fmt="poscar")
                    syms = vs.potcar_symbols
                    elements.update(syms)
                    with open(os.path.join(d, "POTCAR.symbols"), "w",
                              newline="\n") as fh:
                        fh.write("\n".join(syms) + "\n")
                    manifest.append({
                        "cod": r["cod"], "formula": r["formula"], "cls": cls,
                        "d_top": dt, "D": r["D"], "frac": r["frac"],
                        "nsite": r["nsite"], "nmetal": r["nmetal"],
                        "natoms": len(st)})
                    kept += 1
                except Exception:
                    skipped["other"] += 1
                    if os.path.isdir(d):
                        shutil.rmtree(d, ignore_errors=True)
            print("   {:10s} {}D  kept {:>3}".format(cls, dt, kept))

    lib = os.path.join(a.outdir, "_potlib")
    os.makedirs(lib, exist_ok=True)
    missing = []
    for sym in sorted(elements):
        src = os.path.join(a.pot_root, sym, "POTCAR")
        if os.path.exists(src):
            shutil.copy(src, os.path.join(lib, sym))
        else:
            missing.append(sym)

    json.dump(manifest, open(os.path.join(a.outdir, "manifest.json"), "w"),
              indent=1)
    print("\nprepared {} structures; skipped {}".format(
        len(manifest), dict(skipped)))
    print("distinct POTCARs needed: {} (copied to _potlib)".format(len(elements)))
    if missing:
        print("MISSING from the local library: " + ", ".join(missing))
        print("  structures needing them will fail; drop those elements or add "
              "the potentials before submitting.")


if __name__ == "__main__":
    main()
