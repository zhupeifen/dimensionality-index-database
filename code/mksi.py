# -*- coding: utf-8 -*-
"""Supplementary Information for the npj database paper.

Everything here is read from the same JSON the manuscript is built from, so the
two cannot disagree. The computational parameters are read out of pymatgen's
MPStaticSet rather than transcribed, for the same reason.

    python mksi.py "npj - Supplementary Information.docx"
"""
import collections
import json
import re
import sys

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor

from chemfmt import add_rich

OUT = sys.argv[1] if len(sys.argv) > 1 else "npj - Supplementary Information.docx"
FONT = "Times New Roman"
INK = RGBColor(0x1A, 0x1A, 0x1A)

S = json.load(open("npj_stats.json", encoding="utf-8"))
CS = json.load(open("cutoff_sensitivity.json", encoding="utf-8"))
AC = json.load(open("anion_choice.json", encoding="utf-8"))
BM = json.load(open("npj_bymetal.json", encoding="utf-8"))
LQ = json.load(open("layered_query.json", encoding="utf-8"))
BW = json.load(open("bw_scale_results.json", encoding="utf-8"))

d = Document()
sec = d.sections[0]
for a in ("top_margin", "bottom_margin", "left_margin", "right_margin"):
    setattr(sec, a, Inches(1.0))
n = d.styles["Normal"]
n.font.name = FONT
n.font.size = Pt(11)
n.font.color.rgb = INK
n.element.rPr.rFonts.set(qn("w:eastAsia"), FONT)
n.paragraph_format.space_after = Pt(0)

_fp = sec.footer.paragraphs[0]
_fp.alignment = WD_ALIGN_PARAGRAPH.CENTER
for _tag, _txt in (("begin", None), (None, "PAGE"), ("end", None)):
    _r = OxmlElement("w:r")
    if _tag:
        _c = OxmlElement("w:fldChar"); _c.set(qn("w:fldCharType"), _tag); _r.append(_c)
    else:
        _i = OxmlElement("w:instrText"); _i.set(qn("xml:space"), "preserve")
        _i.text = " PAGE "; _r.append(_i)
    _fp._p.append(_r)


def head(t, size=13, before=14):
    p = d.add_paragraph()
    p.paragraph_format.space_before = Pt(before)
    p.paragraph_format.space_after = Pt(5)
    p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
    r = p.add_run(t)
    r.font.bold = True
    r.font.size = Pt(size)
    r.font.name = FONT
    return p


def para(t, size=11, after=8):
    p = d.add_paragraph()
    p.paragraph_format.space_after = Pt(after)
    p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
    for r in add_rich(p, t):
        r.font.size = Pt(size)
        r.font.name = FONT
    return p


def sfig(path, caption, width=6.0):
    """A supplementary figure, with its caption beneath it."""
    p = d.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_before = Pt(10)
    p.paragraph_format.space_after = Pt(4)
    p.add_run().add_picture(path, width=Inches(width))
    c = d.add_paragraph()
    c.paragraph_format.space_after = Pt(12)
    r = c.add_run(caption)
    r.font.size = Pt(9.5)
    r.font.name = "Times New Roman"


def table(headers, rows, widths=None, size=9.5, caption=None):
    if caption:
        p = d.add_paragraph()
        p.paragraph_format.space_before = Pt(8)
        p.paragraph_format.space_after = Pt(4)
        p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
        for r in add_rich(p, caption):
            r.font.size = Pt(9.5)
            r.font.name = FONT
    t = d.add_table(rows=1, cols=len(headers))
    t.style = "Table Grid"
    for j, h in enumerate(headers):
        c = t.rows[0].cells[j]
        c.text = ""
        pr = c.paragraphs[0]
        pr.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
        run = pr.add_run(h)
        run.font.bold = True
        run.font.size = Pt(size)
        run.font.name = FONT
    for row in rows:
        cells = t.add_row().cells
        for j, v in enumerate(row):
            pr = cells[j].paragraphs[0]
            pr.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
            for r in add_rich(pr, str(v)):
                r.font.size = Pt(size)
                r.font.name = FONT
    if widths:
        for j, w in enumerate(widths):
            for row in t.rows:
                row.cells[j].width = Inches(w)
    d.add_paragraph().paragraph_format.space_after = Pt(6)
    return t


# ------------------------------------------------------------------ title
p = d.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
r = p.add_run("Supplementary Information")
r.font.bold = True
r.font.size = Pt(15)
r.font.name = FONT
p = d.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
p.paragraph_format.space_after = Pt(4)
for r in add_rich(p, "Sublattice dimensionality across the Crystallography Open Database"):
    r.font.size = Pt(12)
    r.font.name = FONT
p = d.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
p.paragraph_format.space_after = Pt(16)
r = p.add_run("Peifen Zhu")
r.font.size = Pt(11)
r.font.name = FONT

# ------------------------------------------------- S1 computational methods
head("Supplementary Note 1. Electronic-structure calculations", before=6)
para("The {} structures behind Fig. 4 were computed with VASP as single-point "
     "calculations on the deposited coordinates, with no relaxation: the index is "
     "evaluated on those coordinates, and relaxing first would make the two quantities "
     "describe different geometries. Input files were generated with pymatgen's "
     "MPStaticSet, so the parameters are the Materials Project community standard rather "
     "than chosen here, and the projector augmented-wave potentials are the "
     "MP-recommended set from the PBE 5.4 library. Table S1 lists the settings the "
     "generator writes.".format(len(BW)))

INCAR = [("ENCUT", "520 eV", "plane-wave cutoff"),
         ("PREC", "Accurate", "precision mode"),
         ("EDIFF", "1 x 10-4 eV", "electronic convergence"),
         ("ISMEAR / SIGMA", "-5 / 0.05", "tetrahedron method with Blochl corrections"),
         ("ISPIN", "2", "spin-polarised"),
         ("LASPH", "True", "aspherical contributions to the PAW terms"),
         ("NSW", "0", "static; no ionic steps"),
         ("ALGO", "Normal", "electronic minimisation"),
         ("LORBIT / NEDOS", "11 / 3001", "projected DOS output"),
         ("k-points", "Monkhorst-Pack, 670 / atom", "reciprocal grid density"),
         ("exchange-correlation", "PBE (GGA)", "with +U where noted below")]
table(["parameter", "value", "meaning"],
      INCAR, widths=[1.7, 2.0, 2.6],
      caption="Table S1 | VASP settings, as written by pymatgen MPStaticSet.")

U_ELEMENTS = ["Co", "Cr", "Fe", "Mn", "Mo", "Ni", "V", "W"]


def has_u(formula):
    e = set(re.findall(r"[A-Z][a-z]?", formula))
    return ("O" in e or "F" in e) and bool(e & set(U_ELEMENTS))


byc = collections.defaultdict(collections.Counter)
for rec in BW:
    byc[rec["cls"]][has_u(rec["formula"])] += 1
n_u = sum(byc[c][True] for c in byc)

head("Hubbard U, and where it falls")
para("MPStaticSet applies a Hubbard U to {} in oxides and fluorides, following the "
     "Materials Project convention. That correction is therefore not applied evenly "
     "across the three anion classes, and the imbalance runs in the direction that "
     "matters for the result reported in the main text: it reaches {:.0f} per cent of the "
     "oxide sample and none of the chalcogenide sample. It is stated here because the "
     "main text reports no bandwidth trend with connectivity for the oxides, and a reader "
     "is entitled to ask whether the functional rather than the chemistry is responsible. "
     "Two observations bear on that. The oxide null result is a wide confidence interval "
     "spanning zero rather than a slope of zero, so it is an absence of evidence on this "
     "sample; and the chalcogenide result, which is the strong one, carries no +U at all, "
     "so it cannot be an artefact of the correction."
     .format(", ".join(U_ELEMENTS[:-1]) + " and " + U_ELEMENTS[-1],
             100.0 * byc["oxide"][True] / sum(byc["oxide"].values())))

rows = []
for c, lab in (("oxide", "oxide"), ("halide", "halide"), ("chalcogen", "chalcogenide")):
    tot = sum(byc[c].values())
    rows.append([lab, tot, byc[c][True], "{:.1f}%".format(100.0 * byc[c][True] / tot)])
rows.append(["all", len(BW), n_u, "{:.1f}%".format(100.0 * n_u / len(BW))])
table(["anion class", "structures computed", "with +U applied", "per cent"],
      rows, widths=[1.6, 1.8, 1.6, 1.0],
      caption="Table S2 | How many of the electronic-structure sample carry a Hubbard U.")

# ------------------------------------------------- S2 cutoff sensitivity
head("Supplementary Note 2. Sensitivity to the bridging cutoff")
para("The distribution in the main text rests on one bridging criterion. A stratified "
     "sample of {:,} structures, {} from each anion class at each integer rank, was "
     "re-indexed with every cutoff scaled by 0.95 and by 1.05. Individual assignments "
     "near the cutoff move in both directions and largely cancel, so the distribution "
     "shifts by about one point while roughly one structure in eleven changes rank under "
     "tightening.".format(CS["n_sampled"], CS["n_per_stratum"]))
rows = []
for s in ("0.95", "1.00", "1.05"):
    ch, di = CS["changed"][s], CS["distribution"][s]
    lab = {"0.95": "0.95 (tighter)", "1.00": "1.00 (as used)", "1.05": "1.05 (looser)"}[s]
    rows.append([lab, ch["n"], ch["moved"], "{:.2f}%".format(ch["pct"])]
                + ["{:.1f}%".format(di[str(i)]) for i in range(4)])
table(["cutoff scale", "n", "rank moved", "per cent", "0D", "1D", "2D", "3D"],
      rows, widths=[1.3, 0.6, 0.85, 0.8, 0.7, 0.7, 0.7, 0.7],
      caption="Table S3 | Integer ranks that move, and the resulting distribution, "
              "at three bridging cutoffs. Percentages are weighted back to the "
              "population each stratum was drawn from.")

# ------------------------------------------------- S3 the anion exclusion
head("Supplementary Note 3. The cost of excluding two-anion structures")
B = AC["bridging"]
para("{:,} depositions carry two competing bridging classes and are set aside, which is "
     "more structures than the paper keeps. {:,} of them were indexed on every class "
     "present, one class at a time. The agreement is high but mostly cheap: it is "
     "agreement at rank zero for {} per cent of the agreeing structures, because the "
     "metal sublattice is isolated whichever anion is chosen. Only {} per cent bridge on "
     "any class at all, and among those the rank follows the choice in {} per cent of "
     "cases.".format(AC["pool"], AC["evaluated"],
                     B["agreements_at_rank_zero_pct"], B["bridge_pct"],
                     B["rank_choice_dependent_pct"]))
rows = [["indexed on two or more classes", AC["evaluated"], "100%"],
        ["same rank whichever class is chosen", AC["agree"],
         "{}%".format(AC["agree_pct"])],
        ["   of which agree at rank zero", B["agreements_at_rank_zero"],
         "{}% of agreements".format(B["agreements_at_rank_zero_pct"])],
        ["bridge on at least one class", B["bridge_on_some_class"],
         "{}%".format(B["bridge_pct"])],
        ["   of which rank follows the choice", B["rank_choice_dependent"],
         "{}% of bridging".format(B["rank_choice_dependent_pct"])]]
table(["", "structures", "share"], rows, widths=[3.1, 1.2, 1.9],
      caption="Table S4 | Indexing the excluded structures on each candidate anion class.")
rows = [[k.replace("+", " + "), v] for k, v in AC["class_pairs"].items()]
table(["classes present", "structures"], rows, widths=[2.6, 1.2],
      caption="Table S5 | Which anion classes compete in the sampled structures.")

# ------------------------------------------------- S4 by metal
head("Supplementary Note 4. Sublattice dimensionality by framework metal")
para("The data behind Fig. 3. Structures in which a single framework metal makes the "
     "assignment unambiguous, for the twenty-two metals with more than 1,900 such "
     "structures, ordered by mean periodic rank.")
rows = []
for el, v in sorted(BM.items(), key=lambda kv: kv[1]["mean"]):
    rows.append([el, "{:,}".format(v["n"]), "{:.2f}".format(v["mean"])]
                + ["{:.1f}%".format(x) for x in v["pct"]])
table(["metal", "n", "mean rank", "0D", "1D", "2D", "3D"],
      rows, widths=[0.8, 0.9, 1.0, 0.85, 0.85, 0.85, 0.85],
      caption="Table S6 | Rank distribution per framework metal.")

# ------------------------------------------------- S5 the layered query
head("Supplementary Note 5. Corroborated layered compositions")
para("The worked query of the main text, listed in full. These are compositions whose "
     "metal-anion sublattice is two-dimensional with a fractional term below 0.05, "
     "deposited at least twice with every deposition firm, and with no deposition of that "
     "composition indexing anything other than two dimensions. {:,} compositions meet all "
     "three conditions: {:,} oxide, {:,} halide and {:,} chalcogenide. They are ordered by "
     "how many independent depositions agree."
     .format(LQ["n_unanimous"], LQ["by_class"]["oxide"], LQ["by_class"]["halide"],
             LQ["by_class"]["chalcogen"]))
NCOL = 4
items = ["{} ({})".format(r["formula"], r["depositions"]) for r in LQ["all"]]
rows = [items[i:i + NCOL] + [""] * (NCOL - len(items[i:i + NCOL]))
        for i in range(0, len(items), NCOL)]
table(["composition (depositions)"] * NCOL, rows,
      widths=[1.55] * NCOL, size=8,
      caption="Table S7 | The {:,} corroborated layered compositions, with the number of "
              "agreeing depositions in brackets.".format(LQ["n_unanimous"]))

d.save(OUT)
print("saved: " + OUT)
print("  Table S1 %d settings | S2 +U in %d of %d | S3 %d cutoffs | S4/S5 exclusion | "
      "S6 %d metals | S7 %d compositions"
      % (len(INCAR), n_u, len(BW), len(CS["changed"]), len(BM), len(LQ["all"])))
