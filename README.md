# Sublattice dimensionality across the Crystallography Open Database

Indexed database, code and computed results behind the database paper on sublattice
dimensionality. Structures are identified throughout by Crystallography Open Database
entry number rather than redistributed, so everything here can be rejoined to the source
archive without carrying a copy of it.

The companion methods paper, *A sublattice dimensionality index validated against
deposited metal halides*, defines and validates the index and has its own deposit. This
one is the database-scale application.

## What the index is

For a chosen metal–ligand sublattice, two metal centres are bridged when they share a
ligand within the first coordination shell. Each bridge carries a weight

    w = m · exp[−(d − d0) / λ]

with *m* the number of shared bridging ligands, *d* the metal–metal separation, *d0* the
shortest bridged separation in the structure and λ = 1.0 Å a coupling decay length. The
index is

    D = d_top + w_next / (w_next + w̃)

where `d_top` is the periodic rank of the bridging network — the rank of the lattice
translation subgroup under which the connected component maps onto itself — and the
fractional term measures how close the structure sits to gaining the next dimension.

Two conventions are applied throughout: the halides are taken as a single bridging class,
and the framework metals are taken as a single sublattice.

## The indexed database

`data/cod_dimensionality.jsonl.gz` holds one JSON record per line, 521,901 of them, one
for every entry in the bulk archive. 162,587 carry a computed index; the rest are
recorded with `D` null and are excluded for structural reasons rather than by failure:

| reason | count |
|---|---|
| no framework metal | 205,706 |
| bridging anion class not defined | 148,440 |
| more than 90 metal centres | 2,069 |
| unreadable | 3,099 |

Fields per record:

| field | meaning |
|---|---|
| `cod` | Crystallography Open Database entry number |
| `formula` | formula as deposited |
| `D` | the index; null where not computed |
| `d_top` | integer periodic rank of the bridging network |
| `frac` | the fractional term, in [0,1) |
| `anion_class` | oxide, halide or chalcogen |
| `nsite` | sites in the cell |
| `nmetal` | framework metal centres |
| `nbridge` | bridges in the network |

## Layout

    code/     the indexer and every script that produced a number in the paper
    data/     the indexed database and the computed results, as JSON
    figures/  the rendered figures and the json they are drawn from

### code

| file | what it does |
|---|---|
| `dimindex.py` | reference implementation of the index on a named sublattice |
| `cod_wide_index.py` | database-scale indexer; chooses the anion class per structure |
| `wholenet_scale.py` | whole-network comparison on a stratified random sample |
| `baseline_and_scale.py` | stoichiometric baseline and whole-network ranks |
| `setup_bandwidth_scale.py` | writes the single-point inputs for the stratified sample |
| `analyze_bw_scale.py` | band edges and widths from EIGENVAL, per spin channel |
| `export_polyhedra.py` | coordination polyhedra and repeat directions for the Fig. 1 panels |

### data

| file | contents |
|---|---|
| `cod_dimensionality.jsonl.gz` | the indexed database, one record per structure |
| `npj_stats.json` | totals, skip reasons and the rank distribution by anion class |
| `npj_figdata.json` | panel data for Fig. 2 |
| `npj_bymetal.json` | rank distribution per framework metal, Fig. 3 |
| `npj_bw.json` | bandwidth against the index, Fig. 4 |
| `bw_scale_results.json` | per-structure band widths, gaps and k-point counts |
| `npj_disagreements.json` | compositions whose depositions disagree, Fig. 5 |
| `npj_wholenet.json` | whole-network agreement, summarised |
| `wholenet_scale.json` | the per-structure whole-network comparison |
| `npj_panels.json` | polyhedra, cell edges and lattice for the Fig. 1 structures |
| `baseline_scale.json` | the density and stoichiometric baselines |

### figures

The five rendered figures at 600 d.p.i., together with the json they are drawn from. The
plotting scripts themselves are MATLAB and are not distributed; every number they draw is
in data/ and code/.

## Reproducing the numbers

Requires Python 3.11 with `pymatgen` (2026.5.18) and `numpy`; `ase` (3.26.0) for the
whole-network comparison, `scipy` for the polyhedra. Structures come from a local copy of
the Crystallography Open Database bulk archive.

    python code/cod_wide_index.py <cod cif root> cod_dimensionality.jsonl
    python code/wholenet_scale.py  cod_dimensionality.jsonl wholenet_scale.json
    python code/baseline_and_scale.py cod_dimensionality.jsonl <cache> baseline_scale.json
    python code/setup_bandwidth_scale.py cod_dimensionality.jsonl <rundir>
    python code/analyze_bw_scale.py <rundir> -o bw_scale_results.json
    python code/export_polyhedra.py npj -o npj_panels.json

The single-point electronic-structure step needs VASP and is the only part that does not
run from this archive alone; `setup_bandwidth_scale.py` writes the inputs and
`analyze_bw_scale.py` reads the EIGENVAL files it produces.

## Limits

Cells containing more than 90 metal centres are declined rather than evaluated. About a
third of the archive carries partial occupancies and is indexed on the deposited
coordinates without an ordering approximation; the electronic-structure sample is
restricted to ordered structures for that reason. The index describes the sublattice it
is given and makes no claim about the rest of the structure.

## Licence

Code under the MIT licence; data under CC BY 4.0. See LICENSE.txt.
