# Cross-check of TopoQSPR indices against ChemGraphX

ChemGraphX (Jacob, Clement, Arockiaraj, Jindal, *J. Cheminform.* 18:39, 2026;
MIT licence; <https://github.com/aditya-jindal/chemgraphx>) is an independent
implementation of degree, neighbourhood-degree-sum and distance indices from
the senior co-author's group. It is a consistency check, not ground truth
(`docs/recommendations-review.md`). Ground truth is still the closed forms,
brute force, networkx and Mordred checks in `docs/INDEX_VERIFICATION.md`.

## Result

* **41 TopoQSPR indices have a ChemGraphX counterpart.** On a fixed panel of
  43 graphs (9 closed-form graphs, 10 benzenoid PAHs and the 24 approved
  drugs of `tests/fixtures/indices/drugs.json`) all 41 agree on all 43
  graphs (1,763 comparisons). The largest relative difference is
  1.3e-15 (Randić index), which is floating-point summation order.
  40 indices agree directly. One, the edge Wiener index, agrees after a
  documented convention offset of m(m-1)/2.
* **No TopoQSPR bug was found and no TopoQSPR formula was changed.**
* Every difference found is a naming or convention difference, or a matter
  of how ChemGraphX handles its input. Each one has a minimal example below
  and a test in `tests/test_chemgraphx_crosscheck.py`. Inside the domain
  TopoQSPR accepts (connected, hydrogen-suppressed, at least 3 heavy atoms)
  none of them changes a value, except for two points. ChemGraphX's
  unqualified "PI" is TopoQSPR's `PI_v`. ChemGraphX's edge Wiener index is
  TopoQSPR's `W_e` minus m(m-1)/2.

## Version and reproduction

| Item | Value |
|---|---|
| ChemGraphX commit | `d8fea04ffc90ec4acb93748a79392c45c93adc1d` (branch `main`, cloned 2026-10-02) |
| Working tree | identical to that commit: the blob hashes of all 52 tracked files match the git index |
| Files exercised | `server/graph.py`, `basic_values.py`, `four_columns.py`, `column_ab.py`, `column_cd.py`, `column_e.py`, `column_f.py` (SHA-256 recorded in the fixture) |
| Runtime used for ChemGraphX | python-igraph 0.11.8, numpy 2.2.6, networkx 3.4.2, matplotlib 3.9.2, in a separate virtualenv layered on the TopoQSPR env |
| TopoQSPR side | rdkit 2026.3.6, chembl_structure_pipeline 1.2.4 |

```
# needs a ChemGraphX clone and python-igraph
PYTHONPATH=src python tools/crosscheck_chemgraphx.py --chemgraphx /path/to/chemgraphx --write-fixture
# offline tests (use the recorded fixture)
PYTHONPATH=src python -m pytest -q -m "not network" tests/test_chemgraphx_crosscheck.py
# re-run ChemGraphX and compare with the fixture (skipped without a clone or igraph)
CHEMGRAPHX_DIR=/path/to/chemgraphx PYTHONPATH=src python -m pytest -q tests/test_chemgraphx_crosscheck.py
```

The tool writes `tests/fixtures/chemgraphx/values.json`. That file holds the
graph edge lists, all 76 ChemGraphX index values per graph, and the ChemGraphX
output on the minimal examples. Two runs produce byte-identical files. The
fixture contains only numbers computed by ChemGraphX's MIT-licensed code and
no ChemGraphX source.

### How the identical graph is fed to ChemGraphX

For every panel entry TopoQSPR builds its graph: `standardize(smiles).graph`
for molecules, `MolGraph.from_edges` for closed-form graphs. ChemGraphX then
receives the same graph through three routes:

1. a `.mol` V2000 file (`read_mol_v2`),
2. a `.mol` V3000 file (`read_mol_v3`),
3. the adjacency-list JSON of the text-area route (`create_graph`). The web
   client builds this from a pasted adjacency matrix (`client/src/tabs/HomeTab.jsx`,
   `convertMatrixToList`).

For molecules the molblock is RDKit's for the hydrogen-free standardised
parent, with heteroatoms and kekulised bond orders kept. ChemGraphX ignores
both. For closed-form graphs the molblock is all-carbon with single bonds. The
classes are called as `server/main.py` calls them (`ColumnF(file, textarea)`
and `ColumnCD(file, textarea)`). Values are read from the instance dicts
`degree_indices` / `distance_indices`, which the endpoints only relabel
(`capitalize_dictionary`). On every route ChemGraphX's graph is checked to be
isomorphic to TopoQSPR's (networkx `vf2pp_is_isomorphic`). The three routes
agree to a relative 2.0e-16. The only non-zero gap, in paclitaxel, is one ulp
caused by a different edge order.

## 1. What ChemGraphX computes (as implemented at d8fea04)

### Degree-based block (`four_columns.py:compute_matrices`)

For every edge uv, with vertex weights x_u, x_v, ChemGraphX evaluates 26
edge functions and sums them over E. These are the "Degree Based Indices"
(`ColumnAB`, x = d, the degree) and the "Degree Sum Based Indices"
(`ColumnCD`, x = S, where `S_u = sum_{w ~ u} d_w`, from
`get_weight_edge_list`). Writing a = x_u, b = x_v:

| ChemGraphX key | edge function as implemented | TopoQSPR (x = d) | TopoQSPR (x = S) |
|---|---|---|---|
| `first_zagreb` | a + b | `M1` | no counterpart |
| `second_zagreb` | a b | `M2` | no counterpart |
| `randić` | 1/sqrt(a b) | `R` | `S_R` |
| `atom_bond_connectivity` | sqrt((a + b - 2)/(a b)) | `ABC` | `S_ABC` |
| `harmonic` | 2/(a + b) | `H` | `S_H` |
| `sum_connectivity` | 1/sqrt(a + b) | `SC` | `S_SC` |
| `hyper_zagreb` | (a + b)^2 | `HM` | no counterpart |
| `geometric_arithmetic` | 2 sqrt(a b)/(a + b) | `GA` | `S_GA` |
| `irregularity_measure` | \|a - b\| (Albertson) | no counterpart | no counterpart |
| `sigma` | (a - b)^2 | `sigma` | no counterpart |
| `somber` | sqrt(a^2 + b^2) (Sombor) | `SO` | no counterpart |
| `forgotten` | (a + b)^2 - 2ab = a^2 + b^2 | `F` | no counterpart |
| `symmetric_division_degree` | ((a + b)^2 - 2ab)/(a b) | `SDD` | `S_SDD` |
| `augmented_zagreb` | (a b/(a + b - 2))^3 | `AZI` | no counterpart |
| `bi_zagreb` | a + b + a b | `BM` | no counterpart |
| `tri_zagreb` | a^2 + b^2 + a b | `TM` | no counterpart |
| `geometric_harmonic` | sqrt(a b)(a + b)/2 | `GH` | no counterpart |
| `geometric_bi_zagreb` | sqrt(a b)/(a + b + a b) | `GBM` | `S_GBM` |
| `geometric_tri_zagreb` | sqrt(a b)/(a^2 + b^2 + a b) | `GTM` | no counterpart |
| `harmonic_geometric` | 2/(sqrt(a b)(a + b)) | `HG` | no counterpart |
| `harmonic_bi_zagreb` | 2/((a + b)(a + b + a b)) | `HBM` | no counterpart |
| `harmonic_tri_zagreb` | 2/((a + b)(a^2 + b^2 + a b)) | `HTM` | no counterpart |
| `bi_zagreb_geometric` | (a + b + a b)/sqrt(a b) | `BMG` | no counterpart |
| `bi_zagreb_harmonic` | (a + b + a b)(a + b)/2 | `BMH` | no counterpart |
| `tri_zagreb_geometric` | (a^2 + b^2 + a b)/sqrt(a b) | `TMG` | no counterpart |
| `tri_zagreb_harmonic` | (a + b)(a^2 + b^2 + a b)/2 | `TMH` | `S_TMH` |

Of these, 25 + 8 = 33 have a TopoQSPR counterpart. TopoQSPR has two
neighbourhood indices that ChemGraphX does not compute: `S_RR` (sqrt(S_u S_v))
and `S_ReZ3` (S_u S_v (S_u + S_v)).

### Distance block (`column_e.py`)

The distance matrix D comes from igraph `shortest_paths()` and is stored as
`float32` (`SAVE_MEM`). For an edge e = uv: n_u = #{x : D(x,u) < D(x,v)}.
The vertex-edge distance is D(x, f) = min(D(x,a), D(x,b)) for f = ab, and
m_u = #{f : D(u,f) < D(v,f)}. These are the same definitions TopoQSPR uses
(`GraphContext.edge_split`).

| ChemGraphX key | as implemented | TopoQSPR |
|---|---|---|
| `wiener_index` | sum(D)/2 | `W` |
| `szeged` | sum n_u n_v | `Sz` |
| `edge_szeged` | sum m_u m_v | `Sz_e` |
| `PI` | sum (n_u + n_v) | **`PI_v`** (vertex PI) |
| `edge_PI` | sum (m_u + m_v) | **`PI`** (Khadikar PI) |
| `mostar` | sum \|n_u - n_v\| | `Mo` |
| `edge_mostar` | sum \|m_u - m_v\| | `Mo_e` |
| `edge_wiener` | (1/2) sum over ordered edge pairs (e, f), e = f included, of min(D(a,x), D(a,y), D(b,x), D(b,y)) | `W_e` - m(m-1)/2 (see D2) |
| `vertex_edge_wiener` | (1/2) sum_x sum_f D(x, f) | no counterpart (not registered, see INDEX_VERIFICATION "Not registered") |
| `edge_vertex_szeged` | sum (n_u m_v + n_v m_u)/2 | no counterpart |
| `total_szeged` | Sz + Sz_e + 2 Sz_ev | no counterpart |
| `total_mostar` | sum \|(n_u + m_u) - (n_v + m_v)\| | no counterpart |
| `<X>_plus`, X in {szeged, edge_szeged, PI, edge_PI, mostar, edge_mostar} | sum (d_u + d_v) * (edge term of X) | no counterpart (6 indices) |
| `<X>_mul`, same six X | sum d_u d_v * (edge term of X) | no counterpart (6 indices) |

That makes 24 distance indices, 8 of them with a TopoQSPR counterpart.
TopoQSPR's `J`, `DD` and `Gut` have no ChemGraphX counterpart.

### Entropies (not compared)

For every degree, degree-sum and distance edge function f with total
I = sum_e f(e), ChemGraphX reports -sum_e (f(e)/I) ln(f(e)/I) over edges with
f(e) != 0 (`entropy_sum`, natural log). It returns the string
`"indeterminate"` when I = 0, for example `sigma` on a regular graph. That
gives 26 + 26 + 21 entropies (the 21 distance entropies exclude `wiener_index`,
`vertex_edge_wiener` and `edge_wiener`), so the code totals are 52 degree
indices + 52 degree entropies + 24 distance indices + 21 distance entropies
= 149. This matches the paper's table total of 149, not the "150" in its text.
TopoQSPR's only entropy, `I_deg`, is a vertex-degree entropy in log2 and has
no ChemGraphX counterpart. The `/graph_information` endpoint reports n, m and
the degree and degree-sum edge partitions, which are not indices.

### Input handling

* **.mol V2000** (`read_mol_v2`): atom and bond counts come from fixed
  columns 0-2 and 3-5 of the counts line. Every atom of the atom block is a
  vertex. Every bond line adds an edge, whatever its type: single, double,
  triple and aromatic (type 4) are all treated the same, and duplicates
  collapse in `nx.Graph`. Element, charge, isotope, stereo and coordinates
  are not read. **Hydrogens are not removed**: an explicit H atom in the file
  becomes a vertex (D6).
* **.mol V3000** (`read_mol_v3`): counts come from `M  V30 COUNTS`, bonds from
  the `BEGIN BOND` block (fields 5-6). The handling is otherwise the same as
  V2000. The format is chosen by the regex `V(\d+)` on line 4; any other
  version returns a string, not a graph.
* **.pdb**: the `CONECT` records. **.txt**: one row of neighbours per vertex.
  Neither route was exercised here.
* **Adjacency JSON** (text area): edges only. A vertex with no neighbours is
  not added, so an all-zero row of a pasted matrix disappears.
* Vertices are relabelled 1..n in insertion order. The distance code then
  builds an igraph graph from the **edge list only**, with vertex ids 0..max
  label, and drops row/column 0 (D5).

## 2. Agreement table

Relative difference |ChemGraphX - TopoQSPR| / |TopoQSPR| on the 43 panel
graphs. The pass threshold is 1e-10, the same as the Mordred check.

| ChemGraphX (block.key) | TopoQSPR | relation | graphs | max rel. diff | status |
|---|---|---|---|---|---|
| `degree.first_zagreb` | `M1` | equal | 43 | 0.0e+00 | agree |
| `degree.second_zagreb` | `M2` | equal | 43 | 0.0e+00 | agree |
| `degree.randić` | `R` | equal | 43 | 1.3e-15 | agree |
| `degree.atom_bond_connectivity` | `ABC` | equal | 43 | 1.1e-15 | agree |
| `degree.harmonic` | `H` | equal | 43 | 7.9e-16 | agree |
| `degree.sum_connectivity` | `SC` | equal | 43 | 3.9e-16 | agree |
| `degree.hyper_zagreb` | `HM` | equal | 43 | 0.0e+00 | agree |
| `degree.geometric_arithmetic` | `GA` | equal | 43 | 4.1e-16 | agree |
| `degree.sigma` | `sigma` | equal | 43 | 0.0e+00 | agree |
| `degree.somber` | `SO` | equal | 43 | 9.9e-16 | agree |
| `degree.forgotten` | `F` | equal | 43 | 0.0e+00 | agree |
| `degree.symmetric_division_degree` | `SDD` | equal | 43 | 2.8e-16 | agree |
| `degree.augmented_zagreb` | `AZI` | equal | 43 | 2.0e-16 | agree |
| `degree.bi_zagreb` | `BM` | equal | 43 | 0.0e+00 | agree |
| `degree.tri_zagreb` | `TM` | equal | 43 | 0.0e+00 | agree |
| `degree.geometric_harmonic` | `GH` | equal | 43 | 6.7e-16 | agree |
| `degree.geometric_bi_zagreb` | `GBM` | equal | 43 | 6.3e-16 | agree |
| `degree.geometric_tri_zagreb` | `GTM` | equal | 43 | 8.2e-16 | agree |
| `degree.harmonic_geometric` | `HG` | equal | 43 | 5.6e-16 | agree |
| `degree.harmonic_bi_zagreb` | `HBM` | equal | 43 | 9.5e-16 | agree |
| `degree.harmonic_tri_zagreb` | `HTM` | equal | 43 | 6.8e-16 | agree |
| `degree.bi_zagreb_geometric` | `BMG` | equal | 43 | 6.9e-16 | agree |
| `degree.bi_zagreb_harmonic` | `BMH` | equal | 43 | 0.0e+00 | agree |
| `degree.tri_zagreb_geometric` | `TMG` | equal | 43 | 7.7e-16 | agree |
| `degree.tri_zagreb_harmonic` | `TMH` | equal | 43 | 0.0e+00 | agree |
| `degree_sum.randić` | `S_R` | equal | 43 | 9.1e-16 | agree |
| `degree_sum.atom_bond_connectivity` | `S_ABC` | equal | 43 | 4.4e-16 | agree |
| `degree_sum.harmonic` | `S_H` | equal | 43 | 1.1e-15 | agree |
| `degree_sum.sum_connectivity` | `S_SC` | equal | 43 | 6.6e-16 | agree |
| `degree_sum.geometric_arithmetic` | `S_GA` | equal | 43 | 7.9e-16 | agree |
| `degree_sum.symmetric_division_degree` | `S_SDD` | equal | 43 | 4.1e-16 | agree |
| `degree_sum.geometric_bi_zagreb` | `S_GBM` | equal | 43 | 4.9e-16 | agree |
| `degree_sum.tri_zagreb_harmonic` | `S_TMH` | equal | 43 | 0.0e+00 | agree |
| `distance.wiener_index` | `W` | equal | 43 | 0.0e+00 | agree |
| `distance.szeged` | `Sz` | equal | 43 | 0.0e+00 | agree |
| `distance.edge_szeged` | `Sz_e` | equal | 43 | 0.0e+00 | agree |
| `distance.PI` | `PI_v` | equal | 43 | 0.0e+00 | agree |
| `distance.edge_PI` | `PI` | equal | 43 | 0.0e+00 | agree |
| `distance.mostar` | `Mo` | equal | 43 | 0.0e+00 | agree |
| `distance.edge_mostar` | `Mo_e` | equal | 43 | 0.0e+00 | agree |
| `distance.edge_wiener` | `W_e` | W_e - m(m-1)/2 | 43 | 0.0e+00 | agree |

Panel: P3, P5, P8, C5, C6, C9, K1,4, K4, K2,3; benzene, naphthalene,
anthracene, phenanthrene, tetracene, pyrene, chrysene, triphenylene,
benzo[a]pyrene, coronene; and aspirin through paclitaxel (24 drugs, 9 to 101
heavy atoms; the largest Wiener index is vancomycin's, 52,810). The degree-sum
`S_ABC` agreement also confirms the v0.x numerator defect fixed in v1 (see
INDEX_VERIFICATION, "v0.x defects"): ChemGraphX uses S_u + S_v - 2, like v1.

## 3. Differences, root causes and minimal examples

Each item names the side it concerns, the root cause and the smallest input
that shows it. ChemGraphX values are from commit d8fea04, recorded in the
fixture's `probes`.

### D1. PI naming (convention, values agree)

ChemGraphX's unqualified `PI` is sum (n_u + n_v), the vertex PI index
(TopoQSPR `PI_v`, khalifeh2008pi). Its `edge_PI` is sum (m_u + m_v), Khadikar's
original PI index (TopoQSPR `PI`, khadikar2000). The values agree once mapped.
Joining the two tools' outputs by column name would swap them.
Minimal example: **P3**. TopoQSPR gives `PI` = 2 and `PI_v` = 6; ChemGraphX
gives `PI` = 6 and `Edge PI` = 2. Test: `test_pi_naming_is_swapped_on_p3`.

### D2. Edge Wiener index (convention, offset C(m,2))

TopoQSPR registers W_e = W(L(G)), the Wiener index of the line graph
(iranmanesh2009). For e != f, D_L(e,f) = min endpoint distance + 1. ChemGraphX
uses the minimum endpoint distance itself: adjacent edges are at distance 0,
and the e = f diagonal adds 0. On a connected graph the two differ by exactly
C(m,2) = m(m-1)/2. The minimum-distance convention also occurs in the
literature (often written W_e0). That attribution was not re-checked against a
primary source here, so neither convention is called wrong. Minimal examples:
**P3** gives TopoQSPR 1 (L(P3) = P2) and ChemGraphX 0. **P4** gives 4 and 1.
On all 43 panel graphs, ChemGraphX = `W_e` - C(m,2) exactly. A networkx
computation of the minimum-endpoint sum confirms the identity on the panel.
Tests: `test_edge_wiener_convention_on_paths`,
`test_edge_wiener_offset_identity_on_the_panel`.

### D3. Undefined values (convention)

On an edge with a zero denominator, TopoQSPR returns NaN plus a reason.
ChemGraphX returns numpy `inf` or `nan` and emits a RuntimeWarning, which is
not shown to the user. Minimal example: **K2** (one edge, d = (1, 1)). AZI
divides by d_u + d_v - 2 = 0: TopoQSPR gives NaN ("division by zero ..."),
ChemGraphX gives `inf`. In TopoQSPR's domain (n >= 3, connected) no degree or
degree-sum denominator can be 0, so this never touches the panel. For entropies
ChemGraphX returns the string `"indeterminate"` when the index total is 0.
Test: `test_k2_azi_is_nan_here_and_infinite_in_chemgraphx`.

### D4. Disconnected graphs (ChemGraphX input handling)

TopoQSPR makes every distance index NaN on a disconnected graph, with the
reason "graph has k components". ChemGraphX returns W = `inf` and
vertex-edge / edge Wiener = `inf`. Its Szeged, PI and Mostar values stay finite,
because a vertex in another component is at distance inf from both ends of
an edge and so counts as equidistant. Minimal examples: **2K2** gives
W = inf, Sz = 2, PI(vertex) = 4. **P3 + K2** gives W = inf, Sz = 5,
PI(vertex) = 8, Mo = 2. TopoQSPR rejects such molecules earlier
(`mixture` / `disconnected_graph`). Test:
`test_disconnected_graph_distance_indices`.

### D5. An isolated atom's position changes ChemGraphX's Wiener index (ChemGraphX graph construction)

Root cause: `column_e.py:nx_to_igraph` builds the igraph graph from the edge
list only (`igraph(edges=...)`). igraph then has vertices 0..max(label in an
edge). Row and column 0 are dropped as the unused label 0. If the isolated
vertex has the largest label, it is simply missing from the distance matrix,
and W is finite and computed without it. If it has any other label, it is
present at distance inf, and W = inf. The vertex count (5) and the degree
indices are the same either way. Minimal example: sodium acetate written as a
`.mol` file. With **`CC(=O)[O-].[Na+]`** (Na is atom 5) ChemGraphX gives
W = 9 and vertex-edge Wiener 3. With **`[Na+].CC(=O)[O-]`** (Na is atom 1) it
gives W = inf and vertex-edge Wiener inf. Acetic acid alone gives W = 9.
TopoQSPR standardises both inputs to the acetic acid parent (n = 4, W = 9).
Passing the vertex count to igraph (`igraph(n=..., edges=...)`) would make the
result independent of atom order. This is recorded as an observation and has
not been reported upstream. Test:
`test_salt_counter_ion_position_changes_chemgraphx_wiener`.

### D6. Explicit hydrogens in a .mol file become vertices (ChemGraphX input convention)

ChemGraphX takes the graph exactly as written in the file. The
hydrogen-suppressed values of the paper therefore need a hydrogen-free
`.mol`. TopoQSPR never makes H a vertex (`graph_from_mol`). Minimal examples:
**toluene** without H gives n = 7, M1 = 30 and W = 42 in both tools. With
`AddHs` (15 atoms) ChemGraphX gives n = 15, M1 = 78 and W = 318. One explicit
stereo H, as some drawing tools write at a stereocentre, is enough:
**alanine** (`C[C@H](N)C(=O)O`) gives M1 = 22 (both tools), but with one explicit H on the
stereocentre ChemGraphX gives n = 7 and M1 = 30. Test:
`test_explicit_hydrogens_become_chemgraphx_vertices`.

### D7. float32 accumulation in distance sums (ChemGraphX numerics)

ChemGraphX stores distances and the n_u / m_u counts as `float32`
(`SAVE_MEM`) and sums them in `float32`, so integer results stop being exact
once a sum passes 2^24 = 16,777,216. TopoQSPR uses exact Python integers.
Minimal examples on paths, where W(P_n) = C(n+1, 3):

| graph | exact W | ChemGraphX W | ChemGraphX Sz (exact = W) | ChemGraphX edge Wiener (exact) |
|---|---|---|---|---|
| P300 | 4,499,950 | 4,499,950 | 4,499,950 | 4,410,549 (4,410,549) |
| P600 | 35,999,900 | 35,999,908 | 35,999,900 | 35,641,100 (35,641,099) |
| P1000 | 166,666,500 | 166,666,352 | 166,666,496 | 165,668,560 (165,668,499) |

The largest drug in the panel (vancomycin, W = 52,810) is far below this
limit, so drug-sized graphs are not affected. Test:
`test_float32_accumulation_on_long_paths`.

### Other observations (no effect on the compared values)

* With Flask 3.0.3 and numpy 2.2.6, the `/distance_based` endpoint returned
  HTTP 500 (`TypeError: Object of type float32 is not JSON serializable`)
  when called through Flask's test client. The cause is that `wiener_index`
  stays a `numpy.float32` (`column_e.py`, `np.sum(self.dist_matrix) / 2`),
  while every other value is passed through `float()`. The `/degree_based` and
  `/graph_information` endpoints returned 200. The deployed Lambda service
  was not called; its library versions are unknown.
* `Graph.shortest_paths()` is deprecated in python-igraph 0.11.8 (warning
  "Graph.shortest_paths() is deprecated; use Graph.distances() instead"). Newer igraph
  releases may remove it, so the cross-check pins igraph 0.11.8.
* `server/test.py` calls `ColumnF(file, testing=True)`. `BasicValues.__init__`
  passes its second argument positionally to `Graph(file, textarea_data,
  testing=False)`, so `True` arrives as `textarea_data` and `json.loads(True)`
  raises a TypeError. The web routes are unaffected, because they pass
  `textarea_data` in that position.
* `vertex_edge_wiener` halves sum_x sum_f D(x, f), although vertex-edge pairs
  are not double-counted. Toluene gives the half-integer 28.5. TopoQSPR has no
  counterpart: W_ve is deliberately unregistered because conventions differ.
  The cited definition was not checked here, so this is recorded without a
  verdict.

## 4. Consequences for TopoQSPR

* No formula change. All 41 shared indices are independently reproduced to
  1.3e-15 on 43 graphs. This includes the 12 hybrid indices of
  arockiaraj2023sar, all checked against ChemGraphX for the first time, and 8
  neighbourhood-degree-sum indices.
* When TopoQSPR output is put next to ChemGraphX output, map `PI_v` to
  ChemGraphX "PI", map `PI` to "Edge PI", and subtract C(m,2) from `W_e`.
* To reproduce TopoQSPR values in the ChemGraphX web tool, upload a
  hydrogen-free, single-component `.mol` file (D5, D6).
