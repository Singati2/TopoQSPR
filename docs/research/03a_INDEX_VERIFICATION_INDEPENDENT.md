# 03a. Independent verification of the 104 TopoQSPR indices

Date: 2026-10-13. Branch `research-2026-10-13`. Scripts and raw outputs are in
`docs/research/index_verification/`.

## Result

- **All 104 indices agree with the independent implementation on every defined value.** That is 111,920 value comparisons on 591 graphs. The largest relative difference is 2.3e-15, which is binary64 rounding. The 67 rational-valued indices agree exactly; the 37 irrational ones differ only at that level.
- **Code bugs: 0. Unresolved: 0.** One edge-case **definition difference** (HCM2 on edgeless graphs), which no molecule can reach.
- **Undefined cases.** On 2,373 (graph, index) cells both implementations say "undefined", and the package gives a NaN with a reason every time. The package never returns NaN without a reason and never returns 0 where the value is undefined.
- **Hosoya Z** is exact as a Python `int` through `compute_indices` on every graph, including 9 cohort molecules with Z > 2^53 (largest: CHEMBL1688530, Z = 52,259,511,064,445,879,724,729,840). Through `index_table` (a float64 DataFrame) these values are rounded, with absolute error up to 1.2e9 (FLOATING POINT, see below).
- **Dependencies.** `ReZ1` equals the number of heavy atoms n on every graph without isolated vertices, so it adds nothing as a feature. For molecules with Δ = 4, the 81 degree, R_k and HCM2 indices span only 9 dimensions; for Δ = 3 they span 5. On the 3 ≤ Δ ≤ 4 atlas, all 104 indices together with 1, n and m have rank 38.

## What "independent" means here

`refimpl.py` was written from the formulas in `docs/INDEX_VERIFICATION.md`, from the LaTeX strings in `REGISTRY`, and from the standard definitions summarised below. From `topoqspr.indices` I read only the registry metadata (names, family, LaTeX, ref) and the list of function names. I did not read the bodies of the computation functions. Values come from the public API only:

- `MolGraph.from_edges` with `compute_indices` (the "direct" route, for abstract graphs);
- `prepare_structures` with `index_table` (the "smiles" route, which runs the full ChEMBL standardisation).

Arithmetic in the reference:

- Rational-valued indices use `fractions.Fraction`, so they are exact.
- Square roots and logarithms use `decimal.Decimal` at 60 significant digits.
- Graph distances come from my own BFS.
- The line graph comes from `networkx.line_graph`.
- Hosoya Z is counted two ways: brute-force include/exclude enumeration of matchings (graphs with n ≤ 10), and the vertex recursion Z(G) = Z(G−v) + Σ_{w~v} Z(G−v−w), memoised and multiplicative over components (all graphs, Python ints). The two methods agree on all 461 atlas graphs.

Where the literature has more than one convention, the reference follows the convention fixed in `docs/INDEX_VERIFICATION.md`. Agreement therefore means the code matches its documented definition. It does not re-adjudicate the decisions recorded there. Those conventions are:

- R_k: M_kR(d) = Δ − d + k if k ≤ d, otherwise (Δ − d + k) mod Δ; NaN for k > Δ.
- HCM2: c = Δ − d + 1.
- W_e = W(L(G)).
- n_u(e) counts vertices strictly closer to u. m_u(e) counts edges f = xy with min(D(u,x), D(u,y)) < min(D(v,x), D(v,y)); equidistant edges are excluded.
- J uses μ = m − n + 1.
- I_deg uses log base 2.

## Test set

| Set | Graphs | Route(s) | Notes |
|---|---|---|---|
| (a) Exhaustive: all connected graphs, n = 2..7, Δ ≤ 4 (networkx atlas) | 461 (1+2+6+21+78+353; 995 connected graphs before the Δ ≤ 4 filter) | direct + smiles | Realised as saturated hydrocarbons (RDKit RWMol, C atoms, single bonds). 0 sanitisation failures, 0 standardisation rejections, 0 non-isomorphic round trips. |
| (b) Closed-form families: P_n (n = 2..30), C_n (3..30), stars K_{1,s} (2..30) | 86 | direct; smiles for P_n, C_n and stars with s ≤ 4 | The 14 SMILES with n ≤ 7 duplicate atlas graphs (InChIKey) and are covered there. 1,615 hand-derived closed-form values (M1, M2, F, sigma, R, H, ReZ1, HCM2, W, J, Sz, PI, PI_v, Mo, Mo_e, W_e, DD, Gut, Z = F_{n+1} / L_n / s+1, I_deg) match both the reference (to 1e-40) and the package (1615/1615). |
| (c) Cohort drugs from `main_confirmatory/molecules.csv` | 36 | smiles (+ direct for Z) | Lanatoside C (CHEMBL506569), digoxin (CHEMBL1751), atracurium (CHEMBL1360); the 9 cohort molecules with package Z > 2^53; 24 molecules spread evenly over n_heavy (3..49). My own SMILES→graph is isomorphic to the package graph for all 36. |
| Undefined-value probes | 9 (8 new; `P3` duplicates the closed-form P3) | direct | K2, P3, 2K2, P3 ∪ K1, K1, empty graph, 3K1, K5 (Δ = 4), K6 (Δ = 5) |
| Cohort-wide checks on package values | 1,821 | smiles | Rings, ReZ1 = n, ranks and relations (see Dependencies) |

The atlas used here stops at n = 7 because networkx's graph atlas only goes that far. Drug-size graphs are covered by sets (b) and (c).

## Per-index results

The full table is in `index_verification/out/summary_per_index.csv`, with one row per index giving the number of graphs, comparisons, both-undefined cells, max abs, max rel, the worst graph and the returned Python type. By family:

| Family | Indices | Value comparisons | Both undefined (NaN + reason) | Disagreements | Max rel. diff |
|---|---|---|---|---|---|
| degree | 32 | 35,163 | 5 (AZI on K2 components) | 0 | 2.3e-15 |
| reverse_degree (R1/R2/R3 × 16, HCM2) | 49 | 51,528 | 2,320 (k > Δ) | 3 (HCM2, see below) | 1.9e-15 |
| neighbourhood_degree_sum | 10 | 10,990 | 0 | 0 | 1.6e-15 |
| distance | 11 | 12,044 | 45 (disconnected / empty / J without edges) | 0 | 2.5e-16 |
| matching (Z) | 1 | 1,099 | 0 | 0 | 1.1e-16 (table route only; direct route exact) |
| entropy (I_deg) | 1 | 1,096 | 3 (no edges) | 0 | 1.6e-15 |

- **Integer-valued indices** (M1, M2, F, HM, sigma, ReZ3, BM, TM, S_ReZ3, W, Sz, Sz_e, PI, PI_v, Mo, Mo_e, W_e, DD, Gut, Z, and the integer R_k ones) are returned as Python `int` on the direct route and agree exactly. BMH, TMH, S_TMH and R_k_BMH/TMH can be half-integers; they are returned as floats and agree exactly.
- **Rational, non-integer indices** (H, ISI, mM2, SDD, ReZ1, HBM, HTM, ...) are returned as floats and are correct to within 1 ulp-scale relative error.

### Disagreements, classified

| Index | Class | Minimal witness | Detail |
|---|---|---|---|
| HCM2 | DEFINITION DIFFERENCE (edge case) | K1 (also empty graph, 3K1) | Reference: Kulli's c_v = Δ − d_v + 1 is defined for Δ = 0 and the edge sum is empty, so HCM2 = 0. Package: NaN, reason "modified reverse degree needs Delta(G) >= 1". HCM2 seems to be computed through M_1R, which the package leaves undefined at Δ = 0. Only edgeless graphs reach this, and `standardize` rejects them (fewer than 3 heavy atoms, or disconnected). It is defensible, but it is inconsistent with the other degree-family indices, which return 0 on edgeless graphs. |
| Z (table route only) | FLOATING POINT | CHEMBL1688530 (n = 125), or the alkane P80 (Z = F_81 = 37,889,062,373,143,906) | `compute_indices` returns the exact int. `index_table` stores every index in a float64 DataFrame, so a Z above 2^53 is rounded: absolute error −2 for P80, −14,533 for P100, and up to 1,195,973,136 for CHEMBL1688530 (relative error ≤ 1.1e-16). This affects 9 of 1,821 cohort molecules. It is harmless for modelling, which uses floats anyway, but the documentation's claim that Hosoya Z is returned as an exact int holds only for `compute_indices`. Exported tables are not exact. |

There are no CODE BUG, INPUT NORMALIZATION, SOURCE TYPO or UNRESOLVED cases.

### Undefined-value handling (checked, `out/undefined_probes.csv`)

- **AZI.** NaN on every graph with a (1,1) edge (K2, 2K2), with the reason "division by zero on edge (0,1)".
- **R_k for k > Δ.** NaN with a reason, e.g. R2_\*/R3_\* on K2 and R3_\* on paths and cycles. In the cohort this is exactly the 80 R3 cells of Δ = 2 molecules. R_k for k ≤ Δ is never NaN, which confirms that M_kR is a permutation of 1..Δ.
- **Distance indices.** NaN with the reason "graph has c components" on 2K2, P3 ∪ K1 and 3K1, and "empty graph" for n = 0. J is NaN on any graph without edges. W and the rest are 0 on K1 in both implementations.
- **I_deg.** NaN on graphs without edges. With an isolated vertex (P3 ∪ K1) the package uses the 0·log 0 = 0 limit, the same as the reference.
- **No zero imputation and no NaN without a reason** in any of the 1,099 (graph, route) evaluations.

### Rings

The package's `n_rings` (from `size_table`) equals the cyclomatic number m − n + 1 for all 507 atlas/closed-form SMILES and all 1,821 cohort molecules. Each of these graphs is connected, as standardisation requires. RDKit's default `RingInfo.NumRings()` (symmetrised SSSR) differs from m − n + 1 on 37 cohort molecules, for example CHEMBL1688530 (15 vs 13). The package does not use that count, which is correct.

## Algebraic dependencies

All of these come from the exact reference values. They were found on the atlas and then confirmed on non-atlas graphs and on the package's float values for the whole cohort. Raw output: `out/dependencies.json` and `out/full_panel.json`.

1. **`ReZ1 ≡ n`.** On every graph without isolated vertices, (d_u + d_v)/(d_u d_v) = 1/d_u + 1/d_v, so the sum over edges is Σ_v d_v · (1/d_v) = n.
   - It holds exactly on all 586 connected test graphs and, to 1e-9, on all 1,821 cohort molecules.
   - `docs/INDEX_VERIFICATION.md` says the "redefined first Zagreb index Σ(1/d_u + 1/d_v) equals n ... so it is not registered". It *is* registered, as `ReZ1`, which was added 2026-10-05 for Experiment D. `ReZ1` is missing from the "Linear identities kept on purpose" list.
   - In the same way, R_k_ReZ1 = Σ_v d_v / M_kR(d_v), which depends only on the degree sequence.
2. **Functions of the degree sequence** (constant within all 83 groups of non-isomorphic graphs that share (n, m, degree sequence)): M1, F, ReZ1, R1/R2/R3_M1, R1/R2/R3_F, R1/R2/R3_ReZ1 and I_deg. Each is vertex-additive. For Δ = 4 every one of them except I_deg is an exact rational combination of n, m, M1 and F, because the degree sequence has 4 counts n_1..n_4 and a Vandermonde system ties them to n, 2m, M1 = Σd²n_d and F = Σd³n_d. Examples, all holding on all 792 Δ = 4 cohort molecules:
   - R1_M1 = 10m − M1
   - R1_F = 50m − 10M1 + F
   - R2_M1 = −16n + (140/3)m − 7M1 + (2/3)F
   - R3_ReZ1 = −(25/3)n + (505/18)m − (143/24)M1 + (55/72)F
3. **Universal linear identities**, exact on every connected test graph:
   - HM = F + 2M2
   - sigma = F − 2M2
   - BM = M1 + M2
   - TM = F + M2
   - R1_F = F − 2M2 + 2R1_M2
   - R1_HM = F − 2M2 + 4R1_M2
   - R1_TM = F − 2M2 + 3R1_M2
   - R1_BM = R1_M1 + R1_M2
   - R1_BMH = ½F − M2 + 2R1_M2 + ½R1_ReZ3

   The R1 ones follow from M_1R(d) = Δ + 1 − d. No pair of indices is exactly equal or proportional.
4. **Edge-type collapse (the main redundancy).**
   - Every degree, R_k and HCM2 index has the form Σ_{a ≤ b} N_ab φ(w(a), w(b)), where N_ab is the number of edges joining degrees a and b. For Δ ≤ 4 and no K2 component there are only 9 edge types: (1,2) through (4,4).
   - The exact coefficient matrix of the 32 degree indices has rank 9. Adding n, m, all 48 R_k and HCM2 with Δ = 4 weights still gives rank 9.
   - On the package's own values: 792 Δ = 4 cohort molecules × 81 indices have numerical rank 9 (σ_9/σ_1 = 9e-4, σ_10/σ_1 = 1.2e-16). 1,024 Δ = 3 molecules have rank 5 (5 edge types).
   - Across both Δ strata, the [1, n, m, 104 indices] matrix on the 3 ≤ Δ ≤ 4 atlas has rank 38. Pivots: const, n, m, M1, M2, F, R, RR, SC, H, R1_M1, R1_M2, R1_ReZ3, R1_TMH, R1_R; the 10 S_\* indices; the 11 distance indices; Z; I_deg.
   - All 69 dependent columns hold on every non-atlas 3 ≤ Δ ≤ 4 test graph. They also hold on the package's floats for all 1,816 such cohort molecules (worst relative residual 1.1e-15).
   - In practical terms, the TI feature block carries about 38 linearly independent directions, not 104. Only 14 of them come from the 81 degree/R_k/HCM2 indices together with n and m.
   - The coefficients of most of these relations are large irrational or rational numbers; they are listed in `dependencies.json`. Outside the molecular regime (Δ ≥ 5, or Δ ≤ 2 for R_k), these relations do not hold, and the script reports a witness, e.g. S5.
5. **Cyclomatic identity.** rings = m − n + 1 holds for all connected test graphs and cohort molecules (see Rings).

## Reproduce

From the repository root (about 3.5 minutes; the run also standardises the full 1,821-molecule cohort):

```
~/.venvs/topoqspr/bin/python -I docs/research/index_verification/verify_indices.py \
    ~/topoqspr_runs/main_confirmatory/molecules.csv docs/research/index_verification/out
```

Outputs in `out/`:

- `summary_per_index.csv`: one row per index
- `problems.csv`: every non-agreeing cell
- `closed_forms.csv`
- `undefined_probes.csv`
- `drugs.csv`: graph isomorphism, rings, Z
- `hosoya_large.csv`: Z > 2^53, direct int vs float table
- `dependencies.json`, `full_panel.json`
- `run_log.txt`

## Limits

- Values were checked against the definitions as documented. I did not re-read the primary papers, so whether the documented conventions match each primary source (for example the mod precedence in M_kR, or the W_e convention, where ChemGraphX differs by m(m−1)/2) rests on the existing `docs/INDEX_VERIFICATION.md` evidence.
- The exhaustive set stops at n = 7. Larger graphs are covered only by the families and the 36 drugs (n up to 125).
- The dependency relations with irrational coefficients were found numerically at 60 digits (tolerance 1e-35) and confirmed on held-out graphs. They were not derived symbolically, except for the edge-type rank argument, which is exact by construction.
