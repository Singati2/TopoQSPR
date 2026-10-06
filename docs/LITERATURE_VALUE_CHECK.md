# Independent literature value check (final, 2026-10-05)

`tests/test_literature_values.py` (199 test cases, fixtures in `tests/fixtures/literature/`, each
with DOI and source table) compares TopoQSPR's index values with values printed by seven papers from
groups other than the senior author's. Structures are PubChem connectivity SMILES looked up by the
printed name (the papers print drawings); formulas and partitions are quoted in the fixtures.

Outcome: **no discrepancy attributable to TopoQSPR.** Every printed value that follows from the
paper's own structure or partition reproduces at its printed precision, with these exceptions, each
listed by row in the test file:

| Paper | Checked | Result |
|---|---|---|
| Youssef et al., Sci. Rep. 14 (2024) 20533, FeTe2 network, reverse degree | Table 1 partition (64 (m, n) cases) and Tables 3–5 (reverse-degree indices from the partition) | all reproduce; partition consistent with the stated vertex counts |
| Albadrani et al., Materials 15 (2022) 2889, graphene sheets, reverse degree | Table 3 (full pipeline on named sheets; partition level on the others), Remarks 1–3 closed forms | reproduce, except: 3 printed cells contradict the paper's own partition (C32 GA, C30 SC, C30 H); the C22 row implies a non-integer vertex count (no graph); Remark 2 prints two wrong constants (RR, RReZG1) and Remark 3 misprints the first reverse Zagreb index; 4 cells differ from ours by one unit of the last printed digit (truncation) |
| Rao et al., ACS Omega 9 (2024) 31056, reverse Zagreb indices | Table 3 (full pipeline on named compounds), Theorems 2.1–2.7 on explicit cactus graphs | Table 3 reproduces; two rows print an odd M1 (impossible for any graph); Theorems 2.1, 2.3, 2.5–2.7 reproduce exactly; Theorem 2.2 reproduces except its M2 (printed 4n + 2k² + 6k; the defined graph has 4n + k² + 6k); Theorem 2.4's formulas fit no single graph (its M1 and CM1 are those of a cycle with pendant edges at k − 1 other vertices, its products and hyper indices are not) |
| Abubakar et al., Heliyon 10 (2024) e28260, antituberculosis drugs, neighbourhood degree sum | Table 2 partition (terizidone), Table 3 (full pipeline on the drugs whose printed NM1 equals 2 M2 of the real graph) | reproduce, except 2 printed cells (ethambutol NABC, linezolid NHP); 7 rows print an odd NM1 (impossible, NM1 = 2 M2) and 3 rows belong to a different graph than the named drug; 4 cells within one unit of the last printed digit |
| Almulla et al., Sci. Rep. 16 (2026) 25695, antifungal drugs, neighbourhood degree sum | Tables 2 and 4 (partition and indices, full pipeline) | 16 of 20 drugs reproduce exactly; for 4 drugs the printed partition is not that of the named molecule (two have an odd NM1), while the printed indices do follow from the printed partition |
| Ravi & Desikan, Sci. Rep. 13 (2023), reduced reverse degree | definition (RR = c + 1), Table 2 benzenoids, Table 7 HA-paclitaxel partition and closed forms | Table 2 is labelled RRR₋₁ but prints plain-degree values; the closed forms of Table 7 follow from the printed partition except GA and R, whose printed slopes (93.99, 26.18 per unit) are not those of the partition (94.46, 28.86) |
| Ravi & Chidambaram, Sci. Rep. 15 (2025) 16941, Parkinson drugs, open neighbourhood sum | Table 2 partition (procyclidine) | reproduces |
| Mondal, De, Pal, Eur. Phys. J. Plus 136 (2021) 303, Cayley trees | closed forms of neighbourhood indices | reproduce |

"Within one unit of the last printed digit" (8 cells in two papers) means |ours − printed| ≤ one
unit of the last printed decimal while the printed value is not the correctly rounded one: the
authors truncated, or rounded intermediate values. Such cells are listed explicitly in
`_LAST_DIGIT`; the test refuses a silent tolerance anywhere else.

History: an agent built the check (174 pass, 25 fail) and stopped at a usage limit on 2026-10-05;
the lead author classified the 25 failures the same day (the 18 Rao failures were the Theorem 2.2
M2 misprint and the inconsistent Theorem 2.4, verified against the article's own statements in
Europe PMC; the 7 others were last-digit truncations and one arbitrary tolerance in the test).
