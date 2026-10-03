# Topological index verification (TopoQSPR v1)

This file is the provenance record for every index in
`src/topoqspr/indices.py` (`REGISTRY`, 86 indices) and for every one of the
135 index names of v0.x (`config.py:TOPO_INDICES`). Bibliography keys refer to
`docs/index_refs.bib`. Test names refer to `tests/test_indices.py` (and
`tests/test_graph.py` for standardisation). The tests
`test_every_index_is_documented` and `test_every_legacy_name_is_mapped` fail
if a registry index or a legacy name is missing from this file.

## Conventions

* **Graph.** `topoqspr.graph.standardize` parses the SMILES with RDKit, runs
  the ChEMBL Structure Pipeline (`standardize_mol`, then `get_parent_mol`;
  Bento et al. 2020, `bento2020`), and builds the hydrogen-suppressed graph of
  the parent: heavy atoms are vertices, bonds between heavy atoms are edges.
  Bond order, charge, element and isotope are ignored. The pipeline strips
  isotope labels, so `[2H]` drugs get the graph of the protium parent.
  Input with internal whitespace is rejected as `unparsable`, because RDKit
  would otherwise read `"CCO junk"` as ethanol named "junk".
* **Parent fallback.** `get_parent_mol` keeps *every* fragment when all of
  them are on its salt list, so sodium benzoate, sodium salicylate, choline
  chloride and ammonium stearate used to come out as `disallowed_element` or
  `mixture`, while the free acid or base was accepted. When more than one
  fragment survives, v1 drops fragments without heavy atoms (`[H+]`),
  single-heavy-atom Li/Na/K/Rb/Cs/Mg/Ca/Sr/Ba/F/Cl/Br/I/N ions, and the
  counter-ions in `SUPPLEMENTARY_COUNTERIONS`. Of 24 common pharmaceutical
  counter-ions tested, xinafoate was the only one missing from the ChEMBL
  1.2.4 list. Identical copies are then merged. If one organic fragment
  remains, it is uncharged and becomes the parent, and
  `StandardizationResult.parent_note` records what was dropped. Lithium
  carbonate is not turned into carbonic acid. If several fragments remain
  and every one is on the ChEMBL salt list (choline salicylate), the input
  is rejected as `parent_ambiguous`. Tests:
  `test_salts_and_solvents_are_stripped_to_the_parent`,
  `test_choline_chloride_keeps_the_quaternary_cation`,
  `test_parent_fallback_is_recorded_and_absent_when_unused`,
  `test_lithium_carbonate_is_not_turned_into_carbonic_acid`.
* **Rejections** (`Rejection.code`): `unparsable`, `standardization_failed`,
  `disallowed_element` (default allowed set H, B, C, N, O, F, Si, P, S, Cl,
  Se, Br, I; configurable), `no_carbon`, `mixture` (more than one fragment
  after parent extraction and fallback), `parent_ambiguous` (more than one
  fragment, all on the ChEMBL salt list), `too_few_heavy_atoms` (< 3, configurable),
  `disconnected_graph`. The standardised canonical SMILES and the full
  InChIKey are recorded for accepted structures and for most rejections.
  Tests: `test_rejection_reasons`, `test_salts_and_solvents_are_stripped_to_the_parent`,
  `test_hydrogens_are_never_vertices`, `test_allowed_elements_and_min_heavy_atoms_are_configurable`.
* **Undefined values.** An index that is undefined for a graph is `NaN` and
  `compute_indices` returns the reason in its second dict. Cases: division by
  zero in an edge function (for example AZI on an edge with d_u = d_v = 1),
  a disconnected or empty graph for any distance index, no edges for the
  entropy, k outside 1..Delta(G) for the modified reverse degree (and
  Delta(G) = 0, which includes the empty graph), and Hosoya Z when the
  frontier width exceeds `HOSOYA_MAX_FRONTIER` (see the Hosoya Z decision). A real
  zero (for example sigma of a regular graph, or ABC of K2) stays 0. Nothing is
  imputed. Tests: `test_reasons_are_given_exactly_for_nan_values`,
  `test_division_by_zero_is_nan_with_reason`, `test_disconnected_graph_distance_indices_are_nan`,
  `test_rk_beyond_delta_is_nan_with_reason`, `test_edgeless_graph`,
  `test_rk_is_nan_on_graphs_without_edges_including_empty`,
  `test_hosoya_frontier_cap_gives_nan_with_reason_quickly`.
* **Exactness.** Integer-valued indices (Zagreb-type, W, Sz, PI, Mo, DD, Gut,
  Z, ...) are returned as Python `int`. Hosoya Z is exact whenever it is
  returned; above the frontier-width cap it is NaN with a reason.
* **Invariance and determinism.** `test_invariance_under_atom_renumbering`
  (random RDKit atom orders, 5 drugs, all 86 indices),
  `test_invariance_under_graph_relabelling`,
  `test_values_are_in_registry_order_and_deterministic`.
* **No duplicates.** `test_no_two_indices_are_identical_on_the_drug_panel`
  checks every pair of the 86 indices on 24 approved drugs.
* **Edge-function notation.** Edge-additive indices are
  `sum_{uv in E} phi(x_u, x_v)` with x = d (degree family), S (neighbourhood
  degree sum, `S_u = sum_{w ~ u} d_w`) or `M_kR(d)` (modified reverse degree).
  The same `phi` is shared across families (`indices.PHI`), which is exactly
  the construction of Godlin & Radha 2026, Eqs. (2)-(3), and Arockiaraj et
  al. 2025, Eqs. (2)-(4).

## Verification evidence (summary)

* Closed forms on P_n (n = 3..30), C_n (n = 3..30), stars K_{1,n-1} and K_n:
  every degree-family index has a hand-derived value on C_n and the irregular
  ones also on stars.
* Brute-force enumeration of matchings on 40 random graphs (including
  disconnected ones) for Hosoya Z.
* networkx oracle (all-pairs shortest paths, `wiener_index`, `line_graph`) for
  W, W_e, J, Sz, Sz_e, PI, PI_v, Mo, Mo_e, DD, Gut on 25 random connected
  graphs and 6 drugs.
* Mordred-community 2.0.7 on 24 approved drugs (fixture
  `tests/fixtures/indices/mordred_reference.json`, regenerated by
  `tests/fixtures/indices/make_fixtures.py`; `test_mordred_fixture_is_reproducible`
  recomputes part of it live): ABC, BalabanJ, WPath, Zagreb1, Zagreb2, Xp-1d
  agree with ABC, J, W, M1, M2, R to 1e-10.
* Published worked examples: all 40 values of the hydrochlorothiazide example
  of Godlin & Radha 2026, Table 5 (10 indices for the degree and for k = 1, 2,
  3; rounding to 2 decimals; `test_godlin_hctz_table_is_complete` checks that
  all 40 are encoded) and
  the five gemcitabine ABC values of Arockiaraj et al. 2025 (degree,
  neighbourhood, k = 1, 2, 3) are reproduced.
* Hosoya Z of 5 drugs equals the values of the independent recomputation in
  the expert audit (aspirin 335, donepezil 965,380, sitagliptin 394,051,
  imatinib 79,015,504, atorvastatin 267,895,296).

## Index table

Status `verified` means: definition read in the cited source (or in the
open-access paper named in the bib note when the primary full text was not
accessible), implementation matches it, and the listed tests pass.

| Index | Definition | Primary source | Code function | Verified by (test names) | Status |
|---|---|---|---|---|---|
| `M1` | $\sum_{uv \in E(G)} \left[d_u + d_v\right]$ | gutman1972 | `_deg_func` + `PHI` | `test_cycle_closed_forms_for_every_degree_index`, `test_path_closed_forms`, `test_star_closed_forms`, `test_complete_graph_closed_forms`, `test_overlapping_indices_match_mordred`, `test_algebraic_identities` | verified |
| `M2` | $\sum_{uv \in E(G)} d_u d_v$ | gutman1975 | `_deg_func` + `PHI` | `test_cycle_closed_forms_for_every_degree_index`, `test_path_closed_forms`, `test_star_closed_forms`, `test_complete_graph_closed_forms`, `test_overlapping_indices_match_mordred` | verified |
| `F` | $\sum_{uv \in E(G)} \left[d_u^2 + d_v^2\right]$ | furtula2015 | `_deg_func` + `PHI` | `test_cycle_closed_forms_for_every_degree_index`, `test_path_closed_forms`, `test_star_closed_forms`, `test_algebraic_identities` | verified |
| `HM` | $\sum_{uv \in E(G)} (d_u + d_v)^2$ | shirdel2013 | `_deg_func` + `PHI` | `test_cycle_closed_forms_for_every_degree_index`, `test_algebraic_identities` | verified |
| `R` | $\sum_{uv \in E(G)} \frac{1}{\sqrt{d_u d_v}}$ | randic1975 | `_deg_func` + `PHI` | `test_cycle_closed_forms_for_every_degree_index`, `test_path_closed_forms`, `test_star_closed_forms`, `test_overlapping_indices_match_mordred`, `test_registry_reproduces_godlin_2026`, `test_hydrochlorothiazide_matches_godlin_2026` | verified |
| `RR` | $\sum_{uv \in E(G)} \sqrt{d_u d_v}$ | gutman2014rr | `_deg_func` + `PHI` | `test_cycle_closed_forms_for_every_degree_index`, `test_star_closed_forms_for_irregular_indices` | verified |
| `SC` | $\sum_{uv \in E(G)} \frac{1}{\sqrt{d_u + d_v}}$ | zhou2009 | `_deg_func` + `PHI` | `test_cycle_closed_forms_for_every_degree_index`, `test_star_closed_forms_for_irregular_indices` | verified |
| `H` | $\sum_{uv \in E(G)} \frac{2}{d_u + d_v}$ | fajtlowicz1987 | `_deg_func` + `PHI` | `test_cycle_closed_forms_for_every_degree_index`, `test_star_closed_forms_for_irregular_indices` | verified |
| `GA` | $\sum_{uv \in E(G)} \frac{2\sqrt{d_u d_v}}{d_u + d_v}$ | vukicevic2009 | `_deg_func` + `PHI` | `test_cycle_closed_forms_for_every_degree_index`, `test_star_closed_forms_for_irregular_indices` | verified |
| `ABC` | $\sum_{uv \in E(G)} \sqrt{\frac{d_u + d_v - 2}{d_u d_v}}$ | estrada1998 | `_deg_func` + `PHI` | `test_cycle_closed_forms_for_every_degree_index`, `test_star_closed_forms_for_irregular_indices`, `test_overlapping_indices_match_mordred`, `test_registry_reproduces_godlin_2026`, `test_hydrochlorothiazide_matches_godlin_2026`, `test_gemcitabine_matches_arockiaraj_2025` | verified |
| `AZI` | $\sum_{uv \in E(G)} \left(\frac{d_u d_v}{d_u + d_v - 2}\right)^3$ | furtula2010 | `_deg_func` + `PHI` | `test_cycle_closed_forms_for_every_degree_index`, `test_star_closed_forms_for_irregular_indices`, `test_division_by_zero_is_nan_with_reason` | verified |
| `SO` | $\sum_{uv \in E(G)} \sqrt{d_u^2 + d_v^2}$ | gutman2021 | `_deg_func` + `PHI` | `test_cycle_closed_forms_for_every_degree_index`, `test_star_closed_forms_for_irregular_indices` | verified |
| `ISI` | $\sum_{uv \in E(G)} \frac{d_u d_v}{d_u + d_v}$ | vukicevic2010 | `_deg_func` + `PHI` | `test_cycle_closed_forms_for_every_degree_index`, `test_star_closed_forms_for_irregular_indices`, `test_registry_reproduces_godlin_2026`, `test_hydrochlorothiazide_matches_godlin_2026` (as ReZ2) | verified |
| `SDD` | $\sum_{uv \in E(G)} \frac{d_u^2 + d_v^2}{d_u d_v}$ | vukicevic2010 | `_deg_func` + `PHI` | `test_cycle_closed_forms_for_every_degree_index`, `test_registry_reproduces_godlin_2026`, `test_hydrochlorothiazide_matches_godlin_2026` (as SSD) | verified |
| `sigma` | $\sum_{uv \in E(G)} (d_u - d_v)^2$ | gutman2018sigma | `_deg_func` + `PHI` | `test_cycle_closed_forms_for_every_degree_index`, `test_cycle_closed_forms`, `test_star_closed_forms`, `test_algebraic_identities` | verified |
| `ABS` | $\sum_{uv \in E(G)} \sqrt{\frac{d_u + d_v - 2}{d_u + d_v}}$ | ali2022 | `_deg_func` + `PHI` | `test_cycle_closed_forms_for_every_degree_index`, `test_star_closed_forms_for_irregular_indices` | verified |
| `ReZ3` | $\sum_{uv \in E(G)} d_u d_v (d_u + d_v)$ | ranjini2013 | `_deg_func` + `PHI` | `test_cycle_closed_forms_for_every_degree_index`, `test_registry_reproduces_godlin_2026`, `test_hydrochlorothiazide_matches_godlin_2026` | verified |
| `AG` | $\sum_{uv \in E(G)} \frac{d_u + d_v}{2\sqrt{d_u d_v}}$ | shegehall2015 | `_deg_func` + `PHI` | `test_cycle_closed_forms_for_every_degree_index`, `test_registry_reproduces_godlin_2026`, `test_hydrochlorothiazide_matches_godlin_2026` | verified |
| `HG` | $\sum_{uv \in E(G)} \frac{2}{\sqrt{d_u d_v}\,(d_u + d_v)}$ | arockiaraj2023sar | `_deg_func` + `PHI` | `test_cycle_closed_forms_for_every_degree_index`, `test_registry_reproduces_godlin_2026`, `test_hydrochlorothiazide_matches_godlin_2026` | verified; primary full text not accessed (formula from godlin2026) |
| `GBM` | $\sum_{uv \in E(G)} \frac{\sqrt{d_u d_v}}{d_u + d_v + d_u d_v}$ | arockiaraj2023sar | `_deg_func` + `PHI` | `test_cycle_closed_forms_for_every_degree_index`, `test_registry_reproduces_godlin_2026`, `test_hydrochlorothiazide_matches_godlin_2026` | verified; primary full text not accessed (formula from godlin2026, arockiaraj2025lung) |
| `BMG` | $\sum_{uv \in E(G)} \frac{d_u + d_v + d_u d_v}{\sqrt{d_u d_v}}$ | arockiaraj2023sar | `_deg_func` + `PHI` | `test_cycle_closed_forms_for_every_degree_index`, `test_registry_reproduces_godlin_2026`, `test_hydrochlorothiazide_matches_godlin_2026` | verified; primary full text not accessed (formula from godlin2026) |
| `BMH` | $\sum_{uv \in E(G)} \frac{(d_u + d_v + d_u d_v)(d_u + d_v)}{2}$ | arockiaraj2023sar | `_deg_func` + `PHI` | `test_cycle_closed_forms_for_every_degree_index` | verified; primary full text not accessed (formula from kalaam2024, Eq. 8) |
| `TMH` | $\sum_{uv \in E(G)} \frac{(d_u^2 + d_v^2 + d_u d_v)(d_u + d_v)}{2}$ | arockiaraj2023sar | `_deg_func` + `PHI` | `test_cycle_closed_forms_for_every_degree_index` | verified; primary full text not accessed (formula from arockiaraj2025lung, kalaam2024 Eq. 9) |
| `HCM2` | $\sum_{uv \in E(G)} (c_u c_v)^2,\ c_u = \Delta - d_u + 1$ | kulli2018 | `_rk_func` + `modified_reverse_degree` + `PHI` | `test_star_closed_forms`, `test_cycle_closed_forms_for_every_degree_index` | verified |
| `S_ABC` | $\sum_{uv \in E(G)} \sqrt{\frac{S_u + S_v - 2}{S_u S_v}}$ | ghorbani2010 | `_nbr_func` + `PHI` | `test_cycle_closed_forms_for_every_degree_index`, `test_gemcitabine_matches_arockiaraj_2025`, `test_legacy_defects_are_fixed` | verified (v0.x numerator fixed) |
| `S_GA` | $\sum_{uv \in E(G)} \frac{2\sqrt{S_u S_v}}{S_u + S_v}$ | graovac2011 | `_nbr_func` + `PHI` | `test_cycle_closed_forms_for_every_degree_index`, `test_star_closed_forms_for_irregular_indices` | verified |
| `S_RR` | $\sum_{uv \in E(G)} \sqrt{S_u S_v}$ | mondal2021 | `_nbr_func` + `PHI` | `test_cycle_closed_forms_for_every_degree_index`, `test_star_closed_forms_for_irregular_indices` | verified |
| `S_SC` | $\sum_{uv \in E(G)} \frac{1}{\sqrt{S_u + S_v}}$ | mondal2021 | `_nbr_func` + `PHI` | `test_cycle_closed_forms_for_every_degree_index` | verified |
| `S_ReZ3` | $\sum_{uv \in E(G)} S_u S_v (S_u + S_v)$ | mondal2021 | `_nbr_func` + `PHI` | `test_cycle_closed_forms_for_every_degree_index` | verified |
| `S_R` | $\sum_{uv \in E(G)} \frac{1}{\sqrt{S_u S_v}}$ | mondal2021 | `_nbr_func` + `PHI` | `test_cycle_closed_forms_for_every_degree_index`, `test_neighbourhood_sum_zero_is_nan` | verified |
| `S_SDD` | $\sum_{uv \in E(G)} \frac{S_u^2 + S_v^2}{S_u S_v}$ | mondal2021 | `_nbr_func` + `PHI` | `test_cycle_closed_forms_for_every_degree_index` | verified |
| `S_H` | $\sum_{uv \in E(G)} \frac{2}{S_u + S_v}$ | arockiaraj2025lung | `_nbr_func` + `PHI` | `test_cycle_closed_forms_for_every_degree_index` | verified |
| `S_GBM` | $\sum_{uv \in E(G)} \frac{\sqrt{S_u S_v}}{S_u + S_v + S_u S_v}$ | arockiaraj2025lung | `_nbr_func` + `PHI` | `test_cycle_closed_forms_for_every_degree_index` | verified |
| `S_TMH` | $\sum_{uv \in E(G)} \frac{(S_u^2 + S_v^2 + S_u S_v)(S_u + S_v)}{2}$ | arockiaraj2025lung | `_nbr_func` + `PHI` | `test_cycle_closed_forms_for_every_degree_index` | verified |
| `W` | $\sum_{\{u,v\} \subseteq V} D(u,v)$ | wiener1947 | `_wiener` | `test_path_closed_forms`, `test_cycle_closed_forms`, `test_star_closed_forms`, `test_complete_graph_closed_forms`, `test_distance_indices_match_networkx_on_random_connected_graphs`, `test_distance_indices_match_networkx_on_drugs`, `test_overlapping_indices_match_mordred` | verified |
| `J` | $\frac{m}{\mu + 1} \sum_{uv \in E} (D_u D_v)^{-1/2},\ \mu = m - n + 1$ | balaban1982 | `_balaban_j` | `test_path_closed_forms`, `test_cycle_closed_forms`, `test_distance_indices_match_networkx_on_random_connected_graphs`, `test_distance_indices_match_networkx_on_drugs`, `test_overlapping_indices_match_mordred` | verified (connected graphs only) |
| `Sz` | $\sum_{uv \in E} n_u(e) n_v(e)$ | gutman1994sz | `_szeged` | `test_path_closed_forms` (Sz = W on trees), `test_cycle_closed_forms`, `test_complete_graph_closed_forms`, `test_distance_indices_match_networkx_on_random_connected_graphs`, `test_distance_indices_match_networkx_on_drugs` | verified |
| `Sz_e` | $\sum_{uv \in E} m_u(e) m_v(e)$ | gutman2008sze | `_edge_szeged` | `test_distance_indices_match_networkx_on_random_connected_graphs`, `test_distance_indices_match_networkx_on_drugs` | verified |
| `PI` | $\sum_{uv \in E} [m_u(e) + m_v(e)]$ | khadikar2000 | `_pi_edge` | `test_cycle_closed_forms`, `test_distance_indices_match_networkx_on_random_connected_graphs`, `test_distance_indices_match_networkx_on_drugs` | verified |
| `PI_v` | $\sum_{uv \in E} [n_u(e) + n_v(e)]$ | khalifeh2008pi | `_pi_vertex` | `test_path_closed_forms` (n m on bipartite graphs), `test_distance_indices_match_networkx_on_random_connected_graphs`, `test_distance_indices_match_networkx_on_drugs` | verified |
| `Mo` | $\sum_{uv \in E} \|n_u(e) - n_v(e)\|$ | doslic2018 | `_mostar` | `test_path_closed_forms`, `test_cycle_closed_forms`, `test_complete_graph_closed_forms`, `test_distance_indices_match_networkx_on_random_connected_graphs`, `test_distance_indices_match_networkx_on_drugs` | verified |
| `Mo_e` | $\sum_{uv \in E} \|m_u(e) - m_v(e)\|$ | arockiaraj2019 | `_edge_mostar` | `test_cycle_closed_forms`, `test_distance_indices_match_networkx_on_random_connected_graphs`, `test_distance_indices_match_networkx_on_drugs` | verified |
| `W_e` | $W(L(G)) = \sum_{\{e,f\} \subseteq E} D_{L(G)}(e,f)$ | iranmanesh2009 | `_edge_wiener` | `test_path_closed_forms` (L(P_n) = P_{n-1}), `test_distance_indices_match_networkx_on_random_connected_graphs`, `test_distance_indices_match_networkx_on_drugs` (networkx `line_graph`) | verified |
| `DD` | $\sum_{\{u,v\} \subseteq V} (d_u + d_v) D(u,v)$ | dobrynin1994 | `_degree_distance` | `test_path_closed_forms` (DD = 4W - n(n-1) on trees), `test_distance_indices_match_networkx_on_random_connected_graphs`, `test_distance_indices_match_networkx_on_drugs` | verified |
| `Gut` | $\sum_{\{u,v\} \subseteq V} d_u d_v D(u,v)$ | gutman1994schultz | `_gutman` | `test_path_closed_forms` (Gut = 4W - (2n-1)(n-1) on trees), `test_distance_indices_match_networkx_on_random_connected_graphs`, `test_distance_indices_match_networkx_on_drugs` | verified |
| `Z` | $Z(G) = \sum_{k \ge 0} p(G, k)$ | hosoya1971 | `hosoya_z` | `test_path_closed_forms` (F_{n+1}), `test_cycle_closed_forms` (L_n), `test_star_closed_forms`, `test_complete_graph_closed_forms` (involutions), `test_hosoya_matches_brute_force_on_random_graphs`, `test_benzene_and_naphthalene_hosoya`, `test_hosoya_reference_values_for_drugs`, `test_hosoya_is_exact_beyond_float_precision`, `test_hosoya_speed_on_drug_size_graphs` | verified (v0.x approximation replaced by exact count) |
| `I_deg` | $-\sum_{v \in V} \frac{d_v}{2m} \log_2 \frac{d_v}{2m}$ | cao2014 | `_degree_entropy` | `test_cycle_closed_forms`, `test_complete_graph_closed_forms`, `test_star_closed_forms_for_irregular_indices`, `test_edgeless_graph` | verified |

### Modified reverse degree indices (family `reverse_degree`)

Registered for k = 1, 2, 3 (`RK_VALUES`). Name pattern `R{k}_{phi}`, value
`sum_{uv in E} phi(M_kR(d_u), M_kR(d_v))`. Code: `_rk_func` +
`modified_reverse_degree` + `PHI`. All are verified by
`test_rk_indices_use_the_published_weight_tables` (weights looked up in the
printed Delta = 4 tables of the sources, edge functions re-typed in the test,
three Delta = 4 molecules), `test_rk_beyond_delta_is_nan_with_reason` and
`test_drugs_have_no_undefined_indices_except_rk_beyond_delta`; further tests
per row.

| Index | phi(a, b) | Primary source | Further tests | Status |
|---|---|---|---|---|
| `R1_M1`, `R2_M1`, `R3_M1` | a + b | kalaam2024 Eq. (1); R1_M1 = Kulli CM1 (kulli2018) | `test_star_closed_forms`, `test_cycle_closed_forms_for_every_degree_index`, `test_algebraic_identities` | verified |
| `R1_M2`, `R2_M2`, `R3_M2` | a b | kalaam2024 Eq. (2); R1_M2 = Kulli CM2 (kulli2018) | `test_star_closed_forms`, `test_complete_graph_closed_forms`, `test_cycle_closed_forms_for_every_degree_index` | verified |
| `R1_F`, `R2_F`, `R3_F` | a^2 + b^2 | kalaam2024 Eq. (3) | `test_algebraic_identities` | verified |
| `R1_HM`, `R2_HM`, `R3_HM` | (a + b)^2 ("HZ" in the source) | kalaam2024 Eq. (4); R1_HM = Kulli HCM1 | `test_algebraic_identities` | verified |
| `R1_ReZ3`, `R2_ReZ3`, `R3_ReZ3` | (a + b) a b | kalaam2024 Eq. (5) | `test_registry_reproduces_godlin_2026` (462, 934, 738) | verified |
| `R1_BM`, `R2_BM`, `R3_BM` | a + b + a b | kalaam2024 Eq. (6) | `test_algebraic_identities` | verified |
| `R1_TM`, `R2_TM`, `R3_TM` | a^2 + b^2 + a b | kalaam2024 Eq. (7) (labelled "BM" there, a typo; the text names it tri-Zagreb) | `test_algebraic_identities` | verified |
| `R1_BMH`, `R2_BMH`, `R3_BMH` | (a + b + a b)(a + b)/2 | kalaam2024 Eq. (8) | weight-table test only | verified |
| `R1_TMH`, `R2_TMH`, `R3_TMH` | (a^2 + b^2 + a b)(a + b)/2 | kalaam2024 Eq. (9) | weight-table test only | verified |
| `R1_R`, `R2_R`, `R3_R` | (a b)^(-1/2) | godlin2026 | `test_registry_reproduces_godlin_2026` (8.55, 7.89, 8.55) | verified |
| `R1_ABC`, `R2_ABC`, `R3_ABC` | sqrt((a + b - 2)/(a b)) | godlin2026; arockiaraj2025lung | `test_registry_reproduces_godlin_2026`, `test_gemcitabine_matches_arockiaraj_2025` | verified |
| `R1_ISI`, `R2_ISI`, `R3_ISI` | a b/(a + b) ("ReZ_2" in the source) | godlin2026 | `test_hydrochlorothiazide_matches_godlin_2026` (18.42, 23.39, 20.51), `test_cycle_closed_forms_for_every_degree_index` | verified |
| `R1_RR`, `R2_RR`, `R3_RR` | sqrt(a b) ("G" in the source) | arockiaraj2025lung | weight-table test only | verified |
| `HCM2` | (c_u c_v)^2 with c = Delta - d + 1 | kulli2018 Eq. (2) | `test_star_closed_forms`, `test_cycle_closed_forms_for_every_degree_index` | verified |

## Decision: reverse degree

v0.x used `rd(v) = n + 1 - d(v)` (`ml.py`, block 3; n = number of atoms). I
found no source that defines this. With it, RM1, RM2 and RBM are polynomials
in n, m, M1 and M2 (for example `RM1 = n(n+1)^2 - 4m(n+1) + M1`), so they
mostly measure size.

Kulli 2018 (`kulli2018`, full text read), Section 1:

> "The reverse vertex degree of a vertex v in G is defined as
> c_v = Delta(G) - d_G(v) + 1."

> "The first reverse Zagreb beta index and second reverse Zagreb index [3] of
> a graph G are respectively defined as CM_1(G) = sum_{uv in E(G)} (c_u + c_v),
> CM_2(G) = sum_{uv in E(G)} c_u c_v."

> "We now introduce the first and second reverse hyper-Zagreb indices of a
> graph G as HCM_1(G) = sum_{uv in E(G)} (c_u + c_v)^2,
> HCM_2(G) = sum_{uv in E(G)} (c_u c_v)^2."

Kulli's reference [3] is S. Ediz (TWMS J. App. Eng. Math.); Kalaam et al. 2024
cite Ediz 2015 (`ediz2015`) for "the reverse version of the degree". Neither
Ediz paper was accessible.

**Decision.** v1 uses c_v = Delta - d_v + 1 and the edge-sum forms. Because
c_v = M_1R(d_v) for every vertex (k = 1 <= d_v always), Kulli's CM1, CM2 and
HCM1 are exactly `R1_M1`, `R1_M2` and `R1_HM`; they are not registered a
second time. HCM2 has no R_k counterpart in the R_k sources and is registered
as `HCM2`. The v0.x vertex forms RM1 (sum c^2) and RF (sum c^3) are dropped:
Kulli's indices are edge sums.

## Decision: generalised / modified reverse degree R_k

The generalised reverse degree was proposed by Arockiaraj, Greeni and Kalaam,
Int. J. Quantum Chem. 123 (2023) e27136 (`arockiaraj2023ijqc`). Its abstract
(Crossref) says "we propose a generalized form of reverse degree parameters";
the full text is paywalled and was not read. The definition is restated by
the same group in two open-access papers, which I read in full:

Kalaam, Greeni, Arockiaraj, Front. Chem. 12 (2024) 1470231 (`kalaam2024`),
Section 2:

> "The reverse version of the degree (Ediz, 2015) and its generalized form
> (Arockiaraj et al., 2023b) have received significant attention in recent
> years. It is denoted as M_kR d_v, incorporating a variable parameter k
> (k >= 1) and defined as
> M_kR d_v = Delta(G) - d_v + k : k <= d_v ;
> Delta(G) - d_v + k mod Delta(G) : k > d_v"

and, for the frameworks with maximum degree 4 (Section 3):

> "M_1R d_v = 4 : d_v = 1, 3 : d_v = 2, 2 : d_v = 3, 1 : d_v = 4;
> M_2R d_v = 1 : d_v = 1, 4 : d_v = 2, 3 : d_v = 3, 2 : d_v = 4;
> M_3R d_v = 2 : d_v = 1, 1 : d_v = 2, 4 : d_v = 3, 3 : d_v = 4"

Godlin & Radha, Sci. Rep. 16 (2026) 8889 (`godlin2026`), Eq. (1), writes the
parentheses explicitly:

> "M_l R(d(w_i)) = Delta(G) + l - d(w_i), l <= d(w_i);
> (Delta(G) + l - d(w_i)) mod Delta(G), l > d(w_i)"

and on the range of the parameter:

> "The maximum degrees of the chemical graphs in Figure 1 are 3 or 4;
> therefore, the modified reverse degree index phi^{m_l d} is limited to l
> from 1 to 3"

**Decisions.**

1. `mod` applies to the whole sum `(Delta + k - d) mod Delta`. This is the
   only reading consistent with the printed tables (e.g. Delta = 4, k = 2,
   d = 1: (4 + 2 - 1) mod 4 = 1). Test:
   `test_modified_reverse_degree_tables_from_the_sources`.
2. k is restricted to 1..Delta(G). For 1 <= k <= Delta the map is a
   permutation of {1, ..., Delta}: if k <= d the value is Delta + k - d in
   [k, Delta]; if k > d then Delta < Delta + k - d <= 2 Delta - 1, so the
   residue is in [1, Delta - 1]. It is never 0, so no edge function divides
   by zero and no edge is ever dropped. Test:
   `test_modified_reverse_degree_is_a_permutation_for_k_up_to_delta`.
3. For k > Delta the rule sends every vertex with d = k - Delta to 0, which
   is not a reverse degree, and the sources never use such k (Godlin & Radha
   limit l to 1..3 because their minimum Delta is 3). v1 therefore returns
   `NaN` with the reason "modified reverse degree with k=... is outside
   1..Delta(G)=..." instead of guessing. v0.x silently skipped the
   zero-weight edges: for metformin (Delta = 3) the legacy code returns
   R4_A = 4.0 from 3 of its 8 edges (re-run 2026-10-02); the expert audit
   counts 55% of approved drugs affected at k = 4 (not re-counted here).
4. Registered k values are 1, 2, 3 (the range used in both open sources). For
   drugs with Delta = 2 (rare, chain-like), `R3_*` is NaN with a reason.
5. Caveat for users, not a code defect: the same k is a different relabelling
   of degrees for Delta = 3 and Delta = 4 molecules, so R_k values for k >= 2
   are not comparable across molecules with different Delta. This is a
   property of the published definition.
6. Edge functions: the nine of Kalaam et al. 2024, Eqs. (1)-(9), plus R, ABC
   and ReZ_2 (= ISI) used with M_lR by Godlin & Radha 2026 and G (= RR) used
   by Arockiaraj et al. 2025. All 40 hydrochlorothiazide values and the 5
   gemcitabine values printed in those papers are reproduced.

## Decision: neighbourhood degree sum S(v)

Mondal, De, Pal (arXiv:1906.06660; journal version Mondal, Dey, De, Pal,
Complex Intell. Syst. 7 (2021) 977-996, `mondal2021`):

> "where delta_G(u) is the sum of degrees of all neighboring vertices of
> u in V(G), i.e, delta_G(u) = sum_{v in N_G(u)} d_G(v)"

with ND_1 = sum sqrt(delta_u delta_v), ND_2 = sum 1/sqrt(delta_u + delta_v),
ND_3 = sum delta_u delta_v (delta_u + delta_v), ND_4 = sum 1/sqrt(delta_u
delta_v), ND_5 = sum (delta_u/delta_v + delta_v/delta_u). These are
`S_RR`, `S_SC`, `S_ReZ3`, `S_R`, `S_SDD`. ND_6 is printed as
"[d_G(u) delta_G(u) + d_G(v) + delta_G(v)]" in the preprint, which is
ambiguous (likely a typo for d_u delta_u + d_v delta_v); it is not
registered. ABC_4 and GA_5 (`S_ABC`, `S_GA`) are older S-based indices
(`ghorbani2010`, `graovac2011`; definitions as quoted in Molecules 21 (2016)
821, Eqs. (4)-(5)): ABC_4 = sum sqrt((S_u + S_v - 2)/(S_u S_v)),
GA_5 = sum 2 sqrt(S_u S_v)/(S_u + S_v). `S_H`, `S_GBM`, `S_TMH` follow the
generic construction zeta^{nd}(G) = sum nd_G(r, s) zeta(r, s) of Arockiaraj et
al. 2025, Eq. (3), with zeta in {H, GBM, TMH} from their list. Not registered:
sum (S_u + S_v), which equals 2 M2 for every graph (each vertex v contributes
d_v S_v and sum_v d_v S_v = 2 sum_{uv} d_u d_v).

## Decision: Hosoya Z

Definition (Hosoya 1971): Z(G) = sum_k p(G, k), the number of sets of k
pairwise disjoint edges, with p(G, 0) = 1. v0.x enumerated matchings only for
m <= 18, returned 1 for 19 <= m <= 25 and `1 + m + m(m-1)//4` above 25 (I
re-ran the legacy code: imatinib 452, true value 79,015,504; eicosane 1, true
F_21 = 10,946). v1 counts matchings exactly with a frontier dynamic programme
over a greedy low-width vertex order (state: free processed vertices that still
have unprocessed neighbours). It is exact for every graph and returns a Python
int. The cost grows as 2^w for frontier width w of the greedy order, and
the width is not bounded by the atom count. Timings on this machine:
vancomycin (101 heavy atoms, w = 6) 0.02 s, cyclosporine (85) 0.01 s, a
148-atom peptide 0.01 s, a 150-vertex 37-ring ladder 0.01 s, C60 (w = 11)
0.2 s, a 160-vertex hexagonal sheet (w = 10) 0.8 s, a 10x15 grid (w = 11)
1.4 s, a 12x12 grid (w = 13) 4.7 s, and random 3-regular graphs on 100
vertices (w = 16-18) 8-17 s. An earlier "w <= ~8 for drugs" claim was
untested and is withdrawn (C60 alone is 11). `hosoya_z` therefore computes
w before counting. If w > `HOSOYA_MAX_FRONTIER` (= 12, configurable via
`max_frontier`; `None` disables the cap), it raises `UndefinedIndex`, so Z
is NaN with the reason "frontier width w exceeds the cap" instead of
stalling. Tests: `test_hosoya_frontier_cap_gives_nan_with_reason_quickly`,
`test_hosoya_cap_is_configurable_and_exact_when_lifted`,
`test_hosoya_default_cap_admits_c60` (Z(C60) = 1,417,036,634,543,488).
Hosoya Z is a matching count, not a distance index, so its family is
`matching`.

## Decision: graph entropy

v0.x "GE" is `-sum_v (d_v/2m) log2(d_v/2m)`, a degree-based entropy that
config.py listed under distance-based indices. Source: Cao, Dehmer, Shi 2014
(`cao2014`), the degree-power entropy `I_d^t` with t = 1, built on Dehmer's
information-functional entropy (`dehmer2008`). The full text of Cao et al.
was not accessed; the definition and the base-2 logarithm are as restated in
arXiv:1709.09594, Definition 1 ("Let G be a connected graph ... the t-th graph
entropy based on degrees of G"; "logarithms are always taken to the base
two"). v1 name `I_deg`, family `entropy`. Checked by `log2 n` on regular
graphs and `1/2 + (1/2) log2(2k)` on stars.

## Decision: Balaban J

J = m/(mu + 1) sum_{uv in E} (D_u D_v)^(-1/2) with cyclomatic number
mu = m - n + 1 (Balaban 1982). v0.x used `m - n + 2` for mu + 1, which is
right only for connected graphs, and it did not strip salts. Re-running the
legacy code: metformin hydrochloride gives J = 0 (v1, on the parent:
3.3766) and amlodipine besylate gives 9.41 (v1, on the parent: 2.3787). The
expert audit reports 66 multi-fragment approved drugs with J = 0 (not
re-counted here). v1 rejects mixtures in `standardize` and returns NaN with a
reason for any disconnected graph.

## Linear identities kept on purpose

These indices are exact linear combinations of other registered indices. They
are kept because each has its own source, but they add no information to a
linear model. They are asserted in `test_algebraic_identities`:
`sigma = F - 2 M2`, `HM = F + 2 M2`, `R{k}_HM = R{k}_F + 2 R{k}_M2`,
`R{k}_BM = R{k}_M1 + R{k}_M2`, `R{k}_TM = R{k}_F + R{k}_M2`,
`R1_M1 = sum_v d_v (Delta - d_v + 1) = 2m(Delta + 1) - M1`. Also, on a
bipartite graph `PI_v = n m`, and on trees `Sz = W`, `DD = 4W - n(n-1)`,
`Gut = 4W - (2n-1)(n-1)` (tested on paths). The redefined first Zagreb index
sum (1/d_u + 1/d_v) equals n for every graph without isolated vertices, so it
is not registered.

## v0.x defects confirmed by running the legacy code

Run on 2026-10-02 with `ml.compute_topological_indices` from the repository
root (aspirin, imatinib, n-eicosane):

| Defect | Evidence | v1 |
|---|---|---|
| Hosoya Z approximate for m > 18 | imatinib 452 (exact 79,015,504); eicosane P_20 Z = 1 (exact 10,946) | exact `hosoya_z`; `test_legacy_defects_are_fixed` |
| SS_ABC numerator S_u + S_v + S_u S_v - 2 | aspirin 15.02 vs ABC_4 = 7.51 | `S_ABC`; `test_legacy_defects_are_fixed` |
| SS_AZ, SS_ABS same extra term | code reading (`ml.py` lines with `xs + xy - 2`) | dropped (no source for the S versions) |
| RABC numerator c_u + c_v + c_u c_v - 2 and c = n + 1 - d | code reading | `R1_ABC` |
| DS1 == DS2 | both 284.0 for aspirin, 932.0 for imatinib | DS1 -> `HM`, DS2 dropped |
| DSR, DSH, DSGA identical to SC, H, GA | code reading: same formulas in d_u + d_v | dropped |
| W_v == W, Sz_v == Sz | code reading (`w_v = float(w)`, `sz_v = float(sz)`) | dropped |
| failures filled with 0 (`fillna(0.0)`) | `ml.py` end of `compute_topological_indices` | NaN + reason |
| salts not stripped before graph construction | `ml.py` only calls `RemoveHs`; metformin HCl J = 0 | ChEMBL parent + `mixture` rejection |
| R4_* drops zero-weight edges silently | metformin R4_A = 4.0 (3 of 8 edges kept) | k > Delta is NaN with reason |

## Not registered (pending a checkable source)

The twelve "hybrid geometric-harmonic-Zagreb" descriptors of Arockiaraj et al.
2023 (`arockiaraj2023sar`) include, judging by the v0.x names, GH, GTM, TMG,
HBM, HTM, BMA, TMA as well as HG, GBM, BMG, BMH, TMH. Only the last five could
be checked (their formulas are printed in open-access papers of the same
group and reproduce published values); the other seven stay out until their
formulas are read in the primary paper. The vertex-edge Wiener and Szeged
indices (W_ve, Sz_ve) also stay out: more than one convention exists and the
v0.x choice could not be tied to a source.

## Legacy name mapping

All 135 names of v0.x `config.py:TOPO_INDICES`, in their original order.
51 map to a v1 index (some renamed); 84 are dropped with the reason given.

| v0.x name | v1 name or decision |
|---|---|
| `M1` | `M1` |
| `M2` | `M2` |
| `ABC` | `ABC` |
| `R` | `R` |
| `H` | `H` |
| `F` | `F` (v0.x vertex form sum d^3 equals the edge form sum d_u^2 + d_v^2) |
| `AZI` | `AZI` |
| `GA` | `GA` |
| `SC` | `SC` |
| `BM` | dropped: BM = M1 + M2 identically; no primary source located for the plain-degree bi-Zagreb index (the modified-reverse versions `R{k}_BM` are kept, kalaam2024) |
| `TM` | dropped: TM = F + M2 identically; no primary source located for the plain-degree tri-Zagreb index (`R{k}_TM` kept, kalaam2024) |
| `GH` | dropped: one of the twelve hybrid descriptors of arockiaraj2023sar (per its abstract); full text not accessible and no open source restates the formula, so the definition could not be checked |
| `GBM` | `GBM` |
| `GTM` | dropped: one of the twelve hybrid descriptors of arockiaraj2023sar (per its abstract); full text not accessible and no open source restates the formula, so the definition could not be checked |
| `HG` | `HG` |
| `BMG` | `BMG` |
| `BMH` | `BMH` |
| `TMG` | dropped: one of the twelve hybrid descriptors of arockiaraj2023sar (per its abstract); full text not accessible and no open source restates the formula, so the definition could not be checked |
| `TMH` | `TMH` |
| `SDD` | `SDD` |
| `RM1` | dropped: vertex sum of (n+1-d_v)^2; no source defines it, and it is the size polynomial n(n+1)^2 - 4m(n+1) + M1. Kulli's first reverse Zagreb index is the edge sum `R1_M1` |
| `RM2` | `R1_M2` (Kulli CM2, now with c_v = Delta - d_v + 1) |
| `RABC` | `R1_ABC` (c_v = Delta - d_v + 1; v0.x numerator c_u + c_v + c_u c_v - 2 fixed to c_u + c_v - 2) |
| `RR` | `R1_R` (now with c_v = Delta - d_v + 1) |
| `RH` | dropped: no primary source found for a reverse harmonic index that could be checked |
| `RF` | dropped: vertex sum of (n+1-d_v)^3; no source. The edge form on Kulli degrees is `R1_F` |
| `RGA` | dropped: no primary source found for a reverse GA index that could be checked |
| `RBM` | `R1_BM` (now with c_v = Delta - d_v + 1) |
| `RSDD` | dropped: no primary source found for a reverse SDD index that could be checked |
| `DS1` | `HM` (sum (d_u + d_v)^2 is the hyper-Zagreb index) |
| `DS2` | dropped: exact duplicate of DS1 (v0.x added s*s to both) |
| `DSR` | dropped: exact duplicate of `SC` |
| `DSH` | dropped: exact duplicate of `H` |
| `DSGA` | dropped: exact duplicate of `GA` |
| `W` | `W` |
| `J` | `J` (v0.x was wrong for multi-fragment input; v1 rejects mixtures before computing) |
| `Z` | `Z` (v0.x returned 1 for 19-25 bonds and 1 + m + m(m-1)/4 above 25; v1 is exact) |
| `Sz` | `Sz` |
| `GE` | `I_deg` (degree entropy; moved from the distance family to the entropy family) |
| `W_v` | dropped: exact duplicate of `W` |
| `W_e` | `W_e` |
| `W_ve` | dropped: the vertex-edge Wiener index has more than one definition in the literature (min- and max-distance variants); the intended one could not be pinned to a primary source |
| `Sz_v` | dropped: exact duplicate of `Sz` |
| `Sz_e` | `Sz_e` |
| `Sz_ve` | dropped: as for W_ve, the vertex-edge Szeged definition could not be pinned to a primary source |
| `Mo_v` | `Mo` |
| `Mo_e` | `Mo_e` |
| `PI` | `PI` |
| `Schultz` | `DD` (degree distance, dobrynin1994; Gutman 1994 calls it the Schultz index) |
| `Gutman` | `Gut` |
| `A` | dropped: A = M1/2 exactly |
| `G` | `RR` (reciprocal Randic index, gutman2014rr) |
| `HA` | dropped: no primary source found |
| `SO` | `SO` |
| `ABC_SC` | dropped: sum sqrt((s-2)/p)/sqrt(s) is not the ABS index and has no source; the real ABS index is added as `ABS` |
| `ISI` | `ISI` |
| `sigma` | `sigma` |
| `HBM` | dropped: one of the twelve hybrid descriptors of arockiaraj2023sar (per its abstract); full text not accessible and no open source restates the formula, so the definition could not be checked |
| `HTM` | dropped: one of the twelve hybrid descriptors of arockiaraj2023sar (per its abstract); full text not accessible and no open source restates the formula, so the definition could not be checked |
| `BMA` | dropped: BMA = 2m + 2 ISI identically; formula not verifiable (see GH) |
| `TMA` | dropped: TMA = 2 M1 - 2 ISI identically; formula not verifiable (see GH) |
| `SS_M1` | dropped: sum (S_u + S_v) = 2 M2 identically |
| `SS_M2` | dropped: sum S_u S_v (fifth M2-Zagreb index); primary source not verified |
| `SS_BM` | dropped: = SS_M1 + SS_M2; no source |
| `SS_TM` | dropped: no primary source found |
| `SS_SC` | `S_SC` (ND_2, mondal2021) |
| `SS_GH` | dropped: no primary source found |
| `SS_R` | `S_R` (ND_4, mondal2021) |
| `SS_GBM` | `S_GBM` (arockiaraj2025lung) |
| `SS_A` | dropped: sum (S_u + S_v)/2 = M2 exactly |
| `SS_G` | `S_RR` (ND_1, mondal2021) |
| `SS_GA` | `S_GA` (GA_5, graovac2011) |
| `SS_H` | `S_H` (arockiaraj2025lung) |
| `SS_HG` | dropped: no primary source found for the neighbourhood version |
| `SS_HM` | dropped: no primary source found for the neighbourhood version |
| `SS_HBM` | dropped: no primary source found |
| `SS_HTM` | dropped: no primary source found |
| `SS_SDD` | `S_SDD` (ND_5, mondal2021) |
| `SS_HA` | dropped: no primary source found |
| `SS_SO` | dropped: no primary source found that could be checked for the neighbourhood Sombor index |
| `SS_BMG` | dropped: no primary source found for the neighbourhood version |
| `SS_ABC` | `S_ABC` (ABC_4, ghorbani2010; v0.x numerator S_u + S_v + S_u S_v - 2 fixed to S_u + S_v - 2) |
| `SS_BMH` | dropped: no primary source found for the neighbourhood version |
| `SS_AZ` | dropped: v0.x numerator was wrong (extra S_u S_v term) and no source for a neighbourhood AZI could be checked |
| `SS_BMA` | dropped: = 2m + 2 SS_ISI identically; no source |
| `SS_ISI` | dropped: no primary source found for the neighbourhood version |
| `SS_TMH` | `S_TMH` (arockiaraj2025lung) |
| `SS_ABS` | dropped: v0.x numerator was wrong (extra S_u S_v term) and no source for a neighbourhood ABS index could be checked |
| `SS_TMA` | dropped: no primary source found |
| `SS_sigma` | dropped: = SS_HM - 4 SS_M2 identically; no source |
| `SS_TMG` | dropped: no primary source found |
| `R1_A` | dropped: R1_A = R1_M1 / 2 exactly |
| `R1_G` | `R1_RR` |
| `R1_HA` | dropped: no primary source found |
| `R1_SO` | dropped: no source applies the Sombor function to modified reverse degrees |
| `R1_ABC_SC` | dropped: ABC_SC has no source |
| `R1_ISI` | `R1_ISI` (ReZ_2 in godlin2026) |
| `R1_sigma` | dropped: R1_sigma equals `sigma` identically, since (c_u - c_v)^2 = (d_u - d_v)^2 |
| `R1_HBM` | dropped: one of the twelve hybrid descriptors of arockiaraj2023sar (per its abstract); full text not accessible and no open source restates the formula, so the definition could not be checked |
| `R1_HTM` | dropped: one of the twelve hybrid descriptors of arockiaraj2023sar (per its abstract); full text not accessible and no open source restates the formula, so the definition could not be checked |
| `R1_BMA` | dropped: = 2m + 2 R1_ISI identically; formula not verifiable |
| `R1_TMA` | dropped: = 2 R1_M1 - 2 R1_ISI identically; formula not verifiable |
| `R2_A` | dropped: R2_A = R2_M1 / 2 exactly |
| `R2_G` | `R2_RR` |
| `R2_HA` | dropped: no primary source found |
| `R2_SO` | dropped: no source applies the Sombor function to modified reverse degrees |
| `R2_ABC_SC` | dropped: ABC_SC has no source |
| `R2_ISI` | `R2_ISI` (ReZ_2 in godlin2026) |
| `R2_sigma` | dropped: no source applies sigma to modified reverse degrees |
| `R2_HBM` | dropped: one of the twelve hybrid descriptors of arockiaraj2023sar (per its abstract); full text not accessible and no open source restates the formula, so the definition could not be checked |
| `R2_HTM` | dropped: one of the twelve hybrid descriptors of arockiaraj2023sar (per its abstract); full text not accessible and no open source restates the formula, so the definition could not be checked |
| `R2_BMA` | dropped: = 2m + 2 R2_ISI identically; formula not verifiable |
| `R2_TMA` | dropped: = 2 R2_M1 - 2 R2_ISI identically; formula not verifiable |
| `R3_A` | dropped: R3_A = R3_M1 / 2 exactly |
| `R3_G` | `R3_RR` |
| `R3_HA` | dropped: no primary source found |
| `R3_SO` | dropped: no source applies the Sombor function to modified reverse degrees |
| `R3_ABC_SC` | dropped: ABC_SC has no source |
| `R3_ISI` | `R3_ISI` (ReZ_2 in godlin2026) |
| `R3_sigma` | dropped: no source applies sigma to modified reverse degrees |
| `R3_HBM` | dropped: one of the twelve hybrid descriptors of arockiaraj2023sar (per its abstract); full text not accessible and no open source restates the formula, so the definition could not be checked |
| `R3_HTM` | dropped: one of the twelve hybrid descriptors of arockiaraj2023sar (per its abstract); full text not accessible and no open source restates the formula, so the definition could not be checked |
| `R3_BMA` | dropped: = 2m + 2 R3_ISI identically; formula not verifiable |
| `R3_TMA` | dropped: = 2 R3_M1 - 2 R3_ISI identically; formula not verifiable |
| `R4_A` | dropped: R4_A = R4_M1 / 2 exactly; additionally k = 4 exceeds Delta(G) = 3 for most drugs (v0.x silently dropped the zero-weight edges) |
| `R4_G` | dropped: k = 4 not in the registered range 1..3 (see R_k decision) |
| `R4_HA` | dropped: no primary source found; additionally k = 4 exceeds Delta(G) = 3 for most drugs (v0.x silently dropped the zero-weight edges) |
| `R4_SO` | dropped: no source applies the Sombor function to modified reverse degrees; additionally k = 4 exceeds Delta(G) = 3 for most drugs (v0.x silently dropped the zero-weight edges) |
| `R4_ABC_SC` | dropped: ABC_SC has no source; additionally k = 4 exceeds Delta(G) = 3 for most drugs (v0.x silently dropped the zero-weight edges) |
| `R4_ISI` | dropped: k = 4 not in the registered range 1..3 (see R_k decision) |
| `R4_sigma` | dropped: no source applies sigma to modified reverse degrees; additionally k = 4 exceeds Delta(G) = 3 for most drugs (v0.x silently dropped the zero-weight edges) |
| `R4_HBM` | dropped: one of the twelve hybrid descriptors of arockiaraj2023sar (per its abstract); full text not accessible and no open source restates the formula, so the definition could not be checked; additionally k = 4 exceeds Delta(G) = 3 for most drugs (v0.x silently dropped the zero-weight edges) |
| `R4_HTM` | dropped: one of the twelve hybrid descriptors of arockiaraj2023sar (per its abstract); full text not accessible and no open source restates the formula, so the definition could not be checked; additionally k = 4 exceeds Delta(G) = 3 for most drugs (v0.x silently dropped the zero-weight edges) |
| `R4_BMA` | dropped: = 2m + 2 R4_ISI identically; formula not verifiable; additionally k = 4 exceeds Delta(G) = 3 for most drugs (v0.x silently dropped the zero-weight edges) |
| `R4_TMA` | dropped: = 2 R4_M1 - 2 R4_ISI identically; formula not verifiable; additionally k = 4 exceeds Delta(G) = 3 for most drugs (v0.x silently dropped the zero-weight edges) |
