# 03b. Index definitions checked against their primary sources

Date: 2026-10-13. Scope: the 104 indices in `src/topoqspr/indices.py` (`REGISTRY`).

## Why this is separate from 03a

`03a_INDEX_VERIFICATION_INDEPENDENT.md` showed that the code agrees with a reference implementation written from `docs/INDEX_VERIFICATION.md`. That is **implementation agreement**: the code computes what the documentation says. It does not show that the documentation says what the literature says.

This file checks that missing link: **documented definition ↔ definition in the paper that introduced the index**. I worked from the registry metadata (names, families, LaTeX, `ref`) and `docs/index_refs.bib`. I did not use the repo's own "verified" claims as evidence. Where I could not read a source myself, I say so, even when `INDEX_VERIFICATION.md` reports that someone else read it.

Classification used below:

- **MATCHES PRIMARY**: I read the defining paper (full text, an author's accepted version, or page images), and its formula is the documented one.
- **MATCHES A LATER STANDARD VARIANT**: I read the primary, it differs, and the package follows a later convention that is now standard (cited).
- **DIFFERS**: the documented definition disagrees with the primary in a way that changes values. A witness graph is given.
- **SOURCE INACCESSIBLE**: I could not read the primary. The definition was checked only against restatements, which are named.

## How sources were accessed (and what failed)

Downloads were anonymous HTTP requests: arXiv, Europe PMC REST full text, nature.com PDFs, the MATCH archive (match.pmf.kg.ac.rs), NOPR (niscpr.res.in), the scidar.kg.ac.rs repository, and Crossref metadata. No e-mail address or other identifier was sent. No PDFs are stored in the repo; URLs are listed in the source register at the end.

Access that failed:

- Paywalled (Elsevier, ACS, Springer, Wiley, AIP): Gutman & Trinajstić 1972, Gutman et al. 1975, Randić 1975, Wiener 1947, Balaban 1982, Dobrynin & Kochetova 1994, Gutman 1994 (Schultz), Zhou & Trinajstić 2009, Vukičević & Furtula 2009, Furtula et al. 2010, Furtula & Gutman 2015, Došlić et al. 2018, Arockiaraj et al. 2019 (Crossref abstract only), Arockiaraj, Greeni & Kalaam 2023 (IJQC), Cao et al. 2014, Miličević et al. 2004.
- Returned HTTP 403: Khalifeh et al. 2008 (ScienceDirect) and Hosoya 1971 (OUP).
- Hrčak served a bot challenge instead of the PDF: Vukičević & Gašperov 2010 and Gutman & Ashrafi 2008.
- Host unreachable: Shirdel et al. 2013 (IJMC) and Graovac et al. 2011 (JMNS).
- No online copy found: Fajtlowicz 1987, Ranjini et al. 2013, Shegehalli & Kanabur 2015, Gutman 1994 (Graph Theory Notes NY, Szeged), Khadikar 2000, Ghorbani & Hosseinzadeh 2010.

## Conventions shared by all indices

| Point | Package (from `INDEX_VERIFICATION.md`) | What the sources say |
|---|---|---|
| Graph | Hydrogen-suppressed, simple; bond order ignored | Kulli 2018 (§1: "undirected without loops and multiplied edges"), Gutman 2021 Sombor (§1: "simple graphs ... without directed, weighted or multiple edges, and without self loops") and Arockiaraj 2025 ("hydrogen-suppressed molecular graph") all say the same. Multigraph versions of the Zagreb indices exist in the literature but are not used by any registered source. |
| Degree | d_v = number of heavy-atom neighbours | Same in every source read. |
| Edge-additive form | Σ_{uv∈E} φ(x_u, x_v) with x = d, S or M_kR(d) | Exactly Eqs. (2)–(4) of Arockiaraj et al. 2025 and Eqs. (2)–(3) of Godlin & Radha 2026 (both read). |
| Vertex sums vs edge sums | M1 and F are written as edge sums | Identical for every graph: Σ_v d_v^{α+1} = Σ_{uv}(d_u^α + d_v^α). This identity is stated in Gutman 2021 (MATCH 86, p. 12) and in Mondal et al. (arXiv:1906.06660, p. 2). |

## Per-index results

### Degree family (32)

| Index | Documented φ(d_u, d_v) | Defining source (DOI) | Where | What I read | Class |
|---|---|---|---|---|---|
| M1 | d_u + d_v | Gutman & Trinajstić 1972, 10.1016/0009-2614(72)85099-1 | — | Restatements: Mondal arXiv:1906.06660 p. 2 (Σ d_v² = Σ (d_u + d_v)); Gutman 2021 Table 1 | INACCESSIBLE (restatements agree) |
| M2 | d_u d_v | Gutman et al. 1975, 10.1063/1.430994 | — | Mondal p. 2; Gutman 2021 Table 1 | INACCESSIBLE (restatements agree) |
| F | d_u² + d_v² | Furtula & Gutman 2015, 10.1007/s10910-015-0480-z | — | Mondal p. 2 (Σ d³ = edge form); Gutman 2021 Table 1; Arockiaraj 2023 SAR QSAR p. 570, text | INACCESSIBLE (restatements agree) |
| HM | (d_u + d_v)² | Shirdel et al. 2013, 10.22052/ijmc.2013.5294 | — | Gutman 2021 Table 1 ("first hyper-Zagreb"); Arockiaraj 2023 p. 570 (F and M2 combined 1:2) | INACCESSIBLE (restatements agree) |
| R | (d_u d_v)^{-1/2} | Randić 1975, 10.1021/ja00856a001 | — | Estrada et al. 1998 Eq. (1); Molecules 21 (2016) 821 Eq. (1); Ali et al. ABS Eq. (1). All use exponent −1/2. | INACCESSIBLE (restatements agree) |
| RR | √(d_u d_v) | Gutman, Furtula & Elphick 2014, MATCH 72:617–632 | p. 618, Eq. (5) | **Full text** | MATCHES PRIMARY |
| SC | (d_u + d_v)^{-1/2} | Zhou & Trinajstić 2009, 10.1007/s10910-008-9515-z | — | Mondal p. 2 ("SCI"); Gutman 2021 Table 1 | INACCESSIBLE (restatements agree) |
| H | 2/(d_u + d_v) | Fajtlowicz 1987 (no DOI) | — | Ali et al. ABS Eq. (2); Arockiaraj 2025 list; Gutman 2021 Table 1 | INACCESSIBLE (restatements agree) |
| GA | 2√(d_u d_v)/(d_u + d_v) | Vukičević & Furtula 2009, 10.1007/s10910-009-9520-x | — | Molecules 2016 Eq. (3); Gutman 2021 Table 1; Arockiaraj 2025 list | INACCESSIBLE (restatements agree) |
| ABC | √((d_u + d_v − 2)/(d_u d_v)) | Estrada, Torres, Rodríguez & Gutman 1998, Indian J. Chem. 37A:849–855 | p. 850 Eq. (2); p. 851 Eq. (3) | **Full text** (NOPR scan) | **MATCHES A LATER STANDARD VARIANT** (see below) |
| AZI | (d_u d_v/(d_u + d_v − 2))³ | Furtula, Graovac & Vukičević 2010, 10.1007/s10910-010-9677-3 | — | Gutman 2021 Table 1 | INACCESSIBLE (restatement agrees) |
| SO | √(d_u² + d_v²) | Gutman 2021, MATCH 86:11–16 | p. 14, Eq. (3) | **Full text** | MATCHES PRIMARY (the basic SO; the same paper's SO_red and SO_avg are not used) |
| ISI | d_u d_v/(d_u + d_v) | Vukičević & Gašperov 2010, Croat. Chem. Acta 83:243–260 | — | Gutman 2021 Table 1 ("inverse sum indeg") | INACCESSIBLE (restatement agrees) |
| ReZ1 | (d_u + d_v)/(d_u d_v) | Ranjini, Lokesha & Usha 2013 (no DOI) | — | Godlin & Radha 2026 p. 2 list; Arockiaraj 2025 list | INACCESSIBLE (restatements agree) |
| mM2 | 1/(d_u d_v) | Cited: Miličević et al. 2004, 10.1023/B:MODI.0000047504.14261.2a. The index is from Nikolić et al. 2003, Croat. Chem. Acta 76:113 (the bib note says so). | — | Not read; no restatement read either | INACCESSIBLE (no restatement read; citation points to a later paper) |
| SDD | (d_u² + d_v²)/(d_u d_v) | Vukičević & Gašperov 2010 | — | Mondal p. 2 (Σ [d_u/d_v + d_v/d_u]); Godlin 2026 list ("SSD"); Arockiaraj 2025 | INACCESSIBLE (restatements agree) |
| sigma | (d_u − d_v)² | Gutman, Togan, Yurttas, Cevik & Cangul 2018, MATCH 79:491–508 | Eq. (5), §1 | **Full text** | MATCHES PRIMARY |
| ABS | √((d_u + d_v − 2)/(d_u + d_v)) | Ali, Furtula, Redžepović & Gutman 2022, 10.1007/s10910-022-01403-1 | Defining display after Eq. (1) | **Authors' final version** (scidar.kg.ac.rs) | MATCHES PRIMARY |
| ReZ3 | d_u d_v (d_u + d_v) | Ranjini et al. 2013 | — | Mondal p. 3 ("ReZG3"); Godlin 2026 list; Gutman 2021 Table 1 lists the same φ as the "second Gourava index" | INACCESSIBLE (restatements agree) |
| AG | (d_u + d_v)/(2√(d_u d_v)) | Shegehalli & Kanabur 2015 | — | Gutman 2021 Table 1; Godlin 2026 list | INACCESSIBLE (restatements agree) |
| BM, TM | d_u + d_v + d_u d_v; d_u² + d_v² + d_u d_v | Arockiaraj, Paul, Clement, Tigga, Jacob & Balasubramanian 2023, 10.1080/1062936X.2023.2239149 | pp. 570–571 | **Page images** (`related_papers/…photos.pdf`, pp. 570–572) | MATCHES PRIMARY (naming note below) |
| GH, GBM, GTM, HG, HBM, HTM, BMG, BMH, TMG, TMH | as in `INDEX_VERIFICATION.md` | same paper | pp. 571–572 | **Page images** | MATCHES PRIMARY. All 10 formulas agree symbol by symbol. |

**ABC.** The 1998 paper defines ABC = Σ_e [δ(v_i)δ(v_j)/(2δ(e))]^{-1/2} (Eq. (2), p. 850), with δ(e) = δ(v_i) + δ(v_j) − 2, and rewrites it as Σ [2(δ_i + δ_j − 2)/(δ_i δ_j)]^{1/2} (Eq. (3)). That is √2 times the package's ABC. The package follows the later convention that drops the √2. Ali et al. 2022, p. 2, record the change: "in the original definition of the ABC index a factor √2 was present ... but it was dropped later". Molecules 2016 Eq. (2), Gutman 2021 Table 1 and Mordred all use the √2-free form. Witnesses:

| Graph | Package ABC | 1998 ABC |
|---|---|---|
| P4 | 3/√2 = 2.1213 | 3.0000 |
| K1,3 | √6 = 2.4495 | 3.4641 |
| C4 | 2√2 = 2.8284 | 4.0000 |

The factor is a constant, so correlations and R² are unchanged. Absolute values can only be compared with pre-2008 tables after scaling by √2. The same applies to `R{k}_ABC`. `S_ABC` (ABC₄) was defined without the √2 (Molecules 2016 Eq. (4)), so it is unaffected.

**Naming note (BM, ReZ3, HCM2).** Gutman 2021 Table 1, which I read, lists x + y + xy as the "first Gourava index" and (x + y)xy as the "second Gourava index" (Kulli 2017). So the bi-Zagreb index BM of Arockiaraj et al. 2023 is numerically the first Gourava index, which was published earlier. ReZ3 has the same formula as the second Gourava index. The φ of HCM2, (xy)², is the "second hyper-Zagreb" function of that table. The values are unaffected. Only priority and naming change: the `ref` of BM is a valid source, but it is not the earliest.

### Reverse-degree family (49 = 16 φ × k ∈ {1, 2, 3} + HCM2)

Documented map: M_kR(d) = Δ − d + k if k ≤ d, else (Δ − d + k) mod Δ, with k restricted to 1..Δ. HCM2 uses c = Δ − d + 1.

| Item | Source | What I read | Finding |
|---|---|---|---|
| c_v = Δ − d_v + 1; CM1, CM2 (eq. 1); HCM1, HCM2 = Σ (c_u c_v)² (eq. 2) | Kulli 2018, Ann. Pure Appl. Math. 16:47–51, 10.22457/apam.v16n1a6 | **Full text**, p. 47 | Matches `HCM2`, `R1_M1` (= CM1), `R1_M2` (= CM2) and `R1_HM` (= HCM1). Kulli attributes CM1 and CM2 to Ediz (his ref. [3]), which I did not access. |
| Generalised map M_kR, original | Arockiaraj, Greeni & Kalaam 2023, IJQC 123:e27136, 10.1002/qua.27136 | **Not read** (paywalled) | Inaccessible to me |
| M_kR restated by the same group | Kalaam, Greeni & Arockiaraj 2024, Front. Chem. 12:1470231, §2 | **Full text** (Europe PMC) | Printed as "Δ(G) − d_v + k mod Δ(G) : k > d_v", without parentheses. The Δ = 4 tables in §3 (k = 2: d = 1→1; k = 3: d = 1→2, d = 2→1) only fit the reading with the mod applied to the whole sum. That is the package's reading. |
| Same, with explicit parentheses | Arockiaraj et al. 2025, Sci. Rep. 15:3639, Eq. (1); Godlin & Radha 2026, Sci. Rep. 16:8889, Eq. (1) | **Full text** (PDF) | Both print (Δ(G) + k − d) mod Δ(G), matching the package. |
| φ for M1, M2, F, HZ (= HM), ReZ3, BM, TM, BMH, TMH | Kalaam 2024 Eqs. (1)–(9) | **Full text** | Matches. Eq. (7) is labelled "BM" but named tri-Zagreb, and its formula is TM's. |
| φ for ABC, R, ReZ1, SSD, H, G (= RR), GBM, TMH, GA with M_kR | Arockiaraj 2025, list after Eq. (4) | **Full text** | Matches the package's φ for ABC, R, ReZ1, RR and GA. |
| φ for ABC, R, ReZ1, ReZ2 (= ISI), ReZ3, SSD, AG, HG, BMG, GBM | Godlin & Radha 2026, list after Eq. (3) | **Full text** | Matches. pdftotext scrambled the ReZ1/ReZ2 fractions; ISI = xy/(x + y) was confirmed from Gutman 2021 Table 1 instead. |

Classification:

- **MATCHES PRIMARY (19).**
  - `HCM2` (Kulli 2018).
  - R_k for HM, BM, TM, BMH, TMH × 3 (15). Kalaam 2024 is the earliest source I found for these φ on M_kR. Going by the IJQC list of ten φ quoted in `INDEX_VERIFICATION.md`, they are not in the 2023 paper.
  - `R{1,2,3}_R` (3). The earliest source I read for R on M_kR is Arockiaraj 2025; the registry cites godlin2026, which is later.
- **SOURCE INACCESSIBLE (30).** R_k for M1, M2, F, ReZ3, GA, AG, ABC, RR, ReZ1, ISI × 3. Their original is the IJQC 2023 paper, which I did not read. The definitions agree with three open restatements by the same group, which I did read, and with Kulli 2018 for k = 1 (M1, M2).

Two convention points with numeric consequences:

1. **Precedence of mod.** Read literally without parentheses (Δ − d + (k mod Δ)), K1,3 with k = 2 gives leaf weight 4 instead of 1. Then `R2_M1`(K1,3) = 18 instead of 9. The package's reading is the one forced by Kalaam's own tables and printed explicitly by Arockiaraj 2025 and Godlin 2026.
2. **Domain k ≤ Δ (a documented deviation from the literal formula).** The sources state the parameter as "k ≥ 1" with no upper bound (Kalaam §2; Arockiaraj 2025 Eq. (1)). The package returns NaN when k > Δ. Witness P4 (Δ = 2), k = 3: the formula gives M_3R = (0, 1, 1, 0), so the source formula yields `R3_M1`(P4) = 4, while the package gives NaN. φ with a reciprocal (R, ABC, …) would divide by zero. In practice the sources keep k ≤ Δ ("limited to l from 1 to 3", Godlin 2026). This is a domain restriction, not a different formula. It affects only molecules with Δ < 3 at k = 3.

### Neighbourhood-degree-sum family (10)

S_v = Σ_{u∈N(v)} d_u (open neighbourhood).

| Index | Source | What I read | Class |
|---|---|---|---|
| S_RR, S_SC, S_ReZ3, S_R, S_SDD (ND1–ND5) | Mondal, De & Pal, arXiv:1906.06660, pp. 2–3; journal version Complex Intell. Syst. 7 (2021) 977, 10.1007/s40747-020-00262-0 | **Preprint full text**; journal version not read | MATCHES PRIMARY (preprint). δ_G(u) is defined as the sum over N_G(u), the open neighbourhood, as in the package. ND6 is printed as "[d(u)δ(u) + d(v) + δ(v)]" and is not registered, which is consistent with the repo. |
| S_H, S_GBM, S_TMH | Arockiaraj et al. 2025, Eq. (3) with H, GBM, TMH from their list | **Full text** | MATCHES PRIMARY |
| S_ABC (ABC₄) | Ghorbani & Hosseinzadeh 2010 (no DOI) | Molecules 21 (2016) 821, Eq. (4), with S_u = Σ_{v∈N(u)} d_v | INACCESSIBLE (restatement agrees, no √2) |
| S_GA (GA₅) | Graovac et al. 2011, 10.22061/jmns.2011.461 | Molecules 2016 Eq. (5) | INACCESSIBLE (restatement agrees) |

### Distance family (11), Hosoya Z, degree entropy

| Index | Documented | Source | What I read | Class / convention notes |
|---|---|---|---|---|
| W | Σ_{unordered pairs} D(u,v) | Wiener 1947, 10.1021/ja01193a005 | Not read. Restatements: Molecules 2016 §1 ("path number" = sum of distances between all pairs); arXiv:1809.02358 §1 | INACCESSIBLE. Wiener defined the path number for paraffins, i.e. trees. Extending it to cyclic graphs as the distance-matrix half-sum is the standard later reading, usually credited to Hosoya 1971; I did not read that either. |
| J | m/(μ+1) Σ_{uv∈E} (D_u D_v)^{-1/2}, μ = m − n + 1 | Balaban 1982, 10.1016/0009-2614(82)80009-2 | Not read. Knor et al. arXiv:1701.02716 p. 2: J = m/(m − n + 2) Σ 1/√(w(u)w(v)), w = distance sum | INACCESSIBLE (restatement agrees; m − n + 2 = μ + 1 for connected graphs) |
| Sz | Σ n_u n_v, strictly closer | Gutman 1994, Graph Theory Notes NY 27:9 | Not read. Črepnjak & Tratnik arXiv:1705.02841 §2 (N1, N2 strict); Ashrafi & Mirzargar, MATCH 60 (2008) 897, §1 ("equidistant vertices are not taken into account") | INACCESSIBLE (restatements agree) |
| Sz_e | Σ m_u m_v, with d(f, x) = min over the endpoints of f, strict | Gutman & Ashrafi 2008, Croat. Chem. Acta 81:263 | Not read (Hrčak bot wall). Črepnjak & Tratnik §2 (d(x, e) = min, strict); MATCH 60 (2008) 897 ("Edges equidistant from both ends of an edge are not counted") | INACCESSIBLE (restatements agree) |
| PI (edge) | Σ (m_u + m_v) | Khadikar 2000, Natl. Acad. Sci. Lett. 23:113 | Not read. MATCH 60 (2008) 897 §1, quoting Khadikar's definition; Črepnjak & Tratnik §2 ("PI = PI_e") | INACCESSIBLE (restatements agree) |
| PI_v | Σ (n_u + n_v) | Khalifeh, Yousefi-Azari & Ashrafi 2008, 10.1016/j.dam.2007.08.041 | Not read (403); no restatement read | INACCESSIBLE |
| Mo | Σ \|n_u − n_v\| | Došlić et al. 2018, 10.1007/s10910-018-0928-z | Not read. Tratnik arXiv:1904.04131 abstract and §1 | INACCESSIBLE (restatement agrees) |
| Mo_e | Σ \|m_u − m_v\| | Arockiaraj, Clement & Tratnik 2019, 10.1002/qua.26043 | Crossref abstract only. It introduces "edge-Mostar" but the abstract gives no formula. | INACCESSIBLE (no formula read) |
| W_e | W(L(G)) | Iranmanesh, Gutman, Khormali & Mahmiani 2009, MATCH 61:663–672 | **Full text**: Def. 2a/2b, Eqs. (3)–(4) (p. 666): W_e0 = W(L(G)); Conclusions (p. 670) recommend W_e0 or W_e4 | MATCHES PRIMARY (W_e0 variant) |
| DD | Σ_{pairs} (d_u + d_v) D(u,v) | Dobrynin & Kochetova 1994, 10.1021/ci00021a008 | Not read. arXiv:1809.02358 §1 (same formula, unordered pairs) | INACCESSIBLE (restatement agrees) |
| Gut | Σ_{pairs} d_u d_v D(u,v) | Gutman 1994, 10.1021/ci00021a009 | Not read. arXiv:1809.02358 §1 | INACCESSIBLE (restatement agrees) |
| Z | Σ_k p(G, k), p(G, 0) = 1 | Hosoya 1971, 10.1246/bcsj.44.2332 | Not read (403); no restatement read in this pass | INACCESSIBLE |
| I_deg | −Σ (d/2m) log₂(d/2m) | Cao, Dehmer & Shi 2014, 10.1016/j.ins.2014.03.133 | Not read. Hua et al. arXiv:1709.09594 Def. 1 attributes I_d^t to Cao et al. with base 2. Cambie et al. arXiv:2206.00922 Def. 2 uses the natural log for the same entropy. | INACCESSIBLE; **the log base is contested** in the literature |

## Where the literature disagrees, and which convention the package follows

| Topic | Variants in the literature (sources read) | Package | Witness |
|---|---|---|---|
| ABC normalisation | 1998 original carries √2 (Estrada 1998 Eq. (3)); later dropped (Ali et al. 2022; Gutman 2021 Table 1) | √2 dropped | P4: 2.1213 vs 3.0 |
| Sombor | SO, SO_red, SO_avg (Gutman 2021 Eqs. (3)–(5)) | SO | P4: SO = 2√5 + 2√2 = 7.3006; SO_red = 2 + √2 = 3.4142 |
| Randić exponent | −1/2 (Randić, as restated); general α (RR is α = +1/2, Gutman 2014 p. 619) | R: −1/2; RR: +1/2 | — |
| Edge Wiener | W_e0 = W(L(G)); W_e1 (min endpoint distance) = W_e0 − m(m−1)/2; W_e2, W_e4 (max distance) (Iranmanesh 2009) | W_e0. ChemGraphX uses W_e1, as found in the repo's ChemGraphX cross-check. | P4: 4 vs W_e1 = 1 vs W_e4 = 7. K1,3: 3 / 0 / 6. C4: 8 / 2 / 12. |
| PI naming | "PI" = edge PI (Khadikar; MATCH 60; Črepnjak & Tratnik); PI_v = vertex version (Khalifeh 2008) | PI = edge, PI_v = vertex | P4: PI = 6, PI_v = 12. C4: 8 vs 16. |
| Edge "closer to u" | d(f, u) = min over the endpoints of f (Črepnjak & Tratnik §2) | min, strict | C5: 20 with the min rule. A stricter "both endpoints closer" rule gives 10. That alternative is my own construction for illustration, not taken from a cited paper. On bipartite graphs the two rules agree. |
| W on cyclic graphs vs Sz | Sz was introduced as the cycle-aware extension of Wiener's tree formula. The package keeps both. | W = distance half-sum | C4: W = 8, Sz = 16. Equal on trees. |
| R_k mod precedence | Literal Kalaam text without parentheses vs (Δ + k − d) mod Δ | whole-sum mod | K1,3, k = 2: R2_M1 = 9 vs 18 |
| R_k for k > Δ | Formula stated for all k ≥ 1 | NaN for k > Δ | P4, k = 3: 4 vs NaN |
| Neighbourhood sum | Open N(v) (Mondal; Arockiaraj 2025; Molecules 2016) | open | — |
| Entropy log base | base 2 (Hua et al.) vs natural log (Cambie et al.) | base 2 | K1,3: 1.7925 bits vs 1.2425 nats (constant factor ln 2) |
| Zagreb on multigraphs | Not used by any registered source; all of them use simple graphs | simple hydrogen-suppressed graph | — |

The witness values come from `index_verification/primary_sources/witness_values.py`. It is pure Python and does not import the package.

## Summary

### Counts (104 indices)

| Class | Count | Indices |
|---|---|---|
| MATCHES PRIMARY | 44 | RR, SO, sigma, ABS, and the 12 hybrids (BM, TM, GH, GBM, GTM, HG, HBM, HTM, BMG, BMH, TMG, TMH); HCM2; R_k × {HM, BM, TM, BMH, TMH, R} (18); S_RR, S_SC, S_ReZ3, S_R, S_SDD, S_H, S_GBM, S_TMH; W_e |
| MATCHES A LATER STANDARD VARIANT | 1 | ABC (1998 original × √2) |
| DIFFERS | 0 | No formula disagrees with its primary beyond the ABC constant. One documented domain restriction exists (R_k is NaN for k > Δ, where the literal source formula would return a number). |
| SOURCE INACCESSIBLE | 59 | M1, M2, F, HM, R, SC, H, GA, AZI, ISI, ReZ1, mM2, SDD, ReZ3, AG; R_k × {M1, M2, F, ReZ3, GA, AG, ABC, RR, ReZ1, ISI} (30); S_ABC, S_GA; W, J, Sz, Sz_e, PI, PI_v, Mo, Mo_e, DD, Gut; Z; I_deg |

Of the 59 inaccessible ones, 54 have at least one restatement I read that agrees with the documented formula. Five have none: mM2, PI_v, Mo_e (abstract only), Z, and I_deg (the restatements disagree on the log base).

### Three kinds of evidence, kept separate

| Kind | Status | Where |
|---|---|---|
| (a) Implementation agreement: the code computes the documented formula | Established. 104/104 indices; 111,920 comparisons; 0 code bugs. | `03a_INDEX_VERIFICATION_INDEPENDENT.md` (not redone here) |
| (b) Definition agreement: the documented formula is the primary source's formula | 44 checked against the primary; 1 follows a later standard (ABC); 59 checked only against restatements | this file |
| (c) Values printed in primary papers reproduced | Summarised from `LITERATURE_VALUE_CHECK.md` and the worked-example tests in `INDEX_VERIFICATION.md`; not redone | see below |

What (c) covers:

- `LITERATURE_VALUE_CHECK.md` reports 199 test cases from seven or eight papers by other groups (Youssef 2024, Albadrani 2022, Rao 2024, Abubakar 2024, Almulla 2026, Ravi & Desikan 2023, Ravi & Chidambaram 2025, Mondal/De/Pal 2021 Cayley trees). They are almost all **reverse-degree** and **neighbourhood-degree-sum** indices, plus degree indices computed from bond partitions. Every discrepancy found was traced to the papers, not to the package.
- `INDEX_VERIFICATION.md` adds:
  - Godlin 2026 hydrochlorothiazide (40 values);
  - Arockiaraj 2025 gemcitabine ABC (5 values);
  - IJQC 2023 dexamethasone/molnupiravir (21 values);
  - SAR QSAR 2023 Table 2 (50 of 51 values).
- The distance indices, Hosoya Z and the entropy have **no** value check against primary papers. Their evidence is closed forms (for example Z(P_n) = F_{n+1}, W(P_n) = n(n² − 1)/6) and agreement with networkx, Mordred and ChemGraphX. Those are implementation cross-checks, not literature values.

**What can be claimed.** A definition is supported as correct where (b) and (c) are positive together and consistent with closed forms. That holds for:

- the hybrids;
- R_k for HM, BM, TM, BMH, TMH and R;
- HCM2;
- the eight neighbourhood indices with a primary read;
- RR, SO, sigma and ABS.

For ABC the support is "matches the post-2008 convention". For the classical distance indices, Z, I_deg and most classical degree indices, the definitions rest on consistent secondary restatements and well-known closed forms, not on a primary read.

## Source register (read in this pass)

| Source | URL used | Access |
|---|---|---|
| Kulli 2018 APAM 16:47 | researchmathsci.org/APAMart/apam-v16n1-6.pdf | full text |
| Kalaam et al. 2024 Front. Chem. | Europe PMC REST, PMC11462630 fullTextXML | full text |
| Arockiaraj et al. 2025 Sci. Rep. 15:3639 | nature.com/articles/s41598-025-88044-x.pdf | full text |
| Godlin & Radha 2026 Sci. Rep. 16:8889 | nature.com/articles/s41598-026-37586-9.pdf | full text |
| Mondal, De & Pal preprint | arxiv.org/pdf/1906.06660 | full text (preprint) |
| Molecules 21 (2016) 821 | Europe PMC, PMC6273492 | full text (secondary) |
| Gutman 2021 (Sombor) | MATCH 86 (2021) 11–16, match.pmf.kg.ac.rs | full text |
| Gutman, Furtula & Elphick 2014 (RR) | MATCH 72 (2014) 617–632 | full text |
| Gutman et al. 2018 (sigma) | MATCH 79 (2018) 491–508 | full text |
| Iranmanesh et al. 2009 (edge Wiener) | MATCH 61 (2009) 663–672 | full text |
| Ashrafi & Mirzargar 2008 | MATCH 60 (2008) 897–904 | full text (secondary) |
| Estrada et al. 1998 (ABC) | nopr.niscpr.res.in/bitstream/123456789/40308 | full-text scan |
| Ali et al. 2022 (ABS) | scidar.kg.ac.rs/bitstream/123456789/15555/1/ABS-final.pdf | authors' final version |
| Arockiaraj et al. 2023 SAR QSAR | `related_papers/Arockiaraj_2023_SARQSAR_hybrid_GHZ_descriptors_photos.pdf`, pp. 570–572 | page images |
| Črepnjak & Tratnik | arxiv.org/pdf/1705.02841 | full text (secondary) |
| Tratnik (Mostar) | arxiv.org/pdf/1904.04131 | full text (secondary) |
| Knor et al. (Balaban) | arxiv.org/pdf/1701.02716 | full text (secondary) |
| DD/Gutman methods | arxiv.org/pdf/1809.02358 | full text (secondary) |
| Hua et al. (entropy) | arxiv.org/pdf/1709.09594 | full text (secondary) |
| Cambie et al. (entropy) | arxiv.org/pdf/2206.00922 | full text (secondary) |
