# Scientific research plan (audit 2026-10-12)

Everything below is **new and exploratory** unless it is registered externally *before* the
corresponding analysis is run. Nothing here alters the frozen plan (`docs/EXPERIMENT_PLAN.md`,
tag `plan-freeze-v1`) or any archived result. Each experiment gets its own config, output
directory and dated amendment.

## 1. Mathematical statement: representation-limited empirical error floor

**Setting.** Finite data {(G_i, y_i)}, i = 1..n, y_i ∈ ℝ. Feature map φ: graphs → ℝ^p (here
the 102 frozen indices of the hydrogen-suppressed, unlabeled heavy-atom graph, optionally
concatenated with the 11 composition counts). Equivalence G_i ~ G_j ⇔ φ(G_i) = φ(G_j). The
classes are C_1..C_K, with class means ȳ_C.

**Proposition 1.** For every function f: ℝ^p → ℝ,

  (1/n) Σ_i (y_i − f(φ(G_i)))² ≥ (1/n) Σ_C Σ_{i∈C} (y_i − ȳ_C)² =: L_φ,

and equality holds for f(φ(G)) = ȳ_{C(G)}.

*Proof.* f∘φ is constant on each class, say c_C. The loss splits over classes:
Σ_C Σ_{i∈C} (y_i − c_C)². For each class, Σ_{i∈C}(y_i − c)² = Σ_{i∈C}(y_i − ȳ_C)² + |C|(ȳ_C − c)²
(add and subtract ȳ_C; the cross term vanishes because Σ_{i∈C}(y_i − ȳ_C) = 0). The minimum
over c is at c = ȳ_C. Summing over classes gives the bound. The class-mean function attains it,
and it is well defined because φ takes one value per class. ∎

**Corollary (empirical R²).** If S_yy = Σ_i (y_i − ȳ)² > 0, then for every f,
R²_emp(f) = 1 − Σ(y_i − f)²/S_yy ≤ 1 − n·L_φ/S_yy, the between-class share of variance (the
η² of a one-way ANOVA on the classes).

**What it does not say.**

1. It is in-sample and empirical. It bounds every model's *training* fit on these molecules, and
   so a fortiori any CV estimate that scores the same molecules, *provided* the CV predictor
   is a function of φ alone. It is not a bound on the population risk for unseen molecules.
   That would need assumptions on how y varies within classes in the population.
2. Measurement noise is inside L_φ. A class of replicate measurements of one compound
   contributes its replicate variance. Here collisions are mostly between *different*
   compounds (F9), so L_φ is a lower bound on the representation's error plus noise. It
   says nothing about splitting them.
3. Floating-point equality: values computed from isomorphic graphs are bitwise equal.
   Distinct graphs with equal indices are genuine collisions. The script rounds to 9
   significant digits; with tolerances that coarse, near-collisions (not exact ones) could
   merge. Sensitivity to the tolerance should be reported.
4. Any injective-on-data φ gives L_φ = 0. Adding features can only lower L_φ (refinement of
   the partition). The bound is therefore informative only for *poor* representations, or for
   representations evaluated on isomer-rich data.

**Test results (this audit).** On the four endpoint sets the bound is far from binding
(R² ceilings: TI 0.94–0.997; TI+composition ≥ 0.995; size alone 0.36–0.64). The observed CV
R² of TI models (0.2–0.4) is limited by learning and generalisation, not by collisions. The
proposition is correct but **not a headline result** for these data. It becomes one on
isomer-controlled sets (item E3).

**Small-graph tests** (to add as unit tests): path P3 (ethanol vs dimethyl ether, all 104
indices equal, checked); hexane vs 1-pentanol vs pentylamine vs 1,4-butanediol (already cited in
`main.tex:57`); cycles C_n vs C_n-heteroatom variants; star K_{1,3}; disconnected inputs
(salts are reduced to the parent, so not applicable after standardisation).

**Lean 4.** One compact, meaningful target is Proposition 1 over `Finset` with real-valued
labels: "the class-mean predictor minimises squared error among functions constant on a
setoid's classes". It is about 60–100 lines with Mathlib (`Finset.sum_sub_sq`-style lemmas).
Not attempted in this audit; no Lean toolchain was used. Do not claim it formalised until
`lake build` passes.

## 2. Experiments (in order of information value)

### E1. Composition-conditional gain of TI, with fold-safe transform — RUN in this audit (exploratory)

| Field | Answer |
|---|---|
| Research question | Does ΔR²(TI+composition − composition) > 0 on measured endpoints? Falsifiable: the paired corrected-t interval excludes zero on the positive side. |
| Dataset | The frozen pooled set (ChEMBL_37 approved parents); endpoints Bradley MP (CC0), OPERA logP, AZ logD7.4 (ChEMBL CC BY-SA 3.0), AqSolDB logS (CC0); identity by full InChIKey; n = 711/505/369/351 |
| Target | Measured; caveats F10, F11, F14 |
| Predictor sets | composition (11), TI (102, frozen registry), TI+composition, the same with asinh(TI), rdkit2d raw/asinh, size (2) |
| Split | (a) 5×5 random grouped by InChIKey connectivity group, seeds 0–4; (b) Bemis–Murcko scaffold groups (acyclic singletons). Same registry in both. |
| Model selection | In-fold over {ridge, elastic net, PLS, kNN, SVR, extra trees} with 3-fold inner CV, as frozen |
| Primary metric | Fold-mean ΔR² and ΔMAE, paired on identical outer folds |
| Uncertainty | Corrected repeated-CV t (Nadeau–Bengio/Bouckaert–Frank), 95%; no multiplicity claim (exploratory) |
| Counterexample | ΔR² ≤ 0 with the interval excluding +0.02 would refute "TI adds information beyond composition" for that endpoint |
| Reproducibility | `scripts/sensitivity_composition_asinh.py`; asserts exact reproduction of archived raw sets; outputs `~/topoqspr_runs/audit_2026-10-12_sensitivity/` |

### E2. Similarity-cluster split — NOT RUN (plan)

Butina clustering on Morgan r=2, 2048-bit Tanimoto, threshold 0.6 (decided now, before seeing
any outcome). Clusters become CV groups; same learners and folds otherwise. Report the
per-fold distribution of nearest-training-neighbour similarity to quantify the shift. Compute is
about 1 h per feature-set ladder on 14 cores (from the scaffold run's 55 min).

### E3. Isomer-controlled evaluation — NOT RUN (plan; highest *scientific* value after E1)

Within-formula groups (same molecular formula, ≥ 2 members with measured values) from the
full Bradley and OPERA logP sets, not restricted to approved drugs. Question: do TI explain
*within-formula* variance? This is the classical chemical-graph-theory test (octane isomers),
which drug QSPR drops. It holds size *and* composition fixed by construction. Metric: pooled
within-group R² under leave-one-group-out, plus the Proposition 1 ceiling per group.
Refuted if TI ≤ 0 within-group R² while ECFP/labelled descriptors are > 0.

### E4. Conventional ECFP baseline — NOT RUN (plan)

Morgan r=2, 2048 bits (counts), with ridge and extra trees, inside the same nested CV. Needed
because reviewers in 2026 will ask; it is the standard 2D baseline, cheaper than a GNN. A GNN
(Chemprop D-MPNN) is optional and exploratory only.

### E5. Applicability-domain-stratified reporting — partly available

The archived predictions carry `in_domain`. Report full-population and in-domain metrics side
by side for every set. The paper does this for the ladder already; extend it to the paired ΔR².
In-domain numbers never replace the full-population headline.

### E6. Data card and suspect values

Per endpoint: source, licence (F14), units, censoring counts, replicate-SD drops,
overlap with other endpoints, MNAR selection (only 959 of 1821 drugs have any endpoint; approved
drugs with measured values are not a random sample of chemistry). List suspect values (F10)
without removing them.

## 3. What not to do

* Do not replace the frozen headline with asinh or in-domain numbers. Report them beside it.
* Do not add more indices. F9 and the E1 result decide whether *labelled* information is the gap.
* Do not tune the transform, the AD threshold or the cluster cutoff after seeing ΔR².
* Do not describe Experiment D as showing equivalence (F6).
* Do not call the plan preregistered (F4).
