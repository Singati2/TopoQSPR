# TopoQSPR v1.0 evaluation protocol

This document specifies what `topoqspr.features`, `topoqspr.models`,
`topoqspr.evaluate` and `topoqspr.stats` do. The code is the reference; this
file states the defaults, formulas and output columns so that results can be
read without the source. It implements sections 5–7 of `docs/EXPERIMENT_PLAN.md`
(draft) and sections 5–7, 11, 18 and 19 of `docs/REVIEW_BRIEF.md`.

## 1. Feature sets (the baseline ladder)

All sets are computed from an RDKit `Mol` on the **hydrogen-suppressed graph**
G = (V, E): one vertex per atom with atomic number > 1, one edge per bond
between two such atoms. Explicit H atoms in the input change nothing except
`n_H`.

| Set | Columns | Definition |
|---|---|---|
| `null` | none | M0: ŷ = mean(y_train) (`DummyRegressor`) |
| `size` | `n_heavy`, `n_bonds` | \|V\|, \|E\| |
| `size_rings` | `size` + `n_rings` | n_rings = μ(G) = \|E\| − \|V\| + c(G) (cyclomatic number = SSSR size `len(GetSSSR)`; c = number of connected components. RDKit `CalcNumRings` counts the symmetrised SSSR and can be larger: cubane 6 vs μ = 5) |
| `composition` | `n_C n_N n_O n_S n_P n_F n_Cl n_Br n_I n_other n_H` | element counts of heavy atoms; `n_H` = all H after `AddHs` |
| `rdkit2d` | RDKit `Descriptors.descList` (217 in RDKit 2026.03) minus the target's exclusions | computed on `RemoveHs(mol)` |
| `TI` | supplied by the index module | passed in as a DataFrame and reindexed to the molecules |

**Note on `size_rings`.** Because μ(G) = |E| − |V| + c(G) exactly, ring count
adds information beyond `size` only through c(G) (multi-fragment inputs) for a
linear model. It is kept as the sensitivity model of the plan, but a "gain" of
`size_rings` over `size` in a linear family cannot come from rings.

### rdkit2d exclusions (target leakage)

A descriptor is excluded when it **is** the target, is computed from **the same
per-atom contributions** as the target, or takes the target as an **input**.
Alternative estimators of the same physical quantity are kept, because a
conventional baseline would use them.

| Target (aliases) | Excluded | Why |
|---|---|---|
| `MolWt` (`PC_MolecularWeight`, `PC_ExactMass`, `PC_MonoisotopicMass`) | `MolWt`, `ExactMolWt`, `HeavyAtomMolWt`, `BCUT2D_MWHI`, `BCUT2D_MWLOW`, `qed` | same quantity; BCUT weighted by atomic mass; QED uses MW |
| `LogP` (`MolLogP`; Crippen) | `MolLogP`, `SlogP_VSA1–12`, `BCUT2D_LOGPHI/LOW`, `MolMR`, `SMR_VSA1–10`, `BCUT2D_MRHI/LOW`, `qed` | all derive from Wildman–Crippen atom contributions (same atom typing); QED uses Crippen logP |
| `MolMR` (Crippen) | the Crippen logP and MR families above | same contributions table |
| `TPSA` (`PC_TPSA`) | `TPSA`, `qed` | QED uses PSA. **`LabuteASA` is kept**: it is total approximate surface area from Labute per-atom terms, not Ertl's polar fragment contributions; no other descList entry uses TPSA contributions |
| `HBD` (`PC_HBondDonorCount`) | `NumHDonors`, `NHOHCount`, `qed` | Lipinski donor counts; QED uses HBD |
| `HBA` (`PC_HBondAcceptorCount`) | `NumHAcceptors`, `NOCount`, `qed` | acceptor counts; QED uses HBA |
| `RotBonds` (`PC_RotatableBondCount`) | `NumRotatableBonds`, `qed` | QED uses rotatable bonds |
| `PC_Complexity` | `BertzCT` | PubChem complexity is the Bertz/Hendrickson/Ihlenfeldt formula |
| `PC_XLogP`, `XLogP3` | none | Crippen `MolLogP` is a different, separately fitted estimator |
| `PC_HBD`, `PC_HBA`, `PC_RotatableBonds` | as `HBD`, `HBA`, `RotBonds` | v0.x column names (`config.ML_TARGETS`) for PubChem HBondDonorCount, HBondAcceptorCount, RotatableBondCount |
| `PC_HeavyAtomCount` | `HeavyAtomCount` | same quantity. **The target equals `size.n_heavy` and the sum of the composition heavy-atom counts**, so it is a sanity check, never a QSPR result |
| `PC_AtomStereoCount` | `NumAtomStereoCenters`, `NumUnspecifiedAtomStereoCenters` | PubChem counts defined + undefined stereocentres |
| `PC_FormalCharge` | `NumValenceElectrons`, `Max/MinPartialCharge`, `Max/MinAbsPartialCharge`, `BCUT2D_CHGHI/LO` | NumValenceElectrons subtracts each atom's formal charge; Gasteiger charges start from formal charges (over all atoms incl. H they sum to the net charge) |
| `PC_CovalentUnitCount` | none | no descList entry counts fragments. **`size_rings` encodes it exactly**: c = n_rings − n_bonds + n_heavy |
| `BIO_IC50`, `BIO_Ki`, `BIO_EC50`, `BIO_Kd` | none | measured activities; no descriptor is computed from them |

The BCUT and QED dependencies were checked in RDKit 2026.03
(`rdMolDescriptors.BCUT2D(mol, crippen_contribs)` reproduces
`BCUT2D_LOGPHI/LOW` and `BCUT2D_MRHI/LOW`; `QED.properties` calls
`Crippen.MolLogP`, `MolSurf.TPSA`, MW, HBD, HBA and rotatable bonds).
An unknown target raises `KeyError`: every new target needs an explicit
exclusion decision (`target=None` means "exclude nothing" and must be passed
on purpose). A unit test checks that every excluded name exists in `descList`.

### Failures

A value that cannot be computed or is not finite becomes **NaN**, never 0, and
a record `{"row", "feature", "reason"}` is appended to `df.attrs["failures"]`.
Inside a model, NaN features are imputed with the **training-fold median**
(below); NaN **targets** are not allowed: drop those rows explicitly before
evaluation (`nested_cv` raises).

## 2. Model families (`topoqspr.models`)

Every pipeline (fitted on training rows only):
`SimpleImputer(median)` → `DropConstantAndDuplicate` (drops training-constant
columns and exact duplicate columns, first kept) → `StandardScaler` (only
where scale matters) → model.

| Family | Model | Grid | Notes |
|---|---|---|---|
| `mean` | `DummyRegressor(mean)` | — | null model; used automatically for a 0-column set |
| `linear` | OLS | — | only if n_inner_train ≥ p + 2 (else skipped, reason recorded) |
| `ridge` | `Ridge` | α ∈ {0.01, 0.1, 1, 10, 100} | scaled |
| `elasticnet` | `ElasticNet(max_iter=20000)` | α ∈ {0.001, 0.01, 0.1, 1}, l1_ratio ∈ {0.2, 0.5, 0.8} | scaled; **y standardised** |
| `pls` | `PLSRegression` | n_components ∈ {1, 2, 3, 5, 10}, ≤ min(p, n − 1) | scaled |
| `knn` | `KNeighborsRegressor` | k ∈ {3, 5, 10} (≤ n), uniform/distance | scaled |
| `svr` | RBF `SVR(gamma="scale")` | C ∈ {0.3, 3, 30}, ε = 0.1 | scaled; **y standardised** |
| `rf` | `RandomForestRegressor(200, n_jobs=1)` | max_features ∈ {1.0, 0.33}, min_samples_leaf ∈ {1, 5} | |
| `extratrees` | `ExtraTreesRegressor(200, n_jobs=1)` | as rf | |
| `gbr` | `GradientBoostingRegressor(200)` | learning_rate ∈ {0.05, 0.1}, max_depth ∈ {2, 3} | |
| `xgb` | `XGBRegressor(200, hist, n_jobs=1)` | learning_rate ∈ {0.05, 0.1}, max_depth ∈ {3, 6} | |

Elastic net and SVR standardise y because their penalty and ε-tube are in units
of y (v0.x used the same C for MolWt ≈ 400 and logP ≈ 2). All estimators get
`random_state = 0`. **Every fit is on `sklearn.base.clone`** of a freshly built
estimator; registry prototypes are never fitted (v0.x shared one
`LinearRegression` instance across targets, so a stored model was silently
refitted on the next target).

Default family list (`DEFAULT_FAMILIES`, also the tie-break order): ridge,
elasticnet, pls, linear, knn, svr, rf, extratrees, gbr, xgb.

## 3. Repeated nested cross-validation (`evaluate.nested_cv`)

Defaults: K = 5 outer folds, R = 5 repeats with seeds 0–4, L = 3 inner folds.

For repeat r (seed s_r):

1. Outer split: `KFold(5, shuffle=True, random_state=s_r)`; with `groups`,
   `balanced_group_kfold(groups, 5, s_r)`: groups ordered by size
   (descending, ties in a random order from `default_rng(s_r)`), each placed
   in the fold with the fewest rows so far (Graham's LPT rule). This balances
   **rows**, not group counts. `GroupKFold(shuffle=True)` (scikit-learn ≥ 1.6)
   balances group counts, and on a Murcko-like distribution (one 24-row
   group plus 16 singletons, n = 40) gave 3-row test folds next to a 26-row
   one; LPT gives 24/4/4/4/4. A group larger than n/K still makes one large
   fold; `meta["min_n_test"]`, `meta["max_n_test"]` and the underpowered
   label (below) report it.
2. For outer fold k, on the outer-training rows **only**:
   inner split `KFold(3, shuffle=True, random_state=1000·s_r + k)` (grouped
   if groups are given). For every feasible family and grid point, the inner
   mean squared error MSE_inner = mean over inner folds of
   (1/n_te) Σ (y − ŷ)² is computed (`GridSearchCV(refit=False,
   error_score=nan)`). A grid point that fails in any inner fold is not
   selectable. Each family's best point is recorded; the family with the
   smallest MSE_inner is chosen (ties → first in the family list).
3. A clone of the chosen family with the chosen parameters is fitted on all
   outer-training rows and predicts the outer-test rows **once**.
4. If `ad` is set, an applicability domain is fitted on the outer-training
   rows and flags each outer-test row.

Selection uses MSE, not R², because with equal-size inner folds the mean MSE
ranks configurations exactly as pooled inner R² does, and it is defined even
when an inner test fold has constant y.

The procedure never tunes or selects on outer-test rows; a unit test
instruments a custom estimator and checks that no outer-test row reaches any
fit or inner predict, and that the final fit sees exactly the outer-training
rows (and, when grouped, that no group crosses the outer split).

### Metrics

For a set of predictions with mean ȳ of the true values:
R² = 1 − Σ(y − ŷ)² / Σ(y − ȳ)² (NaN if n < 2 or the denominator is 0),
MAE = mean |y − ŷ|, RMSE = √(mean (y − ŷ)²).

* **Fold metrics** use the outer-test rows of one fold (ȳ = test-fold mean).
* **Pooled per-repeat metrics** use all out-of-fold predictions of one repeat
  (every row once; ȳ = mean of all y). Pooled R² is less noisy than the mean
  of fold R² for small folds; both are reported.

### Output (`NestedCVResult`)

`predictions` — one row per (repeat, molecule):

| column | meaning |
|---|---|
| repeat, fold | outer repeat index (0..R−1) and fold index (0..K−1) in which the row was a test row |
| row | the row label from X's index |
| pos | 0-based row position |
| y_true, y_pred | observed value and out-of-fold prediction |
| in_domain | 1/0 applicability-domain flag for this row in this fold (NaN if `ad=None`, p = 0, or the AD is undefined for that training part) |

`folds` — one row per outer fold:

| column | meaning |
|---|---|
| repeat, fold, seed | split identity |
| n_train, n_test | outer-training / outer-test sizes |
| family, params | the family and hyperparameters chosen by the inner loop (JSON) |
| inner_mse | the winning inner-CV MSE |
| r2, mae, rmse | test-fold metrics |
| ad_coverage | fraction of test rows inside the AD (NaN if not fitted) |
| ad_reason | why the AD was not fitted (`ad=None`, no columns, or the `ValueError` text, e.g. leverage with n ≤ p′ + 1, all columns constant); empty if fitted |
| status, reason | `ok`, or `failed` with the reason (no family fit, refit/predict error, non-finite predictions). An `ok` fold also carries a reason when R² is undefined (`n_test < 2`, constant y_test) or when the null family was substituted (below) |

`repeats` — repeat, seed, n, pooled r2, mae, rmse.
`inner` — for every outer fold and family: status (`ok`/`skipped`/`failed`),
reason, best inner MSE, best params, number of grid points and failed points.
`summary()` — fold mean/SD/median of R², fold-mean MAE/RMSE, mean/SD over
repeats of pooled R² and MAE, `n_failed_folds`, `n_valid_folds` (finite R²),
how often each family was chosen, `underpowered`, and `reason`.
**No NaN is skipped**: if any fold failed or has an undefined R², the fold
statistics are NaN, and if any repeat's pooled R² is undefined the pooled
statistics are NaN; `reason` gives the counts. Averaging only the surviving
folds would select on success. `optimism()` and `nested_r2_statistic()`
inherit this (NaN plus reason), so B_train/B_CV are never computed on a
subset of folds. `meta` — n, p, columns, families, seeds, split settings,
`n_missing_cells` in X, `min_n_test`, `max_n_test`, `n_ad_failed`,
`underpowered`, `underpowered_reason`.

**Uninformative feature sets.** If, in an outer-training part, every column
is all-NaN or constant (nothing would survive imputation and
`DropConstantAndDuplicate`), the fold uses the null family `mean` and records
the reason, so such a set scores as the null model instead of failing.

**Underpowered label.** `meta["underpowered"]` is True when n < `min_n`
(default 50, EXPERIMENT_PLAN §3.3) or the smallest outer test fold has fewer
than `min_n_test` (default 10) rows. The run still completes; an
underpowered result is reported descriptively, with no ML claim.

## 4. Paired comparison of feature sets (`evaluate.compare_feature_sets`)

All sets are evaluated with identical outer folds (folds depend only on n,
seeds and groups; the test suite checks this). For set S against the
reference (default `size`), per outer fold j:

* d_R²_j = R²_S,j − R²_ref,j
* d_MAE_j = MAE_ref,j − MAE_S,j
* d_RMSE_j = RMSE_ref,j − RMSE_S,j

(positive = S is better). The headline ΔR²(TI − size) is the mean of d_R²
over all R·K folds with the corrected repeated-CV t interval (§5).

**Decision rule when the size model is exact.** If y is an exact function of
size (synthetic control), R²_size ≈ 1 and d_R² is a tiny negative number with
near-zero variance (ridge spreads weight over collinear TI columns), so the CI
can exclude 0 on the negative side (observed: d = −1·10⁻⁵, CI
[−1·10⁻⁵, −0], p = 0.001). "CI covers 0" is then the wrong criterion. The
rule that matters is **non-superiority**: the upper CI bound of d_R² must not
exceed a small margin (tested: ≤ 0.01). With noise in y (the realistic case),
"CI covers 0" holds and is tested separately.

`summary` columns: set, p (number of input columns), r2_fold_mean,
r2_pooled_mean, mae_fold_mean, n_failed_folds, n_valid_folds, underpowered,
reason, and for every non-reference
set d_r2_mean, d_r2_ci_low, d_r2_ci_high, d_r2_p, d_mae_mean, d_mae_ci_low,
d_mae_ci_high, d_mae_p, and d_r2_q (BH across the sets in this call). The
study-level BH correction across disease × target pairs (EXPERIMENT_PLAN §7)
is applied by the caller to the collected `d_r2_p` values.

`paired` columns: set, repeat, fold, d_r2, d_mae, d_rmse, n_train, n_test.

## 5. Statistics (`topoqspr.stats`)

### Corrected repeated k-fold t-test

For J = R·K paired differences d_j with training size n₁ and test size n₂:

d̄ = (1/J) Σ d_j, s² = (1/(J−1)) Σ (d_j − d̄)²,
se = √((1/J + n₂/n₁) · s²), t = d̄ / se, df = J − 1,
CI = d̄ ± t₁₋α/₂,J−1 · se, two-sided p from t(J − 1).

The n₂/n₁ term (Nadeau & Bengio 2003) corrects the variance for overlapping
training sets; Bouckaert & Frank (2004) showed that applied to repeated k-fold
CV with df = R·K − 1 it gives replicable tests. For K = 5, n₂/n₁ = 1/4, so the
interval is √(1 + J/4) ≈ 2.7 times wider than the naive one at J = 25. With
grouped (unequal) folds n₂/n₁ = mean(n_test)/mean(n_train). If s² = 0 the CI
collapses to d̄ and p = 0 (d̄ ≠ 0) or 1. NaN differences propagate by default.

### Benjamini–Hochberg

For m non-missing p-values sorted p₍₁₎ ≤ … ≤ p₍ₘ₎:
q₍ᵢ₎ = min_{j ≥ i} min(1, m·p₍ⱼ₎/j); reject iff q ≤ α (default 0.05).
NaN p-values stay NaN, are never rejected and by default do not count towards m;
`benjamini_hochberg(p, q, m=M)` uses a fixed family size M ≥ #finite p instead (an
untestable member then acts as p = 1; used by the confirmatory table, §11.2; the same
`m=` for `benjamini_yekutieli`, where c(M) uses M). Tested against the
15-p-value example of Benjamini & Hochberg (1995), which rejects four.

### Spearman correlation

ρ from `scipy.stats.spearmanr` (pairwise-complete rows); CI by Fisher z:
tanh(atanh ρ ± z₁₋α/₂ · se) with se = √(1.06/(n − 3)) (Fieller, Hartley &
Pearson 1957; default) or √((1 + ρ²/2)/(n − 3)) (Bonett & Wright 2000).
Undefined (n < 4, constant input) → NaN with a reason. Report ρ, CI, q (BH),
n — not stars.

### Y-scrambling (permutation test)

s_obs = score(X, y), s_b = score(X, π_b(y)) for B permutations from
`numpy.random.default_rng(seed)`; p = (1 + #{s_b ≥ s_obs}) / (1 + B)
(Phipson & Smyth 2010). The score function must run the **whole** pipeline,
model selection included (Rücker et al. 2007); `evaluate.nested_r2_statistic`
returns the mean pooled nested R². Default B = 100.

### Applicability domain

Fitted on outer-training rows: median imputation, constant columns dropped,
standardisation with training mean/SD.

* `knn` (default): d(x) = mean Euclidean distance to the k = 5 nearest
  training rows (leave-one-out for training rows); in domain iff
  d(x) ≤ 95th percentile of training d. Works for p ≥ n and collinear TIs.
* `leverage`: h(x) = x̃ᵀ (X̃ᵀX̃)⁺ x̃ with X̃ = [1, Z]; warning threshold
  h* = 3(p′ + 1)/n (Tropsha, Gramatica & Gombar 2003). Refused when
  n ≤ p′ + 1, because every training leverage is then 1.

Report performance inside vs outside the domain from `predictions.in_domain`.

## 6. v0.x reproduction and Experiment C (`legacy_protocol`, `optimism`)

`legacy_protocol(X, y)` reproduces `ml.run_ml_qspr` of v0.x: rows with any
non-finite feature or target are dropped; for each of the ten v0.x families
(`models.LEGACY_MODELS`, identical classes, grids and seed 42, checked against
`ml.py` by a test) `GridSearchCV(Pipeline([scaler, model]), cv=KFold(3,
shuffle, 42), scoring="r2")` is fitted on **all rows**, then
`cross_val_score` on `KFold(5, shuffle, 42)` gives the fold-mean R² and MAE,
and the family with the largest mean R² is reported. This is not nested: the
rows used to tune and to choose the family are later the test folds. Training
R² of the winning pipeline (refitted on, and predicting, all rows) is added.
Differences from v0.x: `n_jobs = 1` (identical seeded results) and every fit
is cloned (the v0.x shared-instance bug does not affect the CV numbers of a
single target, only the stored models).

`optimism(legacy, nested)`:
B_train = R²_train − R²_nested and B_CV = R²_CV(v0.x) − R²_nested, with
R²_nested the mean outer-fold R² (the same kind of estimate as the v0.x fold
mean); `*_pooled` variants use the mean pooled R². On pure-noise data the test
suite checks that nested R² ≤ 0 on average while B_CV > 0 and B_train > 0.2.

## 7. Groups (scaffold split)

`scaffold_groups(mols)` returns RDKit `MurckoScaffoldSmiles` (Bemis & Murcko
1996). Acyclic molecules have an empty scaffold; by default each gets its own
group (`acyclic="singleton"`), since they share no ring framework; use
`"shared"` to put them in one group. Pass the result as `groups=` to
`nested_cv` / `compare_feature_sets`; the inner loop is grouped too.

## 8. Determinism

All splitters and estimators are seeded; families run with `n_jobs = 1`;
results are sorted by (repeat, pos); parameters are serialised with sorted
keys. Running twice, or with `n_jobs > 1` over outer folds, gives identical
frames (tested).

## 9. Runtime

Measured on mpcrlab12 (16 cores, load average ≈ 22 from other jobs) for one
outer fold of a drug-sized set (n_train = 123, p = 200 rdkit2d-like columns),
inner 3-fold selection plus refit, `n_jobs = 1`: ridge 1.4 s, elastic net
6.4 s, PLS 2.0 s, OLS skipped (p ≥ n), kNN 3.2 s, SVR 1.2 s, RF 86 s,
extra trees 49 s, GBR 90 s, XGBoost 98 s: ≈ 340 s per outer fold, so
≈ 2.4 h per 5 × 5 nested CV per feature set, target and disease with all ten
families (an independent run under load ≈ 25 measured ≈ 420 s per fold).
The four tree ensembles take ≈ 95 % of it. Use `n_jobs` > 1 (outer folds in
parallel; results identical) and budget the full grid before freezing the
plan; one full-family real-data smoke run is required before freezing.

## 10. Defaults at a glance

| Setting | Default |
|---|---|
| outer CV | 5-fold × 5 repeats, seeds 0–4 (grouped: row-balanced LPT) |
| inner CV | 3-fold, seed 1000·s_r + k |
| selection criterion | inner mean MSE, ties → family order |
| estimator seed | 0 |
| AD | kNN, k = 5, 95th percentile |
| CI | corrected repeated-CV t, α = 0.05 |
| multiple testing | BH, q = 0.05 |
| `stats.y_scrambling` (stand-alone helper) | B = 100 |
| protocol stage (§11; only with `config.protocol`) | |
| Y-randomised control | K = 200 permutations, seed 0, 95th percentile of the null |
| whole-procedure paired bootstrap | B = 200 resamples, seed 0; planned for the pooled set only (EXPERIMENT_PLAN §7); at `v1.0-dev` enabled by default for every cohort, so per-disease configs set `bootstrap.enabled = false` (PIPELINE.md §9) |
| practical threshold | ΔR² ≥ 0.05 |
| confirmatory FDR | BH (BY sensitivity), q = 0.05, fixed m |
| underpowered | n < 50 or min n_test < 10 |

## References

* Nadeau C, Bengio Y (2003) Inference for the generalization error. *Machine Learning* 52:239–281. doi:10.1023/A:1024068626366
* Bouckaert RR, Frank E (2004) Evaluating the replicability of significance tests for comparing learning algorithms. *PAKDD 2004*, LNCS 3056:3–12. doi:10.1007/978-3-540-24775-3_3
* Benjamini Y, Hochberg Y (1995) Controlling the false discovery rate. *J R Stat Soc B* 57:289–300. doi:10.1111/j.2517-6161.1995.tb02031.x
* Fieller EC, Hartley HO, Pearson ES (1957) Tests for rank correlation coefficients. I. *Biometrika* 44:470–481. doi:10.1093/biomet/44.3-4.470
* Bonett DG, Wright TA (2000) Sample size requirements for estimating Pearson, Kendall and Spearman correlations. *Psychometrika* 65:23–28. doi:10.1007/BF02294183
* Phipson B, Smyth GK (2010) Permutation p-values should never be zero. *Stat Appl Genet Mol Biol* 9:39. doi:10.2202/1544-6115.1585
* Rücker C, Rücker G, Meringer M (2007) y-Randomization and its variants in QSPR/QSAR. *J Chem Inf Model* 47:2345–2357. doi:10.1021/ci700157b
* Tropsha A, Gramatica P, Gombar VK (2003) The importance of being earnest. *QSAR Comb Sci* 22:69–77. doi:10.1002/qsar.200390007
* Tropsha A (2010) Best practices for QSAR model development, validation, and exploitation. *Mol Inform* 29:476–488. doi:10.1002/minf.201000061
* Varma S, Simon R (2006) Bias in error estimation when using cross-validation for model selection. *BMC Bioinformatics* 7:91. doi:10.1186/1471-2105-7-91
* Cawley GC, Talbot NLC (2010) On over-fitting in model selection and subsequent selection bias in performance evaluation. *J Mach Learn Res* 11:2079–2107.
* Krstajic D, Buturovic LJ, Leahy DE, Thomas S (2014) Cross-validation pitfalls when selecting and assessing regression and classification models. *J Cheminform* 6:10. doi:10.1186/1758-2946-6-10
* Alexander DLJ, Tropsha A, Winkler DA (2015) Beware of R²: simple, unambiguous assessment of the prediction accuracy of QSAR and QSPR models. *J Chem Inf Model* 55:1316–1322. doi:10.1021/acs.jcim.5b00206
* Graham RL (1969) Bounds on multiprocessing timing anomalies. *SIAM J Appl Math* 17:416–429. doi:10.1137/0117039
* Bemis GW, Murcko MA (1996) The properties of known drugs. 1. Molecular frameworks. *J Med Chem* 39:2887–2893. doi:10.1021/jm9602928
* Wildman SA, Crippen GM (1999) Prediction of physicochemical parameters by atomic contributions. *J Chem Inf Comput Sci* 39:868–873. doi:10.1021/ci990307l
* Ertl P, Rohde B, Selzer P (2000) Fast calculation of molecular polar surface area. *J Med Chem* 43:3714–3717. doi:10.1021/jm000942e
* Bickerton GR et al. (2012) Quantifying the chemical beauty of drugs. *Nat Chem* 4:90–98. doi:10.1038/nchem.1243

Every DOI above was checked against Crossref on 2026-10-02 (Cawley & Talbot has no DOI; JMLR volume and pages from the journal).


## Benjamini–Yekutieli (sensitivity analysis)

`stats.benjamini_yekutieli(pvals, q)` applies the BY step-up procedure, valid under arbitrary dependence (Benjamini & Yekutieli 2001): BH applied to c(m)·p with c(m) = Σ_{i=1..m} 1/i. It is reported next to BH for the confirmatory family because those endpoints share drugs and CV folds (EXPERIMENT_PLAN §7). Verified against a brute-force implementation of the definition (`test_benjamini_yekutieli_matches_definition`).

## Sanity check: pure-noise targets (synthetic, run 2026-10-03)

n = 80 molecules, p = 60 i.i.d. standard-normal features, y independent standard normal; families
ridge + random forest; one outer repeat (seed 0). Script: `tools/noise_sanity.py`.

| Data seed | Nested CV R² (pooled) | Nested CV R² (fold mean) | v0.x protocol CV R² | v0.x training R² |
|---|---|---|---|---|
| 0 | −0.16 | −0.23 | −0.30 | 0.84 |
| 1 | −0.06 | −0.09 | −0.25 | 0.85 |
| 2 | −0.08 | −0.13 | −0.22 | 0.84 |

Nested CV correctly reports no signal. The v0.x cross-validated score is not optimistic in this
setting (its tuning leak is small relative to split-to-split variance, as the audit also found), but
the **in-sample R² of a random forest on pure noise is about 0.85**: in-sample fit statistics, as
reported by most published disease-specific TI-QSPR studies, carry no evidence of predictive power.

## 11. Confirmatory protocol (pipeline protocol stage)

Implements the remaining parts of EXPERIMENT_PLAN §7 inside `pipeline.run_study`
(configuration and output files: PIPELINE.md §9). Code: `pipeline.permutation_control`,
`pipeline.paired_bootstrap`, `pipeline.ad_report`, `pipeline.confirmatory_table`;
`stats.permutation_pvalue`, `stats.percentile_interval`, `stats.group_bootstrap`.
Tests: `tests/test_protocol.py` (synthetic data only). **The plan is not frozen**: the
family, threshold, K and B below are defaults that the frozen plan must confirm.

### 11.1 Confirmatory family

`protocol.family` lists (target, set, reference) triples, or `"auto"` = every class-3
target × (TI vs `evaluation.reference`), the primary endpoint of EXPERIMENT_PLAN §7. The
reference must be `evaluation.reference` (paired differences exist only against it). The
table uses the results of one learner, `protocol.learner` (default `selected`). The expanded
family, the learner, the threshold and the rules below are written to the run record
(`protocol`).

### 11.2 Confirmatory table (`confirmatory.csv`)

One row per family member. From the main analysis of that learner: n, `d_r2_mean`
(mean of the R·K paired fold differences R²_set − R²_ref), its corrected repeated-CV t CI
and p (§5), and the same for d_MAE = MAE_ref − MAE_set. Then, **across the family**.

**Family membership and m** (`pipeline.FAMILY_RULE`, written to the run record as
`protocol.family_membership`). EXPERIMENT_PLAN §7 restricts the confirmatory family to
targets with ≥ 50 measured values, and an underpowered result carries no ML claim (§3).
So a pre-specified member is left out of the multiplicity family (`in_family` = False) when
its result is underpowered (n < `evaluation.min_n` or smallest outer test fold <
`evaluation.min_n_test`) or its target was not evaluated and has fewer than `min_n` values;
it keeps its d_R² estimate and CI (descriptive) but gets NaN q, `informative_*` = None,
`practically_significant` = None and the reason "underpowered: not in confirmatory
family". Membership depends only on n and the fold sizes, never on an outcome. Every other
member counts in m = `family_m`, **including** a member whose p is undefined (e.g. an
undefined fold R²) or that has no result: it is never rejected (equivalent to p = 1 in the
step-up), keeps NaN q and says why in `reason`. An evaluation failure therefore cannot
shrink m and loosen the correction for the other members (`stats.benjamini_hochberg(...,
m=)`, `stats.benjamini_yekutieli(..., m=)`).

| column | definition |
|---|---|
| `d_r2_q_bh` | `stats.benjamini_hochberg(d_r2_p)` (§5) |
| `d_r2_q_by` | `stats.benjamini_yekutieli(d_r2_p)` (sensitivity under arbitrary dependence: the endpoints share drugs and folds) |
| `informative_bh` / `_by` | q ≤ `protocol.q` (0.05) **and** d_r2_mean > 0: the CI excludes 0 on the positive side after FDR control (EXPERIMENT_PLAN §7) |
| `practically_significant` | d_r2_mean ≥ `protocol.practical_threshold` (default 0.05, the plan's minimal meaningful gain; to be confirmed by the co-authors) |
| `ad_share_outside_set` / `_reference` | mean over outer folds of the share of outer-test molecules outside the kNN AD (§11.4) |
| `perm_k`, `perm_p_set`, `perm_uninformative_set` (+ `_reference`) | §11.3 |
| `boot_*` | §11.5 |
| `reason` | everything that makes a value missing or qualifies it (target not evaluated, failed folds, underpowered, failed bootstrap resamples) |

A row whose p is NaN keeps NaN q values and None flags and counts towards m (see above).
A missing permutation, bootstrap or AD value is explained in `reason` (component disabled,
set not in `protocol.permutation.feature_sets`, AD not fitted, failed permutations, ...),
and the run record says `enabled: false` for a disabled permutation or bootstrap.

### 11.3 Y-randomised control (`permutation.csv`, `permutation_summary.csv`)

For every target of the family and every set in it (both the set and the reference, or
`protocol.permutation.feature_sets`), the target is permuted K times (default **K = 200**,
EXPERIMENT_PLAN §7; Rücker et al. 2007) and the **whole** nested procedure is rerun with
the learner's families, the same outer seeds, inner rule and groups:

* permutation b of a target uses `numpy.random.default_rng([seed, b])`
  (`pipeline.permutation_indices`), the same for every set of that target, so sets are
  compared on identical nulls. Without CV groups it is `rng.permutation(n)`. **With CV
  groups** (connectivity groups in the random split, scaffold groups in the scaffold
  split) it is a restricted (block) permutation: whole groups' y blocks are exchanged
  among groups of equal size, and the rows of singleton groups and of groups whose size is
  unique are permuted row-wise among themselves. A row-wise permutation would break the
  within-group dependence of y that the real data have, and the test would be conservative:
  on a synthetic null with 18 pairs (identical X within a pair, y = u_g + 0.1ε, y
  independent of X; OLS, 3 outer folds, K = 19, 60 data sets) the row-wise mean p was
  0.737 (P(p ≤ 0.25) = 0.10) against 0.525 expected; the block permutation gives 0.544
  (P(p ≤ 0.25) = 0.23). The run record gives the scheme per target;
* X, its rows and the CV group labels are not moved, so outer and inner folds are those of
  the real analysis (they depend only on n, seeds and groups);
* statistic: mean over repeats of the pooled out-of-fold R² (NaN if any outer fold
  failed); the observed value is the main analysis's `r2_pooled_mean` of the same learner
  (`protocol.permutation.learner` may name another learner from `evaluation.learners`, to
  save time; the observed value then comes from that learner);
* p = (1 + #{b: s_b ≥ s_obs}) / (1 + K_valid) (Phipson & Smyth 2010); failed permutations
  are counted (`n_failed`) and excluded from K_valid;
* `uninformative` = s_obs ≤ the 95th percentile (`numpy.quantile`, linear) of the valid
  null values: the plan's rule for reporting a feature set as uninformative.

Calibration is tested two ways: `stats.permutation_pvalue` on exchangeable draws
(K = 19, 4,000 replicates: P(p ≤ a) within 0.025 of a for a ∈ {0.05, 0.1, 0.25, 0.5}),
and the nested procedure itself on six pure-noise data sets (OLS, K = 9: mean p in
[0.25, 0.85], at most two p ≤ 0.1) and on a strong-signal data set (p = 1/(K + 1),
not uninformative).

### 11.4 Applicability-domain coverage (`ad.csv`, `ad_summary.csv`)

Written for every run with `evaluation.ad` set (every target, learner and set, not only
the family). Per outer fold: number of outer-test molecules inside / outside the kNN AD
fitted on that fold's training part (§5), share outside, and R²/MAE on the in-domain and
on the outside test rows. Per (target, learner, set): mean and max share outside over
folds, pooled share, and the mean over repeats of the pooled per-repeat R² and MAE on all
rows (primary), on in-domain rows and on outside rows (**secondary**, EXPERIMENT_PLAN §7).
An undefined metric (fewer than two rows, constant y) is NaN, never 0; a fold without an
AD (the zero-column `null` set) makes the summary shares NaN with the reason.

### 11.5 Whole-procedure paired bootstrap (`bootstrap.csv`)

Sensitivity analysis for the CI of the primary ΔR²(TI − size) and ΔMAE (CV intervals can
undercover; Bates, Hastie & Tibshirani 2023), with B resamples (default B = 200).
The plan restricts it to the **pooled** set (EXPERIMENT_PLAN §7): at per-disease n the
resampled test folds are too small for R² to be defined (PIPELINE.md §10). At `v1.0-dev`
(e35dc4a) `protocol.bootstrap.enabled` defaults to true for every cohort, so the frozen
per-disease configs must set it to false; a code change in preparation (branch `fx-code`)
makes it default to on only for `cohort.source = "pooled"`. For b = 0..B−1:

1. **Resample units, not rows.** The units are the connectivity groups (stereoisomers,
   racemate/enantiomer, isotopologues; PIPELINE.md §1) and every other molecule on its
   own; G units are drawn with replacement with `default_rng([seed, b])`
   (`stats.group_bootstrap`), taking all rows of each drawn unit.
2. **Rerun the identical procedure** on the resample: repeated nested CV of the set and
   the reference with the same outer seeds, inner rule and families. Every copy of an
   original molecule keeps its original CV group label (the main analysis's group, e.g.
   scaffold, or the unit), so copies and partners never straddle a fold boundary
   (a duplicated molecule in training and test would inflate both R²).
3. Record d_R² = fold mean of R²_set − R²_ref and d_MAE likewise (the estimands of §11.2).

Row `resample = −1` reruns the procedure on the original data with the same families
(equal to the main analysis when the families are the confirmatory learner's). The
interval is the percentile interval [Q(α/2), Q(1 − α/2)] (`numpy.quantile`, linear) with
α = `evaluation.alpha`; **if any resample failed, the interval is NaN with the reason**
rather than computed from the successful resamples (that would condition on success).

**Cost and the fast learner subset.** One resample costs two full nested CVs, so B = 200
costs 400 nested CVs per family row. `protocol.bootstrap.families` may name a faster
subset (e.g. `["ridge", "pls", "knn"]` or `["ridge"]`) than the confirmatory learner; the
run record then labels the families as a "fast subset, cost-aware sensitivity", and the
bootstrap describes the variability of that cheaper procedure, not of the confirmatory
learner (`boot_d_r2_original` gives its original-data estimate next to `d_r2_mean`).

Tested on synthetic data with a known effect (y = s + t + ε, set {s, t} vs reference {s},
population ΔR² = 1/3, n = 120, OLS, B = 15): the percentile interval covers 1/3 and
excludes 0; for an added pure-noise column it covers 0.

### 11.6 Determinism

All randomness is seeded (`default_rng([seed, b])` per permutation and resample; the nested
CV is deterministic, §8). Permutations and resamples run in parallel with
`runtime.n_jobs` and give identical rows. A re-run of a protocol run record reproduces
`confirmatory.csv`, `permutation*.csv`, `bootstrap.csv` and the AD files byte for byte
(`test_protocol_is_deterministic`).

* Bates S, Hastie T, Tibshirani R (2023, online; 2024, issue 546) Cross-validation: what does it estimate and how well does it do it? *J Am Stat Assoc* 119:1434–1445. doi:10.1080/01621459.2023.2197686
* Benjamini Y, Yekutieli D (2001) The control of the false discovery rate in multiple testing under dependency. *Ann Stat* 29:1165–1188. doi:10.1214/aos/1013699998
References of §11 checked against Crossref on 2026-10-03.
