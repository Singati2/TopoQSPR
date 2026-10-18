#!/bin/bash
# doc 14 r2 full run, frozen at 1457bc9 (owner approval 2026-10-18)
set -euo pipefail
R=${RUN:-$HOME/topoqspr_runs}; C=$R/code_1457bc9; PY="nice -n 10 $R/venv/bin/python -I"
export TOPOQSPR_CACHE=$R/cache OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
OUT=$R/out_1457bc9; mkdir -p $OUT
cd $C
echo "== restore $(date -Is)"; $PY docs/research/scripts/restore_opera_values.py docs/research/data/E3 $R/e3_restored
echo "== audit $(date -Is)";   $PY docs/research/scripts/audit_environment.py docs/research/data/E3 $R/e3_restored $OUT/00_audit
echo "== build $(date -Is)";   $PY docs/research/scripts/fixed_learner.py build $OUT/fixed_learner
echo 1457bc9358fad88e449a13a218d60d1f0f985077 > $OUT/fixed_learner/APPROVED_FULL_RUN
echo "== run $(date -Is)";     $PY docs/research/scripts/fixed_learner.py run $OUT/fixed_learner --mode full --n-jobs 8
$PY docs/research/scripts/fixed_learner.py status $OUT/fixed_learner --mode full
echo "== analyse $(date -Is)"; $PY docs/research/scripts/fixed_learner.py analyse $OUT/fixed_learner --mode full --restored $R/e3_restored
echo "== robustness $(date -Is)"; $PY docs/research/scripts/robustness.py $R/e3_restored $OUT/fixed_learner $OUT/04_robustness --mode full --n-jobs 8
echo "== DONE $(date -Is)"
