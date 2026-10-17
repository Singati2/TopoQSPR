"""Unit tests for the Colab reproducibility modules in docs/research/scripts/ (no network, no data)."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "docs/research/scripts"))
import fixed_learner as FL  # noqa: E402
import patent_census as PC  # noqa: E402
import repro_common as RC  # noqa: E402
import reproduce_e3 as R3  # noqa: E402


def test_classify_categories():
    d = FL.DELTA
    assert FL.classify(0.03, 0.08, 0.025, 0.07) == "superior_by_more_than_delta"
    # positive 95% CI and point estimate >= delta is NOT a gain of at least delta
    assert FL.classify(0.005, 0.06, 0.01, 0.055) == "superior_not_shown_beyond_delta"
    assert FL.classify(-0.01, 0.01, -0.008, 0.008) == "equivalent_within_delta"
    assert FL.classify(-0.03, 0.03, -0.025, 0.025) == "inconclusive"
    assert FL.classify(-0.08, -0.03, -0.07, -0.025) == "inferior_by_more_than_delta"
    assert FL.classify(-0.03, -0.001, -0.028, -0.004) == "inferior_not_shown_beyond_delta"
    assert d == 0.02


def test_units_full_count_and_unique_ids():
    U = FL.units("full")
    n_int = 2 * len(FL.FAMILIES) * len(FL.ARMS) * len(FL.SEEDS_INTERNAL) * 5
    n_ext = len(FL.ARMS) * (len(FL.FAMILIES) - 1 + len(FL.TRAIN_SEEDS_EXTERNAL["extratrees"]))
    assert len(U) == n_int + n_ext == 760
    ids = [FL.unit_id(u) for u in U]
    assert len(set(ids)) == len(ids)
    # sharding partitions the units exactly
    parts = [{FL.unit_id(u) for j, u in enumerate(U) if j % 4 == i} for i in range(4)]
    assert set().union(*parts) == set(ids) and sum(map(len, parts)) == len(ids)


def test_smoke_and_bench_units_never_full_size():
    assert 0 < len(FL.units("smoke")) < len(FL.units("bench")) < len(FL.units("full"))
    assert {u["context"] for u in FL.units("bench")} == {"mp_E3b", "external"}


def test_arm_frames_union_order():
    C = pd.DataFrame({"comp.C": [1., 2.], "comp.H": [3., 4.]})
    E = pd.concat([pd.DataFrame({"e12": [1., 0.], "maxdeg": [3., 4.]}), C], axis=1)
    T = pd.DataFrame({"ti.W": [5., 6.]})
    arms = FL.arm_frames({"composition": C, "edgetypes_composition": E, "TIa": T,
                          "TIa_composition": pd.concat([T, C], axis=1)})
    assert list(arms["C+E+T"].columns) == ["e12", "maxdeg", "ti.W", "comp.C", "comp.H"]
    assert list(arms["C+T"].columns) == ["ti.W", "comp.C", "comp.H"]


def test_bh_matches_reference():
    p = np.array([0.01, 0.04, 0.03, 0.2])
    # sorted 0.01, 0.03, 0.04, 0.2 -> 0.04, min(0.06, 0.0533), 0.0533, 0.2
    np.testing.assert_allclose(RC.bh(p), [0.04, 0.16 / 3, 0.16 / 3, 0.2])


def test_corrected_t_matches_package():
    from topoqspr import stats
    rng = np.random.default_rng(0)
    d = rng.normal(0.05, 0.02, 15); n1 = np.full(15, 800); n2 = np.full(15, 200)
    m, lo, hi, p = R3.corrected_t(d, n1, n2)
    ref = stats.corrected_ttest(pd.Series(d), pd.Series(n1), pd.Series(n2), alpha=0.05)
    assert np.allclose([m, lo, hi, p], [ref["mean"], ref["ci_low"], ref["ci_high"], ref["p"]])


SDF = """mol-1
  RDKit          2D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0
    1.5000    0.0000    0.0000 O   0  0
  1  2  1  0
M  END
> <MeltingPoint>
-114.1

> <Source>
patent US0000000

$$$$
"""


def test_census_molblocks_drop_every_data_field(tmp_path):
    f = tmp_path / "x.sdf"
    f.write_text(SDF * 2)
    recs = list(PC.molblocks(f))
    assert len(recs) == 2
    for i, title, block in recs:
        assert title == "mol-1"
        assert block.rstrip().endswith("M  END")
        assert "MeltingPoint" not in block and "-114.1" not in block and "patent" not in block
    from rdkit import Chem
    assert Chem.MolToSmiles(Chem.MolFromMolBlock(recs[0][2])) == "CO"


def test_checks_fail_closed(tmp_path):
    ck = RC.Checks("t")
    ck.require("ok", True)
    ck.close("near", 1.0, 1.0 + 1e-12, 1e-9)
    ck.require("bad", False, 1, 2)
    with pytest.raises(SystemExit):
        ck.finish(tmp_path)
    assert (tmp_path / "checks_t.csv").exists()


def test_holm():
    np.testing.assert_allclose(FL.holm([0.01, 0.04, 0.03]), [0.03, 0.06, 0.06])
