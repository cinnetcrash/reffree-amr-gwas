"""Regression tests for the helper scripts. Pure-Python deps only."""
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
EXAMPLE = ROOT / "example"

from newick import parse                       # noqa: E402
from mobility_score import sankoff              # noqa: E402
from kinship_from_tree import mrca_depth_matrix  # noqa: E402

TREE = ("((c1a:0.01,c1b:0.01,c1c:0.01,c1d:0.01):0.20,"
        "(c2a:0.01,c2b:0.01,c2c:0.01,c2d:0.01):0.20,"
        "(c3a:0.01,c3b:0.01,c3c:0.01,c3d:0.01):0.20);")


def test_kinship_is_phylogenetic_covariance():
    names, K = mrca_depth_matrix(parse(TREE))
    assert np.allclose(K, K.T)                      # symmetric
    assert np.allclose(np.diag(K), 0.21)            # root-to-tip on the diagonal
    i, j, k = names.index("c1a"), names.index("c1b"), names.index("c2a")
    assert K[i, j] > K[i, k]                         # within-clade shares more ancestry
    assert abs(K[i, k]) < 1e-9                       # cross-clade shares only the root


def test_sankoff_distinguishes_clonal_from_mobile():
    root = parse(TREE)
    assert sankoff(root, {"c2a", "c2b", "c2c", "c2d"})[1] == 1   # one full clade
    assert sankoff(root, {"c2a", "c2b", "c2c"})[1] == 1          # clade minus a loss
    assert sankoff(root, {"c1b", "c2a", "c2d", "c3a", "c3d"})[1] == 5  # scattered
    assert sankoff(root, {"c1a"})[1] == 1                        # singleton


def test_patterns_threshold(tmp_path):
    f = tmp_path / "pat.txt"
    f.write_text("a\na\nb\nc\nb\n")          # 3 unique patterns
    out = subprocess.check_output(
        [sys.executable, str(SCRIPTS / "patterns_threshold.py"),
         str(f), "--print-threshold"]).decode().strip()
    assert abs(float(out) - 0.05 / 3) < 1e-9


def test_triage_classifies_known_novel_and_lineage_suspect(tmp_path):
    import pandas as pd
    mob = tmp_path / "mobility.tsv"
    tbl = tmp_path / "candidates.tsv"
    fa = tmp_path / "novel.fasta"
    subprocess.check_call(
        [sys.executable, str(SCRIPTS / "mobility_score.py"),
         "--tree", str(EXAMPLE / "core.tree"),
         "--pyseer", str(EXAMPLE / "unitigs.pyseer"),
         "--significant", str(EXAMPLE / "significant.txt"), "--out", str(mob)])
    subprocess.check_call(
        [sys.executable, str(SCRIPTS / "triage_novelty.py"),
         "--significant", str(EXAMPLE / "significant.txt"),
         "--blast", str(EXAMPLE / "blast6.tsv"), "--mobility", str(mob),
         "--out-table", str(tbl), "--out-fasta", str(fa),
         "--min-identity", "80", "--min-coverage", "0.8"])
    d = pd.read_csv(tbl, sep="\t")
    assert (d["classification"] == "KNOWN_AMR").sum() == 1
    novel = d[d["classification"] == "NOVEL_CANDIDATE"]
    assert len(novel) == 2
    assert novel["lineage_suspect"].sum() == 1           # the clade-tracking one
    mobile = novel[~novel["lineage_suspect"]]
    assert (mobile["independent_gains"] >= 3).all()      # scattered = mobile
