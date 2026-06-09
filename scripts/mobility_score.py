#!/usr/bin/env python3
"""
mobility_score.py — turn significant unitigs into a MOBILOME read.

For each significant unitig we take the set of genomes that carry it (from the
unitig-caller .pyseer file) and map presence/absence onto the core-genome
phylogeny. Small-parsimony (Sankoff) then gives the minimum number of state
changes and, assuming absence is ancestral, the number of INDEPENDENT
ACQUISITIONS.

Interpretation
--------------
  independent_gains = 1            -> single acquisition, clonally inherited = VERTICAL
  independent_gains >> 1           -> acquired separately across the tree = MOBILE / horizontal
  gains_per_carrier -> 1.0         -> every carrier is its own acquisition = maximally mobile
  consistency_index (1/changes)    -> low = high homoplasy

Two uses:
  (1) Surveillance: high-gain unitigs flag a determinant jumping hosts, which
      clone-based surveillance never surfaces.
  (2) Confounding filter: a unitig that is GWAS-significant but tracks one clade
      perfectly (gains = 1, high tree-consistency) is a candidate residual
      lineage artifact, not a mobile causal element.

Usage:
    python mobility_score.py --tree core.tree --pyseer unitigs.pyseer \
        --significant significant.txt --out mobility.tsv
"""
import argparse
import hashlib
import sys
from newick import parse


def uid_of(seq):
    return hashlib.md5(seq.encode()).hexdigest()[:12]


def read_significant_seqs(path):
    seqs = []
    with open(path) as fh:
        header = fh.readline()  # variant af filter-pvalue lrt-pvalue beta ...
        for line in fh:
            if not line.strip():
                continue
            seqs.append(line.split("\t")[0].strip())
    return seqs


def read_carriers(pyseer_path, wanted):
    """Map unitig sequence -> set(carrier genome names), for wanted sequences."""
    wanted = set(wanted)
    carriers = {}
    with open(pyseer_path) as fh:
        for line in fh:
            if "|" not in line:
                continue
            left, right = line.rstrip("\n").split("|", 1)
            seq = left.strip()
            if seq not in wanted:
                continue
            names = {tok.split(":")[0] for tok in right.split() if tok}
            carriers[seq] = names
    return carriers


INF = float("inf")


def sankoff(root, present):
    """Binary small-parsimony via Sankoff DP (unit cost). Correct for
    multifurcating trees, unlike the Fitch shortcut. `present` = carrier leaves.
    Returns (parsimony_changes, independent_gains).

    Backtrack assumes ancestral ABSENCE (root prefers 0) and, on ties, keeps the
    parent state — so a single clade acquisition with subsequent losses is scored
    as ONE gain, not many. Gains = count of 0->1 edges in that reconstruction."""
    cost = {}   # node -> [cost_if_absent, cost_if_present]
    for n in root.postorder():
        if n.is_leaf():
            obs = 1 if n.name in present else 0
            cost[n] = [0.0 if obs == 0 else INF, 0.0 if obs == 1 else INF]
        else:
            c = [0.0, 0.0]
            for s in (0, 1):
                tot = 0.0
                for ch in n.children:
                    tot += min(cost[ch][sp] + (0 if sp == s else 1) for sp in (0, 1))
                c[s] = tot
            cost[n] = c
    changes = int(min(cost[root]))
    # backtrack
    assign = {}
    gains = 0
    for n in root.preorder():
        if n.parent is None:
            assign[n] = 0 if cost[n][0] <= cost[n][1] else 1
        else:
            ps = assign[n.parent]
            # choose child state minimising (subtree cost + transition); tie -> ps
            opts = [(cost[n][s] + (0 if s == ps else 1), 0 if s == ps else 1, s)
                    for s in (0, 1)]
            opts.sort(key=lambda t: (t[0], t[1]))
            assign[n] = opts[0][2]
            if ps == 0 and assign[n] == 1:
                gains += 1
    return changes, gains


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", required=True)
    ap.add_argument("--pyseer", required=True, help="unitig-caller .pyseer file")
    ap.add_argument("--significant", required=True,
                    help="pyseer association output filtered to significant hits")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    root = parse(open(args.tree).read())
    tree_tips = {lf.name for lf in root.leaves()}

    sig_seqs = read_significant_seqs(args.significant)
    carriers = read_carriers(args.pyseer, sig_seqs)

    missing = 0
    rows = []
    for seq in sig_seqs:
        present_all = carriers.get(seq, set())
        present = present_all & tree_tips
        missing += len(present_all - tree_tips)
        nc = len(present)
        if nc == 0 or nc == len(tree_tips):
            changes, gains, ci, gpc = 0, 0, float("nan"), float("nan")
        else:
            changes, gains = sankoff(root, present)
            ci = 1.0 / changes if changes else float("nan")
            gpc = gains / nc
        rows.append((uid_of(seq), nc, gains, gpc, changes, ci, seq))

    rows.sort(key=lambda r: (-(r[2] if r[2] == r[2] else -1), -r[1]))
    with open(args.out, "w") as out:
        out.write("uid\tn_carriers\tindependent_gains\tgains_per_carrier\t"
                  "parsimony_changes\tconsistency_index\tunitig_seq\n")
        for uid, nc, gains, gpc, changes, ci, seq in rows:
            out.write(f"{uid}\t{nc}\t{gains}\t{gpc:.3f}\t{changes}\t"
                      f"{ci:.3f}\t{seq}\n")

    if missing:
        sys.stderr.write(f"[mobility] {missing} carrier-genome labels not found "
                         f"in tree tips (ignored).\n")
    sys.stderr.write(f"[mobility] scored {len(rows)} significant unitigs "
                     f"against {len(tree_tips)} tree tips.\n")


if __name__ == "__main__":
    main()
