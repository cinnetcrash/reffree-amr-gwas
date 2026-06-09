#!/usr/bin/env python3
"""
kinship_from_tree.py — build a pyseer --similarity (LMM kinship) matrix from a
core-genome phylogeny, without needing pyseer's bundled scripts.

The kinship used by pyseer's LMM is the phylogenetic covariance under a
Brownian model: K_ij = shared evolutionary path length from the root to the
most recent common ancestor of tips i and j (K_ii = root-to-tip distance).
This is exactly what `scripts/phylogeny_distance.py --lmm` produces; this is a
dependency-light, auditable reimplementation.

Usage:
    python kinship_from_tree.py core_genome.tree > phylogeny_K.tsv

Feed the output to pyseer with:  --similarity phylogeny_K.tsv  --lmm

Note: pyseer expects a MIDPOINT-ROOTED tree for the LMM kinship. Root it first
(e.g. with `gotree reroot midpoint` or ete3) if your tree is unrooted.
"""
import sys
import numpy as np
from newick import parse, root_to_node_depths


def mrca_depth_matrix(root):
    leaves = list(root.leaves())
    names = [lf.name for lf in leaves]
    depth = root_to_node_depths(root)

    # path from each leaf to root (as a set of nodes, ordered)
    paths = {}
    for lf in leaves:
        chain, n = [], lf
        while n is not None:
            chain.append(n)
            n = n.parent
        paths[lf] = chain  # leaf -> ... -> root

    n = len(leaves)
    K = np.zeros((n, n))
    for a in range(n):
        set_a = set(paths[leaves[a]])
        for b in range(a, n):
            # MRCA = first node on leaf_b's path that is also on leaf_a's path
            mrca = next(node for node in paths[leaves[b]] if node in set_a)
            val = depth[mrca]
            K[a, b] = K[b, a] = val
    return names, K


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: kinship_from_tree.py TREE.nwk > K.tsv")
    root = parse(open(sys.argv[1]).read())
    names, K = mrca_depth_matrix(root)
    # write tab-separated matrix with row/column labels (pandas index_col=0 format)
    out = sys.stdout
    out.write("\t" + "\t".join(names) + "\n")
    for i, name in enumerate(names):
        out.write(name + "\t" + "\t".join(f"{v:.6g}" for v in K[i]) + "\n")


if __name__ == "__main__":
    main()
