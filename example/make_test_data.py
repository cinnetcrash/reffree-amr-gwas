#!/usr/bin/env python3
"""Generate small, realistically-formatted inputs to exercise the helper scripts.

Run:  python example/make_test_data.py
Writes the example inputs next to this file.
"""
import hashlib
import os
import random

random.seed(3)
HERE = os.path.dirname(os.path.abspath(__file__))
p = lambda name: os.path.join(HERE, name)
uid = lambda s: hashlib.md5(s.encode()).hexdigest()[:12]
rseq = lambda n: "".join(random.choice("ACGT") for _ in range(n))

# --- core phylogeny: 3 clades x 4 tips ---
def clade(prefix):
    tips = ",".join(f"{prefix}{x}:0.01" for x in "abcd")
    return f"({tips}):0.20"
tree = f"({clade('c1')},{clade('c2')},{clade('c3')});"
open(p("core.tree"), "w").write(tree)

# --- unitigs with known carriage patterns ---
unitigs = {
    "KNOWN_mobile":   (rseq(40), ["c1a", "c2b", "c3c"]),                # known AMR, scattered
    "NOVEL_mobile":   (rseq(40), ["c1b", "c2a", "c2d", "c3a", "c3d"]),  # novel, scattered
    "NOVEL_vertical": (rseq(40), ["c2a", "c2b", "c2c", "c2d"]),         # novel, = one clade
}

# unitig-caller .pyseer format:  SEQ | g1:1 g2:1 ...
with open(p("unitigs.pyseer"), "w") as fh:
    for _, (seq, carriers) in unitigs.items():
        fh.write(seq + " | " + " ".join(f"{c}:1" for c in carriers) + "\n")

# pyseer LMM association output (already filtered to "significant")
with open(p("significant.txt"), "w") as fh:
    fh.write("variant\taf\tfilter-pvalue\tlrt-pvalue\tbeta\tbeta-std-err\t"
             "variant_h2\tnotes\n")
    pvals = {"KNOWN_mobile": 2e-15, "NOVEL_mobile": 8e-12, "NOVEL_vertical": 3e-9}
    betas = {"KNOWN_mobile": 4.1, "NOVEL_mobile": 3.3, "NOVEL_vertical": 2.0}
    for name, (seq, carr) in unitigs.items():
        af = len(carr) / 12
        fh.write(f"{seq}\t{af:.3f}\t1e-6\t{pvals[name]:.0e}\t{betas[name]}\t"
                 f"0.5\t0.2\t\n")

# synthetic BLAST vs AMR db (as if blastn already ran)
# outfmt: 6 qseqid sseqid pident length qlen slen evalue bitscore stitle
with open(p("blast6.tsv"), "w") as fh:
    k = uid(unitigs["KNOWN_mobile"][0])
    fh.write(f"{k}\tblaNDM-1\t99.0\t40\t40\t813\t1e-18\t75\t"
             "blaNDM-1 carbapenem-hydrolysing class B beta-lactamase\n")
    # a weak, low-identity spurious hit for the novel-mobile unitig -> must be filtered
    nm = uid(unitigs["NOVEL_mobile"][0])
    fh.write(f"{nm}\tsome_gene\t70.0\t22\t40\t900\t0.3\t20\tunrelated hypothetical\n")

# patterns file (with duplicate patterns) for the Bonferroni threshold
with open(p("kmer_patterns.txt"), "w") as fh:
    fh.write("\n".join(["aabbaa", "aabbaa", "ababab", "ccddcc", "ababab"]) + "\n")

print("Wrote example inputs to", HERE)
for name, (seq, carr) in unitigs.items():
    print(f"  {name:14s} uid={uid(seq)}  carriers={len(carr)}")
