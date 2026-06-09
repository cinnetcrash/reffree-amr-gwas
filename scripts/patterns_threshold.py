#!/usr/bin/env python3
"""
patterns_threshold.py — Bonferroni threshold for a unitig/k-mer GWAS.

Many unitigs share an identical presence/absence pattern across samples, so the
effective number of independent tests is the number of UNIQUE patterns, not the
number of unitigs. pyseer writes these with --output-patterns; this counts the
unique ones and returns the Bonferroni-corrected significance threshold.
(Equivalent to pyseer's bundled scripts/count_patterns.py.)

Usage:
    python patterns_threshold.py kmer_patterns.txt
    python patterns_threshold.py kmer_patterns.txt --alpha 0.05 --print-threshold
"""
import argparse
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("patterns")
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--print-threshold", action="store_true",
                    help="print ONLY the numeric threshold (for piping into awk)")
    args = ap.parse_args()

    patterns = set()
    with open(args.patterns) as fh:
        for line in fh:
            p = line.strip()
            if p:
                patterns.add(p)

    n = len(patterns)
    thr = args.alpha / n if n else float("nan")
    if args.print_threshold:
        print(f"{thr:.9g}")
    else:
        sys.stdout.write(f"Unique patterns: {n}\n")
        sys.stdout.write(f"Bonferroni threshold (alpha={args.alpha}): {thr:.3g}\n")


if __name__ == "__main__":
    main()
