#!/usr/bin/env python3
"""
unitigs_to_fasta.py — write significant unitigs to FASTA with stable uids.

The uid is md5(sequence)[:12], identical to mobility_score.py, so BLAST output
(qseqid = uid) joins cleanly with the mobility table and the association table.

Usage:
    python unitigs_to_fasta.py significant.txt > sig_unitigs.fasta
"""
import hashlib
import sys


def uid_of(seq):
    return hashlib.md5(seq.encode()).hexdigest()[:12]


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: unitigs_to_fasta.py significant.txt > out.fasta")
    with open(sys.argv[1]) as fh:
        fh.readline()  # header
        seen = set()
        for line in fh:
            if not line.strip():
                continue
            seq = line.split("\t")[0].strip()
            uid = uid_of(seq)
            if uid in seen:
                continue
            seen.add(uid)
            sys.stdout.write(f">{uid}\n{seq}\n")


if __name__ == "__main__":
    main()
