#!/usr/bin/env python3
"""
triage_novelty.py — the bridge from "pyseer found significant unitigs" to
"here are candidate NOVEL resistance determinants, ranked, with a mobility read."

Inputs
  --significant  pyseer association output, filtered to significant hits
  --blast        BLAST of significant unitigs vs a known-AMR database, outfmt:
                 6 qseqid sseqid pident length qlen slen evalue bitscore stitle
  --mobility     mobility.tsv from mobility_score.py
Classification (per significant unitig)
  KNOWN_AMR        : passes identity+coverage to a gene in the AMR database
  NOVEL_CANDIDATE  : significant, survives the LMM, but matches no known AMR gene
Flags
  lineage_suspect  : NOVEL_CANDIDATE but acquired ~once (tracks one clade) -> the
                     signal may be residual population structure, not a mobile
                     causal element. De-prioritise / inspect.
A NOVEL_CANDIDATE that is BOTH significant AND high-mobility (many independent
acquisitions) is the highest-value target: an unexplained determinant that is
demonstrably jumping between hosts.

These are CANDIDATES, not calls. Unitigs tag variants in linkage, so follow up by
assembling the flagged unitigs into contigs, annotating (e.g. bakta/abricate),
checking synteny with mobile elements, and ideally phenotypic validation.

Usage:
    python triage_novelty.py --significant significant.txt --blast blast6.tsv \
        --mobility mobility.tsv --out-table candidates.tsv \
        --out-fasta novel_candidates.fasta --min-identity 80 --min-coverage 0.8
"""
import argparse
import hashlib
import math
import sys
import pandas as pd


def uid_of(seq):
    return hashlib.md5(seq.encode()).hexdigest()[:12]


BLAST_COLS = ["qseqid", "sseqid", "pident", "length", "qlen", "slen",
              "evalue", "bitscore", "stitle"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--significant", required=True)
    ap.add_argument("--blast", required=True)
    ap.add_argument("--mobility", required=True)
    ap.add_argument("--out-table", required=True)
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--min-identity", type=float, default=80.0)
    ap.add_argument("--min-coverage", type=float, default=0.8,
                    help="min fraction of the unitig covered by the AMR hit")
    args = ap.parse_args()

    # --- significant associations ---
    sig = pd.read_csv(args.significant, sep="\t")
    pcol = "lrt-pvalue" if "lrt-pvalue" in sig.columns else sig.columns[3]
    sig = sig.rename(columns={"variant": "unitig_seq", pcol: "pvalue"})
    sig["uid"] = sig["unitig_seq"].map(uid_of)
    keep = ["uid", "unitig_seq", "pvalue"]
    if "beta" in sig.columns:
        keep.append("beta")
    sig = sig[keep].drop_duplicates("uid")

    # --- BLAST vs known-AMR database ---
    try:
        bl = pd.read_csv(args.blast, sep="\t", header=None, names=BLAST_COLS)
    except (pd.errors.EmptyDataError, FileNotFoundError):
        bl = pd.DataFrame(columns=BLAST_COLS)
    if len(bl):
        bl["coverage"] = bl["length"] / bl["qlen"]
        good = bl[(bl["pident"] >= args.min_identity) &
                  (bl["coverage"] >= args.min_coverage)].copy()
        good = good.sort_values("bitscore", ascending=False).drop_duplicates("qseqid")
        amr = good.set_index("qseqid")[["sseqid", "pident", "coverage", "stitle"]]
    else:
        amr = pd.DataFrame(columns=["sseqid", "pident", "coverage", "stitle"])

    # --- mobility ---
    mob = pd.read_csv(args.mobility, sep="\t")[
        ["uid", "n_carriers", "independent_gains", "gains_per_carrier",
         "consistency_index"]]

    df = sig.merge(mob, on="uid", how="left")
    df = df.merge(amr, left_on="uid", right_index=True, how="left")

    df["classification"] = df["sseqid"].apply(
        lambda x: "KNOWN_AMR" if isinstance(x, str) else "NOVEL_CANDIDATE")

    # lineage-confounding suspicion: novel but acquired ~once
    df["lineage_suspect"] = ((df["classification"] == "NOVEL_CANDIDATE") &
                             (df["independent_gains"].fillna(0) <= 1))

    # heuristic priority for novel candidates:
    #   significance  +  mobility (log gains)  +  effect size
    def priority(r):
        if r["classification"] != "NOVEL_CANDIDATE" or r["lineage_suspect"]:
            return float("nan")
        sig_term = -math.log10(max(r["pvalue"], 1e-300))
        gains = r["independent_gains"] if pd.notna(r["independent_gains"]) else 1
        mob_term = math.log2(gains + 1)
        eff = abs(r["beta"]) if "beta" in r and pd.notna(r["beta"]) else 0.0
        return round(sig_term + 3 * mob_term + eff, 3)

    df["priority_score"] = df.apply(priority, axis=1)

    cols = ["uid", "classification", "lineage_suspect", "priority_score",
            "pvalue", "n_carriers", "independent_gains", "gains_per_carrier",
            "consistency_index", "sseqid", "pident", "coverage", "stitle"]
    if "beta" in df.columns:
        cols.insert(5, "beta")
    cols.append("unitig_seq")
    out = df[cols].sort_values(
        ["classification", "priority_score"], ascending=[True, False],
        na_position="last")
    out.to_csv(args.out_table, sep="\t", index=False)

    # FASTA of novel candidates (mobile ones first)
    nov = df[df["classification"] == "NOVEL_CANDIDATE"].sort_values(
        "priority_score", ascending=False, na_position="last")
    with open(args.out_fasta, "w") as fh:
        for _, r in nov.iterrows():
            tag = "lineage_suspect" if r["lineage_suspect"] else "mobile_candidate"
            fh.write(f">{r['uid']} {tag} gains={r['independent_gains']} "
                     f"p={r['pvalue']:.2g}\n{r['unitig_seq']}\n")

    n_known = (df["classification"] == "KNOWN_AMR").sum()
    n_nov = (df["classification"] == "NOVEL_CANDIDATE").sum()
    n_mob = int(((df["classification"] == "NOVEL_CANDIDATE") &
                 (~df["lineage_suspect"])).sum())
    sys.stderr.write(
        f"[triage] {len(df)} significant unitigs | known-AMR: {n_known} | "
        f"novel candidates: {n_nov} (mobile: {n_mob}, "
        f"lineage-suspect: {n_nov - n_mob})\n")


if __name__ == "__main__":
    main()
