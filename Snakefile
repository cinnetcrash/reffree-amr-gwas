# =====================================================================
# Reference-free AMR GWAS + mobilome triage
# unitig-caller -> pyseer (LMM) -> significance -> BLAST vs known-AMR
#               -> novelty triage + phylogenetic mobility score
# =====================================================================
# Run:  snakemake --cores 8 --use-conda
# Edit config.yaml for your paths.

configfile: "config.yaml"

OUT = config["outdir"]
BLAST_PROG = "blastn" if config["amr_db_type"] == "nucl" else "blastx"


rule all:
    input:
        f"{OUT}/candidates.tsv",
        f"{OUT}/novel_candidates.fasta",


# --- absolute paths to assemblies, one per line (unitig-caller input) ---
rule refs_list:
    output:
        f"{OUT}/input.txt",
    params:
        d=config["assemblies_dir"],
        ext=config.get("assembly_ext", "fa"),
    shell:
        "ls -1 {params.d}/*.{params.ext} | xargs -n1 readlink -f > {output}"


# --- build the population graph and call unitigs in pyseer format ---
rule unitig_call:
    input:
        f"{OUT}/input.txt",
    output:
        f"{OUT}/unitigs.pyseer",
    threads: config["threads"]
    params:
        pref=f"{OUT}/unitigs",
    shell:
        "unitig-caller --call --refs {input} --out {params.pref} --threads {threads}"


# --- pyseer reads the variant file gzip-only; keep the plain copy for mobility ---
rule gzip_unitigs:
    input:
        f"{OUT}/unitigs.pyseer",
    output:
        f"{OUT}/unitigs.pyseer.gz",
    shell:
        "gzip -kf {input}"


# --- LMM kinship from the core phylogeny (Brownian covariance) ---
# Faithful alternative: pyseer's own scripts/phylogeny_distance.py --lmm tree.
rule kinship:
    input:
        config["tree"],
    output:
        f"{OUT}/phylogeny_K.tsv",
    shell:
        "python scripts/kinship_from_tree.py {input} > {output}"


# --- the GWAS: linear mixed model, unitigs as the variant input ---
rule pyseer_lmm:
    input:
        pheno=config["phenotype"],
        unitigs=f"{OUT}/unitigs.pyseer.gz",
        K=f"{OUT}/phylogeny_K.tsv",
    output:
        assoc=f"{OUT}/associations.txt",
        patterns=f"{OUT}/patterns.txt",
    threads: config["threads"]
    log:
        f"{OUT}/pyseer.log",  # heritability h^2 is printed here
    shell:
        "pyseer --lmm --phenotypes {input.pheno} --kmers {input.unitigs} "
        "--similarity {input.K} --output-patterns {output.patterns} "
        "--cpu {threads} > {output.assoc} 2> {log}"


# --- Bonferroni threshold from the number of unique presence/absence patterns ---
rule threshold:
    input:
        f"{OUT}/patterns.txt",
    output:
        f"{OUT}/threshold.txt",
    shell:
        "python scripts/patterns_threshold.py {input} --alpha {config[alpha]} "
        "--print-threshold > {output}"


# --- keep unitigs below the threshold (lrt-pvalue is column 4 in LMM output) ---
rule filter_significant:
    input:
        assoc=f"{OUT}/associations.txt",
        thr=f"{OUT}/threshold.txt",
    output:
        f"{OUT}/significant.txt",
    shell:
        r"""thr=$(cat {input.thr}); """
        r"""awk -v t="$thr" 'NR==1 || ($4!="" && $4+0 < t)' {input.assoc} > {output}"""


rule sig_fasta:
    input:
        f"{OUT}/significant.txt",
    output:
        f"{OUT}/sig_unitigs.fasta",
    shell:
        "python scripts/unitigs_to_fasta.py {input} > {output}"


# --- known-AMR database for the "is this already known?" test ---
rule makeblastdb:
    input:
        config["amr_db"],
    output:
        touch(f"{OUT}/.amrdb.ready"),
    params:
        t=config["amr_db_type"],
    shell:
        "makeblastdb -in {input} -dbtype {params.t}"


rule blast_vs_amr:
    input:
        fa=f"{OUT}/sig_unitigs.fasta",
        db=config["amr_db"],
        ready=f"{OUT}/.amrdb.ready",
    output:
        f"{OUT}/blast6.tsv",
    threads: config["threads"]
    params:
        prog=BLAST_PROG,
    shell:
        '{params.prog} -query {input.fa} -db {input.db} -num_threads {threads} '
        '-max_target_seqs 5 -evalue 1e-6 '
        '-outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore stitle" '
        '> {output} || true'


# --- phylogenetic mobility (independent acquisitions) per significant unitig ---
rule mobility:
    input:
        tree=config["tree"],
        unitigs=f"{OUT}/unitigs.pyseer",
        sig=f"{OUT}/significant.txt",
    output:
        f"{OUT}/mobility.tsv",
    shell:
        "python scripts/mobility_score.py --tree {input.tree} "
        "--pyseer {input.unitigs} --significant {input.sig} --out {output}"


# --- join everything: classify known vs novel, flag lineage artifacts, rank ---
rule triage:
    input:
        sig=f"{OUT}/significant.txt",
        blast=f"{OUT}/blast6.tsv",
        mob=f"{OUT}/mobility.tsv",
    output:
        table=f"{OUT}/candidates.tsv",
        fasta=f"{OUT}/novel_candidates.fasta",
    params:
        mi=config["min_identity"],
        mc=config["min_coverage"],
    shell:
        "python scripts/triage_novelty.py --significant {input.sig} "
        "--blast {input.blast} --mobility {input.mob} "
        "--out-table {output.table} --out-fasta {output.fasta} "
        "--min-identity {params.mi} --min-coverage {params.mc}"
