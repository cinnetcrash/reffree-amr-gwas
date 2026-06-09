# Reference-free AMR GWAS + mobilome triage

A Snakemake pipeline that finds **resistance determinants that are not in any
database**, and tells you which of them are **mobile** (jumping between hosts)
rather than clonally inherited.

```
assemblies + binary phenotype + core phylogeny
        │
        ├─ unitig-caller ───────────────► unitigs (pyseer format)
        ├─ core phylogeny ──► kinship K (Brownian covariance)
        │
        ▼
   pyseer --lmm  (mixed model, lineage-corrected)
        │
        ▼
   Bonferroni threshold (unique patterns) ─► significant unitigs
        │
        ├─ BLAST vs known-AMR DB  ──────────┐
        ├─ phylogenetic mobility (parsimony)┤
        ▼                                   ▼
                 triage  ─►  candidates.tsv + novel_candidates.fasta
```

## The idea

Curated tools (ResFinder, CARD, AMRFinderPlus) can only report what is already
in their databases, so a genuinely novel determinant is invisible until the
database catches up. Instead, let the **phenotype** tell you which sequence
features matter: a unitig GWAS (pyseer LMM) flags every unitig associated with
resistance, whether or not anyone has named it.

That alone is noisy, for two reasons this pipeline addresses directly:

1. **Population structure.** A unitig that merely marks a resistant lineage will
   look causal. The LMM (`--lmm` with a phylogeny-derived kinship) corrects for
   this, and the mobility score below provides a second, independent check.
2. **"Is it actually new?"** A significant unitig is only interesting as a *novel*
   determinant if it does **not** match a known AMR gene. The triage step BLASTs
   every significant unitig against a known-AMR database and splits them into
   `KNOWN_AMR` vs `NOVEL_CANDIDATE`.

Then it switches the unit of analysis from the organism to the **element**: for
each significant unitig it maps presence/absence onto the core phylogeny and uses
small-parsimony (Sankoff) to count how many *independent* times it was acquired.
One acquisition = clonal/vertical; many = mobile/horizontal.

### Mobility doubles as a confounding filter

This is the part worth internalising. A unitig that is GWAS-significant but was
acquired **once** and sits inside a single clade (`independent_gains = 1`, high
tree-consistency) is exactly what residual lineage confounding looks like — it
gets flagged `lineage_suspect` and de-prioritised. A unitig that is significant
**and** acquired independently across the tree is both more likely to be causal
*and* harder to explain as a lineage artifact. The highest-value hit is therefore
`NOVEL_CANDIDATE` + not `lineage_suspect` + high `independent_gains`: an
unexplained determinant demonstrably moving between hosts.

## Install

```bash
conda env create -f environment.yml
conda activate reffree-gwas
```

## Inputs (set in `config.yaml`)

- `assemblies_dir` — one assembled genome per file (`*.fa`).
- `phenotype` — TSV, `samples <tab> phenotype` with a header; phenotype is binary 0/1. Sample names must match assembly basenames and tree tip labels.
- `tree` — a **midpoint-rooted** core-genome phylogeny in Newick. Build it however you normally do (e.g. core alignment → IQ-TREE), and ideally remove recombination first (Gubbins/ClonalFrameML) so the kinship reflects vertical descent.
- `amr_db` — a known-AMR reference FASTA. Nucleotide (ResFinder/CARD → `amr_db_type: nucl`, blastn) or protein (AMRFinderPlus reference proteins → `amr_db_type: prot`, blastx).

## Run

```bash
snakemake --cores 8 --use-conda
```

A worked toy example (tiny synthetic inputs) lives in `example/` — regenerate it
with `python example/make_test_data.py` and point a config at those files, or just
inspect `example/candidates.tsv` to see the expected output shape.

## Outputs (`results/`)

- `associations.txt` — full pyseer LMM output. The estimated heritability h² is in `pyseer.log`; a near-zero h² means the phenotype has little additive genetic signal and you should be cautious.
- `significant.txt` — unitigs below the Bonferroni threshold.
- `candidates.tsv` — the answer table. Key columns:
  - `classification` — `KNOWN_AMR` or `NOVEL_CANDIDATE`.
  - `lineage_suspect` — `True` if a novel candidate was acquired ~once (possible residual confounding).
  - `priority_score` — heuristic rank for mobile novel candidates (significance + mobility + effect size); `NaN` for known or lineage-suspect hits.
  - `independent_gains`, `gains_per_carrier`, `consistency_index` — the mobility read.
  - `sseqid` / `stitle` / `pident` / `coverage` — the matched known-AMR gene, if any.
- `novel_candidates.fasta` — sequences of the novel candidates, mobile ones first, headers tagged `mobile_candidate` / `lineage_suspect`.

## Follow-up (these are candidates, not calls)

Unitigs tag variants *in linkage* with the causal change; a significant unitig
may be linked rather than itself causal, and a single unitig is short. So for a
promising `NOVEL_CANDIDATE`:

- assemble the flagged unitigs into longer contigs and annotate (bakta / abricate / RGI on the contig, not the unitig);
- check synteny — does it sit on/near a plasmid backbone or MGE (MOB-suite, integrons)? That corroborates the mobility signal;
- confirm the phenotype association holds out-of-sample, and ideally validate functionally.

## Notes

- The kinship script (`scripts/kinship_from_tree.py`) is a dependency-light
  reimplementation of the Brownian phylogenetic covariance. To use pyseer's own
  version instead: `python /path/to/pyseer/scripts/phylogeny_distance.py --lmm
  core.tree > results/phylogeny_K.tsv`.
- The fixed-effects route (`mash dist | square_mash` → `scree_plot_pyseer` →
  `--distances`) is a valid alternative to the LMM for very large sets; the LMM
  is recommended for unitigs and is what this pipeline uses.
- Tool command-line flags occasionally change between releases (unitig-caller and
  pyseer especially). If a rule errors, check `unitig-caller --help` /
  `pyseer --help` against the `shell:` lines in the `Snakefile`.
- The helper scripts (`scripts/*.py`) are plain Python (stdlib + numpy/pandas) and
  can be run standalone, outside Snakemake, on any equivalent inputs.

## Tests

The helper scripts are covered by a small pytest suite (kinship properties,
parsimony-based mobility on clonal vs scattered patterns, the Bonferroni
threshold, and the triage classification) that runs in CI on every push:

```bash
pip install numpy pandas pytest
pytest -q
```

Example fit: a foodborne-pathogen collection (e.g. *Salmonella* with
plasmid-borne resistance) is an ideal target — the LMM separates the clonal
background from the cargo, and the mobility score is exactly what distinguishes a
plasmid/integron-borne determinant from a lineage marker.
