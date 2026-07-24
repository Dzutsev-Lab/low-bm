# Analysis Run Model

The processing workflow remains a Snakemake DAG. Analyses are intentionally
manual, explicit post-processing steps because different datasets need different
statistical and biological follow-up.

Use `./low-bm analysis` as the stable entry point instead of calling analysis
scripts directly. The underlying scripts remain available for transition and
debugging.

Standard analyses consume one resolved endpoint: `project.compiled_physeq` and,
when sequence lookup is needed, `project.compiled_asv_fasta`. That endpoint can
come directly from one processing batch or from `low-bm meta compile-phyloseq`.

## Configs

Initialize the default local analysis config:

```bash
./low-bm analysis init
```

This creates `config/local/analysis.yaml` from
`config/templates/analysis.yaml`. Edit the shared `project:` section once, then
edit or remove per-analysis sections as needed.

Focused configs are useful for one-off follow-up analyses. For example:

```bash
cp config/templates/analysis/blast-top-abundance.yaml config/local/analysis_blast_top_abundance.yaml
```

## Running Steps

Validate the config:

```bash
./low-bm analysis validate --analysis-config config/local/analysis.yaml
./low-bm analysis validate differential-abundance blast-confirmation --analysis-config config/local/analysis.yaml
```

Dry-run selected steps:

```bash
./low-bm analysis run abundance-barplots ordination \
  --analysis-config config/local/analysis.yaml \
  --dry-run
```

Run selected steps in order:

```bash
./low-bm analysis run abundance-barplots ordination \
  --analysis-config config/local/analysis.yaml
```

By default, the CLI runs scripts through the named conda/mamba environments
declared in `workflow/envs/R-tools-env.yaml` and `workflow/envs/bio-tools-env.yaml`.
If you already activated the right environment yourself, use:

```bash
./low-bm analysis run ordination --analysis-config config/local/analysis.yaml --env-mode direct
```

Each invocation writes `analysis_metadata.json`, `analysis_commands.txt`, and
step logs under `analysis_logs/`.

For multiple processing batches, compile first:

```bash
./low-bm meta compile-phyloseq --analysis-config config/local/meta.yaml
```

Then point `project.compiled_physeq` at `CompPhyseq.RData` and
`project.compiled_asv_fasta` at `MergedASV.fasta` in `config/local/analysis.yaml`.

## Example Workflows

Single-endpoint exploratory plots:

```bash
./low-bm analysis run abundance-barplots ordination \
  --analysis-config config/local/analysis.yaml
```

Multi-batch exploratory plots:

```bash
./low-bm meta compile-phyloseq --analysis-config config/local/meta.yaml
./low-bm analysis run abundance-barplots ordination --analysis-config config/local/analysis.yaml
```

Differential abundance followed by confirmatory BLAST:

```bash
./low-bm analysis run differential-abundance blast-confirmation \
  --analysis-config config/local/analysis.yaml
```

`blast-confirmation` expands to:

```text
blast-candidates -> blast-search -> blast-plots
```

Top-abundance BLAST without differential abundance:

```bash
./low-bm analysis run blast-confirmation \
  --analysis-config config/local/analysis_blast_top_abundance.yaml
```

Manual or reviewed taxa can use the same BLAST path by setting:

```yaml
blast_confirmation:
  candidate_source: "taxa_file"
  candidate_taxa_file: "Exp_Output/my_project/SelectedTaxa/reviewed_taxa.tsv"
  candidate_comparisons: ["ReviewedTaxa"]
```

The taxa file is a tab-separated table with these columns:

```text
source	comparison	taxa_level	taxon	selection_metric	selection_reason
manual	ReviewedTaxa	Genus	g__ExampleGenus		reviewed candidate
```
