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

By default, the CLI creates and reuses project-local conda/mamba environments
from the YAML files declared in `workflow/envs/R-tools-env.yaml` and
`workflow/envs/bio-tools-env.yaml`. These managed prefixes live under
`.low-bm/analysis/envs`, are named with an environment-file hash, and are run
with `mamba/conda run --prefix` so they do not depend on global conda env-name
lookup.

```bash
./low-bm analysis run abundance-barplots ordination \
  --analysis-config config/local/analysis.yaml
```

On HPC systems where you already maintain shared analysis envs and want to use
those instead of the project-managed copies, use explicit prefixes:

```bash
./low-bm analysis run abundance-barplots ordination \
  --analysis-config config/local/analysis.yaml \
  --env-mode prefix \
  --manager /data/taylorng/conda/bin/mamba \
  --r-env-prefix /data/taylorng/conda/envs/low-bm-r-tools
```

For BLAST confirmation, include the bio-tools prefix as well:

```bash
./low-bm analysis run blast-confirmation \
  --analysis-config config/local/analysis.yaml \
  --env-mode prefix \
  --r-env-prefix /data/taylorng/conda/envs/low-bm-r-tools \
  --bio-env-prefix /data/taylorng/conda/envs/low-bm-bio-tools
```

The same values can be supplied through `LOW_BM_R_TOOLS_PREFIX` and
`LOW_BM_BIO_TOOLS_PREFIX`. If you want the older global env-name behavior, use
`--env-mode named`. If you already activated the right environment yourself,
use:

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

## ANCOM-BC2 Differential Abundance

The config-driven differential-abundance step now uses ANCOM-BC2 exclusively.
Set `methods: ["ANCOMBC2"]` in new configs. Existing `ANCOMBC` values are
accepted as a migration alias, but new output directories and filenames use
`ANCOMBC2`.

ANCOM-BC2 differs from the old project wrapper in three practical ways:

- `fix_formula` replaces `formula`; `formula` is still accepted as an alias.
- The `group` variable may have two or more levels.
- Multi-level outputs can include primary/reference contrasts, global tests,
  all pairwise tests, Dunnett-type tests against a reference, and ordered trend
  tests.

Required comparison fields:

```yaml
differential_abundance:
  methods: ["ANCOMBC2"]
  comparisons:
    - name: "ControlVsTreatment"
      sample_filter:
        SampleType: ["Control", "Treatment"]
      factor_levels:
        SampleType: ["Control", "Treatment"]
      fix_formula: "SampleType"
      group: "SampleType"
      tests: ["primary"]
```

Use `factor_levels` to control model order. The first configured level is the
reference unless `reference_level` is set. Do not rely on alphabetical order for
scientific comparisons.

Three-level global and Dunnett-type example:

```yaml
differential_abundance:
  comparisons:
    - name: "DoseGlobalDunnett"
      sample_filter:
        DoseGroup: ["Vehicle", "Low", "High"]
      factor_levels:
        DoseGroup: ["Vehicle", "Low", "High"]
      fix_formula: "DoseGroup + ProcessingBatch"
      group: "DoseGroup"
      reference_level: "Vehicle"
      tests: ["global", "dunnet"]
```

Ordered trend example:

```yaml
differential_abundance:
  comparisons:
    - name: "DoseTrend"
      sample_filter:
        DoseGroup: ["Vehicle", "Low", "High"]
      ordered_levels: ["Vehicle", "Low", "High"]
      fix_formula: "DoseGroup"
      group: "DoseGroup"
      tests: ["trend"]
      trend_patterns: ["increasing", "decreasing"]
```

Trend tests require `ordered_levels`; the code intentionally refuses to infer a
trend order. Built-in trend patterns are `increasing` and `decreasing`. Advanced
users can provide `trend_control` directly, matching ANCOM-BC2's `contrast`,
`node`, and `B` controls.

Optional ANCOM-BC2 settings can be placed globally under
`differential_abundance` or overridden per comparison:

```yaml
differential_abundance:
  p_adj_method: "holm"
  pseudo_sens: true
  prv_cut: 0.10
  lib_cut: 1000
  s0_perc: 0.05
  struc_zero: true
  neg_lb: true
  n_cl: 1
  mdfdr_control:
    fwer_ctrl_method: "holm"
    B: 100
```

The main result file is:

```text
<output_dir>/ANCOMBC2/<comparison>/<trialID>_<comparison>_ANCOMBC2Results.tsv
```

It is a long table with `test` and `contrast` columns. Per-family files are also
written beside it: `ANCOMBC2Primary.tsv`, `ANCOMBC2Global.tsv`,
`ANCOMBC2Pairwise.tsv`, `ANCOMBC2Dunnett.tsv`, and `ANCOMBC2Trend.tsv` when the
corresponding output exists.

Important result fields:

- `log2FoldChange`: ANCOM-BC2 log fold changes converted to log2 scale for
  compatibility with existing plots and thresholds.
- `diff_abn`: ANCOM-BC2 differential-abundance call before sensitivity filtering.
- `passed_ss`: whether the taxon passed pseudo-count sensitivity analysis.
- `diff_robust`: sensitivity-robust ANCOM-BC2 call when available.
- `significance`: project-level call used by BLAST and plots. It uses
  `diff_robust` when present, otherwise `diff_abn`/`padj`, and still applies
  `lfc_cutoff` for contrast-level rows.
- `struc0`: structural-zero indicator. For legacy two-level contrasts this
  remains `group1`/`group2`; for multi-level comparisons it lists affected
  levels.

BLAST candidate selection can use all significant ANCOM-BC2 rows or a specific
slice:

```yaml
blast_confirmation:
  DA_method: "ANCOMBC2"
  DA_comparisons: ["DoseGlobalDunnett"]
  DA_result_test: "global"       # optional
  DA_result_contrast: null       # optional, useful for pairwise/dunnet
```

Meta differential abundance requires one contrast-level result with
`log2FoldChange` and `se`. Use `primary`, `pairwise`, or `dunnet`, and set
`DA_result_contrast` when a comparison has multiple contrasts. `global` and
`trend` summaries are rejected for meta-analysis because they do not represent a
single batch-level effect size.

Migration guide:

```text
Old field/output                 ANCOM-BC2 replacement
methods: ["ANCOMBC"]             methods: ["ANCOMBC2"]  # old value still accepted
formula                          fix_formula
coefficient                      inferred from group levels; no longer needed
structural_zero_groups           inferred from group levels; no longer needed
exactly two group levels          two or more group levels supported
ANCOMBC/<comparison>/...          ANCOMBC2/<comparison>/...
```

Patient duplicate handling for inferential analyses:

```yaml
differential_abundance:
  patient_duplicate_policy:
    action: "keep"     # keep, drop, error

survival_analysis:
  patient_duplicate_policy:
    action: "collapse" # collapse, drop, error
```

Differential-abundance duplicate handling is applied only to comparisons whose
formula includes `PatientID`; duplicate units are `PatientID + group`. Survival
duplicate units are `PatientID + sample_strata_col` when strata are enabled, or
`PatientID` when unstratified. Use `drop` to remove every sample in duplicated
units before modeling, or `error` to stop and review the metadata. Each run
writes a `PatientDuplicatePolicy.tsv` audit file beside the analysis outputs.

Config-driven XGBoost classification:

```bash
./low-bm analysis run xgboost --analysis-config config/local/analysis.yaml
```

Each binary classifier is defined under `xgboost_classification.models`. The
top-level XGBoost settings act as defaults, and per-model fields override them.
Set `split_group_col: null` on a model only when sample-level train/test
splitting is intended.

```yaml
xgboost_classification:
  batch_adj_covar: "Hospital"
  batch_adj_method: "ComBat"
  split_group_col: "PatientID"
  models:
    - name: "TumorVsNontumor"
      plot_title: "Tumor vs Nontumor"
      sample_filter:
        SampleType: ["Tumor", "Nontumor"]
      target:
        column: "SampleType"
        negative: ["Nontumor"]
        positive: ["Tumor"]
    - name: "HCCvsiCC"
      plot_title: "HCC vs iCC Tumors"
      sample_filter:
        SampleType: ["Tumor"]
      target:
        column: "TumorType"
        negative: ["iCC"]
        positive: ["HCC"]
```

The former `class_factors: [PatientSample, TumorType]` presets are no longer
recognized; encode those choices as explicit model entries like the examples
above.

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
