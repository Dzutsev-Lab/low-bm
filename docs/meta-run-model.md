# Meta-Analysis Run Model

`low-bm meta` contains steps that combine or compare outputs from multiple
processing batches. These steps prepare analysis-ready endpoints or compare
batch-level result tables; they are separate from ordinary `low-bm analysis`
steps.

## Initialize Config

```bash
./low-bm meta init
```

This creates `config/local/meta.yaml` from `config/templates/meta.yaml`.

## Compile Phyloseq And ASV FASTA

```bash
./low-bm meta compile-phyloseq --analysis-config config/local/meta.yaml
```

Like `low-bm analysis`, meta steps use project-managed conda/mamba environments
from `workflow/envs/*.yaml` by default. If you already maintain a shared R
analysis env on an HPC system, run meta steps by explicit prefix:

```bash
./low-bm meta compile-phyloseq \
  --analysis-config config/local/meta.yaml \
  --env-mode prefix \
  --manager /data/taylorng/conda/bin/mamba \
  --r-env-prefix /data/taylorng/conda/envs/low-bm-r-tools
```

The compiler reads `meta_compile.batch_table`, `meta_compile.trial_list`, or
`meta_compile.physeqs`, then writes:

```text
<project.output_dir>/CompPhyseq.RData
<project.output_dir>/MergedASV.fasta
```

Sample metadata in the compiled phyloseq includes:

```text
SourceTrialName
SourceTrialID
SourcePhyseqPath
SourceOrder
```

Use those columns downstream as model covariates, plotting groupings, or audit
fields when a compiled endpoint spans multiple processing batches.

## Differential-Abundance Meta-Analysis

```bash
./low-bm meta differential-abundance --analysis-config config/local/meta.yaml
```

This wraps the existing two-result-table meta-analysis. It remains separate
from standard `low-bm analysis run differential-abundance`, which runs DA on one
phyloseq endpoint.

## Standard Analysis Handoff

After compilation, update `config/local/analysis.yaml`:

```yaml
project:
  compiled_physeq: "Exp_Output/<compiled_trial_name>/CompPhyseq.RData"
  compiled_asv_fasta: "Exp_Output/<compiled_trial_name>/MergedASV.fasta"
```

Then run ordinary analyses:

```bash
./low-bm analysis run abundance-barplots ordination differential-abundance \
  --analysis-config config/local/analysis.yaml
```
