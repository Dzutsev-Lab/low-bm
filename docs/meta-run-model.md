# Meta-Analysis Run Model

`low-bm meta` contains steps that combine or compare outputs from multiple
processing batches. These steps prepare analysis-ready endpoints or compare
batch-level result tables; they are separate from ordinary `low-bm analysis`
steps.

Processing outputs are pre-micRoclean, Kraken-annotated ASV endpoints.
micRoclean decontamination is an optional meta step so pre/post endpoints can be
compared deliberately.

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

## Optional micRoclean Decontamination

```bash
./low-bm meta decontaminate-phyloseq --analysis-config config/local/meta.yaml
```

This step reads one explicit phyloseq endpoint from
`meta_decontamination.input_physeq` and writes a decontaminated phyloseq plus
micRoclean reports under `meta_decontamination.output_dir`. It does not merge
or compile inputs on its own.

The default settings match the former processing-stage micRoclean behavior:
`SampleType == NegativeControl` marks controls, `ProcessingBatch` is used when
present, and the micRoclean research goal is `biomarker`. Keep controls in the
input object; filter controls later in downstream analyses if needed.

Running this on a compiled multi-trial endpoint is a sensitivity/QC choice. Use
it only when controls, sample-type labels, and processing-batch labels are
comparable across the compiled input.

For explicit HPC prefixes, provide the micRoclean environment separately:

```bash
./low-bm meta decontaminate-phyloseq \
  --analysis-config config/local/meta.yaml \
  --env-mode prefix \
  --manager /data/taylorng/conda/bin/mamba \
  --microclean-env-prefix /data/taylorng/conda/envs/low-bm-microclean
```

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
