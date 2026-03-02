# Analysis Pipeline

This directory contains the reproducible analysis pipeline for comparing political constraint in human GSS responses and LLM-generated survey responses.

## Entry Point

Run the full analysis from the repository root:

```bash
Rscript analysis/scripts/master.R
```

`master.R` assumes the current working directory is the project root. It reads synthetic responses from `generation/synthetic_data/year_2024/`, writes intermediate model outputs to `analysis/output/`, and writes figures to `analysis/viz/`.

## Script Map

- `scripts/0.config.R`: shared configuration, helper functions, bootstrap counts, and question metadata.
- `scripts/1.data_shaper.R`: harmonizes GSS and model response files into the shapes expected downstream.
- `scripts/2.polychor_bootstrap.R`: runs the multiple-imputation and bootstrap workflow.
- `scripts/v.common_utils.R`: helper functions shared across visualization scripts.
- `scripts/v1_a.mvn_plot.R`: multivariate-normal draws and pairwise visualization outputs.
- `scripts/v1_b.saturn_plot.R`: faceted Saturn plots comparing model constraint to the GSS baseline.
- `scripts/v1_c.saturn_animation.R`: optional animated Saturn plot workflow.
- `scripts/v2_a.constraint_stats.R`: summary constraint metrics.
- `scripts/v2_b.constraint_stats_delta.R`: delta-style comparisons against the GSS baseline.
- `scripts/v3.missing_dimensions.R`: missing-dimension and unexplained-variance analysis.
- `scripts/v4.kappa.R`: cross-run agreement diagnostics.
- `scripts/v5.combine_stability_plots.R`: combines stability figures.
- `scripts/v6_rayleigh_missingness_redesign.R`: redesigned missingness visualization.

## Inputs

- `../generation/synthetic_data/year_2024/*.csv`: one response file per model.
- `../generation/data/gss2024_dellaposta_extract.rds`: processed GSS extract used as the human baseline.
- `GSS_PC_explain.json`: metadata used for principal-component interpretation.

## Outputs

- `output/{model}-{year}/`: per-model `.rds` objects and harmonized data.
- `output/*.csv`: aggregate summary tables used for downstream reporting.
- `viz/`: PDFs and GIFs produced by the visualization scripts.

## Notes For Replication

- The default bootstrap settings live in `scripts/0.config.R` (`B = 500`, `B_MI = 30`) and materially affect run time.
- The repository currently includes precomputed outputs. A release-oriented remote copy should make a deliberate choice about whether to keep those artifacts in Git, move them to Git LFS, or host them in an archival store such as OSF or Zenodo.
- If you modify any script parameters, document the change in the top-level replication notes before publishing the result.
