# Replication Release Plan

## Strategy

This repository uses a **full-artifact, single-repository** release strategy. All code, synthetic data, pre-computed bootstrap outputs, and visualizations are tracked in git and pushed to GitHub as a single self-contained package.

This means:
- No external data archive (Zenodo/OSF) is needed for replication
- `git clone` gives users everything required to reproduce the analysis
- The repository is large (~1.2 GB) but all individual files are within GitHub's 100 MB per-file limit

## What is included

| Directory | Contents | Size |
|-----------|----------|------|
| `generation/scripts/` | Python + R scripts for GSS extraction, persona generation, and LLM querying | < 1 MB |
| `generation/data/` | Processed GSS extract, persona CSV (gitignored: raw `.dta` file) | ~1 MB |
| `generation/synthetic_data/year_2024/` | 30 LLM response CSVs (~7–8 MB each) | ~217 MB |
| `analysis/scripts/` | R analysis and visualization scripts | < 1 MB |
| `analysis/output/` | Per-model bootstrap `.rds` files (31 models × 5 files each) | ~840 MB |
| `analysis/viz/` | Publication-ready PDFs and animation GIF | ~191 MB |
| `analysis/output/*.csv` | Lightweight summary tables | < 5 MB |

## What is not included

- `generation/data/gss7224_r1.dta` — the 566 MB GSS cumulative data file is gitignored per NORC redistribution terms. Users must download it directly from https://gss.norc.org/get-the-data/stata.html if they want to re-run the data extraction or persona generation steps. It is **not needed** to reproduce the statistical analysis.

## Replication paths

**Quickest (analysis only):**
```bash
git clone <repo-url>
cd polreason
Rscript analysis/scripts/master.R
```
Uses pre-committed synthetic responses and pre-computed outputs.

**Full regeneration (requires GSS file + OpenRouter API key):**
```bash
# Extract GSS and generate personas
cd generation/scripts
Rscript 00a_create_gss_extract_multiyear.R --year 2024
Rscript 00b_generate_personas.R

# Query LLMs
export OPENROUTER_API_KEY="your-key"
python 01_generate_synthetic_GSS.py --year 2024 --all-models --personas 1000

# Run analysis
cd ../..
Rscript analysis/scripts/master.R
```
