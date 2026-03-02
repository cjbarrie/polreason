# Political Reasoning and Constraint in Large Language Models

This repository contains code and data for analyzing political constraint patterns in Large Language Model (LLM) responses to survey questions, based on synthetic personas derived from General Social Survey (GSS) data.

## Replication

This repository includes all code, synthetic data, pre-computed bootstrap outputs, and visualizations needed to reproduce the analysis in full. The quickest path is:

```bash
git clone <repo-url>
cd polreason
Rscript analysis/scripts/master.R
```

See `analysis/README.md`, `generation/README.md`, and `docs/replication_release_plan.md` for full details.

## Project Overview

This project investigates how different LLMs exhibit patterns of political constraint when answering survey questions as synthetic personas. Building on della Posta's (2020) framework for measuring belief constraint, we:

1. **Generate synthetic personas** from real GSS respondent demographics
2. **Query multiple LLMs** with comprehensive attitudinal survey questions
3. **Analyze constraint patterns** using polychoric correlations, principal components analysis, and bootstrap methods
4. **Compare LLM behavior** to human survey responses across 28+ different models

The analysis covers 52 survey questions spanning culture-war issues (abortion, immigration, civil liberties) and non-culture-war domains (institutional confidence, economic outlook, social trust).

## Repository Structure

```
polreason/
├── generation/              # Data generation pipeline
│   ├── scripts/            # Python and R scripts for creating synthetic data
│   │   ├── 00a_create_gss_extract_multiyear.R  # Extract GSS data by year
│   │   ├── 00b_generate_personas.R             # Create natural language personas
│   │   ├── 01_generate_synthetic_GSS.py        # Query LLMs (OpenRouter)
│   │   ├── 02_generate_synthetic_GSS_gpt5safe.py  # Alternative querying script
│   │   └── test_models.py                      # Test model availability
│   ├── data/               # Small data files (personas, extracts)
│   │   ├── gss2024_personas.csv                # Natural language personas (823KB)
│   │   ├── gss2024_dellaposta_extract.rds      # Processed GSS data (164KB)
│   │   └── gss2024_variable_summary.csv        # Variable coverage summary
│   └── synthetic_data/     # LLM-generated survey responses
│       └── year_2024/      # 30 CSV files, one per model
│
├── analysis/               # Statistical analysis pipeline
│   ├── README.md           # Analysis-specific guide
│   ├── scripts/            # R analysis scripts
│   │   ├── master.R                            # Main orchestration script
│   │   ├── 0.config.R                          # Configuration and helpers
│   │   ├── 1.data_shaper.R                     # Data preparation
│   │   ├── 2.polychor_bootstrap.R              # Bootstrap analysis
│   │   ├── v.common_utils.R                    # Shared visualization utilities
│   │   ├── v1_a.mvn_plot.R                     # Multivariate normal scatter plots
│   │   ├── v1_b.saturn_plot.R                  # Saturn plots (faceted, Option C)
│   │   ├── v1_c.saturn_animation.R             # Saturn animation (Q95 → Q5 GIF)
│   │   └── v2*.R, v3*.R                        # Additional visualization scripts
│   ├── output/             # Model results (31 directories)
│   │   ├── gss-2024/                           # Human GSS baseline
│   │   ├── anthropic_claude-sonnet-4.5-2024/   # Example model results
│   │   └── ...                                 # Other model outputs
│   ├── viz/                # Visualizations (35+ directories/files)
│   │   ├── constraint_violins_2024/            # Constraint comparisons
│   │   ├── mvn_2024/                           # Multivariate normal plots
│   │   └── [model-name]/                       # Per-model visualizations
│   └── GSS_PC_explain.json # GSS question metadata
│
├── requirements.txt        # Python dependencies
├── r-requirements.txt      # R package list
├── docs/                   # Release planning and extra documentation
├── .gitignore             # Git ignore rules
└── README.md              # This file
```

## Quick Start

If you want to reproduce the analysis from the already-committed synthetic responses, the default path is:

```bash
Rscript analysis/scripts/master.R
```

Run that command from the repository root after installing the R packages listed below.

## Installation

### Prerequisites

- **Python 3.8+** for data generation
- **R 4.0+** for statistical analysis
- **GSS Cumulative Data File** (see Data Requirements below)

### Python Environment

```bash
# Create virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

### R Packages

```r
# Install all required packages
install.packages(c(
  "data.table", "lavaan", "miceRanger", "ggplot2",
  "grid", "gridExtra", "grDevices", "scales",
  "ggnewscale", "irr", "haven", "dplyr",
  "stringr", "optparse"
))

# Optional: for Saturn animation (v1_c.saturn_animation.R)
install.packages(c("gganimate", "gifski"))
```

## Data Requirements

### GSS Cumulative Data File

The GSS cumulative data file (`gss7224_r1.dta`, 565MB) is **included locally** in `generation/data/` but **excluded from GitHub** via `.gitignore` due to size constraints.

**If cloning from GitHub**, you'll need to download the file separately:

1. **Download**: Visit https://gss.norc.org/get-the-data/stata.html
2. **File needed**: GSS 1972-2024 Cumulative Data (Release 1) in Stata format
3. **Filename**: `gss7224_r1.dta`
4. **Location**: Place in `generation/data/gss7224_r1.dta`

**Note**: If you only want to run the analysis pipeline (not generate new synthetic data), you can skip this step since the synthetic responses are already included in `generation/synthetic_data/`.

The GSS cumulative `.dta` file should remain an external download in any public replication release unless its redistribution terms clearly permit bundling.

## Usage

### Option 1: Run Analysis Only (Recommended for Exploration)

If you just want to reproduce the analysis using existing synthetic data:

```bash
# Navigate to polreason root
cd polreason/

# Run the master analysis script
Rscript analysis/scripts/master.R
```

This will:
- Load existing synthetic responses from `generation/synthetic_data/`
- Perform polychoric correlation analysis with bootstrap
- Generate constraint statistics
- Create visualizations in `analysis/viz/`

### Option 2: Generate New Synthetic Data

To create new synthetic survey responses from LLMs:

#### Step 1: Create GSS Extract (requires GSS data file)

```bash
# Navigate to the generation scripts directory
cd generation/scripts

# Extract 2024 data
Rscript 00a_create_gss_extract_multiyear.R --year 2024

# Generate natural language personas
Rscript 00b_generate_personas.R

# Return to project root
cd ../..
```

#### Step 2: Query LLMs

```bash
# Navigate to the generation scripts directory (if not already there)
cd generation/scripts

# Set your OpenRouter API key
export OPENROUTER_API_KEY="your-key-here"

# Query all models with 1000 personas (expensive! ~$100-500 depending on models)
python 01_generate_synthetic_GSS.py --year 2024 --all-models --personas 1000

# Or query specific models
python 01_generate_synthetic_GSS.py --year 2024 --models "anthropic/claude-sonnet-4.5,openai/gpt-5" --personas 1000

# Return to project root
cd ../..
```

See `python 01_generate_synthetic_GSS.py --help` for all options.

#### Step 3: Run Analysis

```bash
# From polreason root directory
Rscript analysis/scripts/master.R
```

## Key Findings & Outputs

### Constraint Metrics

The analysis produces several key metrics of political constraint:

1. **PC1 Variance Explained**: How much variance the first principal component captures
2. **Effective Dependence (D_e)**: Average absolute polychoric correlation
3. **Missing Dimensions**: Analysis of variance unexplained by PC1
4. **Cohen's Kappa**: Agreement between different runs of the same persona

### Visualizations

- **Saturn plots** (`v1_b.saturn_plot.R`): Publication-ready faceted visualization comparing LLM constraint to human (GSS) baseline
  - Faceted layout: each quantile (Q25, Q50, Q75, Q90) shown in separate panel
  - Implements "Option C" highlighting: only models with significant constraint difference (Δ ≥ 0.10) vs GSS are colored
  - Non-highlighted LLMs shown as transparent gray "spaghetti" for context
  - GSS displayed as bold black contours for easy comparison
  - Reference circle (ρ=0) shows independence baseline
  - Parameters: `delta_min` (threshold), `top_n_total` (limit highlights)

- **Saturn animation** (`v1_c.saturn_animation.R`): Animated GIF cycling through quantile levels (Q95 → Q5)
  - Shows how constraint contours evolve from tightest (Q95) to weakest (Q5) correlations
  - 19 frames covering full quantile range in 5% increments
  - Top-N most constrained models highlighted throughout animation
  - Smooth transitions with cubic easing
  - Requires: `gganimate` and `gifski` packages
  - Optional in `master.R` (uncomment to generate, takes 2-5 minutes)

- **Constraint violins**: Compare constraint levels across models and education groups
- **Polychoric correlation matrices**: Triangle plots showing pairwise belief correlations
- **MVN scatter plots**: Multivariate normal draws from correlation matrices
- **Missing variance plots**: Rayleigh distribution analysis of unexplained variance

All visualizations are saved as PDFs in `analysis/viz/`.

## Models Included

The repository includes results for 28+ LLMs across various families:

- **OpenAI**: GPT-5, GPT-5-mini, GPT-4o-mini, GPT-OSS-120b
- **Anthropic**: Claude Sonnet 4.5, Claude 3.7 Sonnet
- **Google**: Gemini 2.5 Flash/Lite, Gemma 3 12B
- **DeepSeek**: v3, v3.1, v3.2
- **Meta**: Llama 3.1/3.3/4
- **Mistral**: Large, Medium, Small, Nemo
- **Others**: Qwen, Kimi, Grok, GLM, Cohere, AI21, Allen AI, and more

Plus **GSS-2024** human baseline for comparison.

## Citation

If you use this code or data, please cite:

> Barrie, C. & Cerina, R. (2026). Synthetic personas distort the structure of human belief systems. *SocArXiv*. https://doi.org/10.31235/osf.io/n7fq8_v1

A `CITATION.cff` file is included for software citation managers.

Related work:
- della Posta, D. (2020). "Pluralistic Collapse: The 'Oil Spill' Model of Mass Opinion Polarization." *American Sociological Review*, 85(3), 507-536.

## Acknowledgments

- General Social Survey (GSS) data from NORC at the University of Chicago
- LLM API access via OpenRouter
- Built on della Posta's constraint measurement framework

## Troubleshooting

### "Please run this script from the polreason/ root directory"

Make sure you're in the `polreason/` directory when running R scripts:
```bash
cd polreason/
Rscript analysis/scripts/master.R
```

### Missing GSS data file

Download `gss7224_r1.dta` from https://gss.norc.org/get-the-data/stata.html and place in `generation/data/`.

### API rate limits

The OpenRouter API has rate limits. The script includes retry logic with exponential backoff. For large runs, consider:
- Running overnight
- Using the `--max-workers` parameter to reduce concurrency
- Splitting across multiple days

### Memory issues in R

The bootstrap analysis can be memory-intensive. If you encounter issues:
- Close other applications
- Reduce `B` (bootstrap iterations) in `analysis/scripts/0.config.R`
- Run models individually instead of the full batch

## Development

To modify or extend this project:

1. **Add new models**: Edit `POPULAR_MODELS` in `generation/scripts/01_generate_synthetic_GSS.py`
2. **Add new questions**: Edit `GSS_QUESTIONS_*` dictionaries in the same file
3. **Modify analysis**: See `analysis/scripts/0.config.R` for parameters (bootstrap iterations, etc.)
4. **Add visualizations**: Create new scripts following the `v*.R` pattern

## Repository Metadata

- **Created**: December 2024
- **Data Year**: 2024 (GSS wave)
- **Models**: 28+ LLMs + human baseline
- **Survey Items**: 52 questions (30 culture-war, 22 non-culture-war)
- **Personas**: 1000 synthetic respondents per model
- **Local Repository Size During Audit**: ~2.7 GB
