# clean workspace
rm(list = ls())
# clean garbadge
gc()
# clear graphics device
# dev.off()
# set decimals to digits instead of scientific
options(scipen = 999)
# set timeout limit (laxed)
options(timeout = 10000)
# ensure warnings are not converted into errors
options(warn = 1)
# set work directory - assumes script is run from polreason/ root
# If running from elsewhere, set working directory to polreason/ first
if (!dir.exists("analysis") || !dir.exists("generation")) {
  stop("Please run this script from the polreason/ root directory")
}

# # --- Source Configs ----------------------------------------------------- # #
YEAR  <- 2024

DIR_SCRIPTS <- 'analysis/scripts'
source(file = file.path(DIR_SCRIPTS,'0.config.R'))

# # --- Define Experiment Details ------------------------------------------ # #
rater_files <- list.files(
  path    = sprintf("generation/synthetic_data/year_%d", YEAR),
  pattern = "\\.csv$",
  full.names = FALSE
)

RATERrs_list <- c("gss", tools::file_path_sans_ext(rater_files))

for(r in seq_along(RATERrs_list)){

  RATER <- RATERrs_list[r]

  DIR_LLM <- sprintf('generation/synthetic_data/year_%d/%s', YEAR, RATER)
  
  DIR_OUT <- file.path( BASE_OUT_DIR, paste0(as.character(RATER), "-", as.character(YEAR)))
  if (!dir.exists(DIR_OUT)) dir.create(DIR_OUT, recursive = TRUE)
  
  source(file = file.path(DIR_SCRIPTS,'1.data_shaper.R'))
}

# I want to detach the different modules to make running if it stops at any point easier
# Imputation + Bootstrap loop

for(r in seq_along(RATERrs_list)){
  
  RATER <- RATERrs_list[r]
  
  if(RATER %in% c('arcee-ai_trinity-mini','minimax_minimax-m2') ){next}
  
  DIR_OUT <- file.path(BASE_OUT_DIR, paste0(as.character(RATER), "-", as.character(YEAR)))
  
  run_from_scratch <- FALSE
  
  source(file = file.path(DIR_SCRIPTS,'2.polychor_bootstrap.R'))
  
}

# Source viz utils
source(file = file.path(DIR_SCRIPTS,'v.common_utils.R'))

# Plotting draws from MVN loop
source(file = file.path(DIR_SCRIPTS,'v1_a.mvn_plot.R'))

# Saturn plot: Global constraint visualization via MAC (faceted)
source(file = file.path(DIR_SCRIPTS,'v1_b.saturn_plot.R'))

# Saturn animation: Cycling through quantiles Q95 → Q5 (optional, takes ~2-5 min)
# Uncomment to generate animated GIF:
# source(file = file.path(DIR_SCRIPTS,'v1_c.saturn_animation.R'))

# Plot constraint statistics
source(file = file.path(DIR_SCRIPTS,'v2_a.constraint_stats.R'))
source(file = file.path(DIR_SCRIPTS,'v2_b.constraint_stats_delta.R'))
source(file = file.path(DIR_SCRIPTS,'v3.missing_dimensions.R'))

