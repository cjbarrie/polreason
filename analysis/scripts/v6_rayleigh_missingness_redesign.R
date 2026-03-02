################################################################################
# 99_rayleigh_missingness_redesign.R
#
# Redesigned Rayleigh missing-structure figure:
#   - fill centered at null (ratio = 1) on log2 scale
#   - fewer PCs: keep major PCs, pool tail into "Other PCs"
#   - rater ordering by PC1 over-realization (reveals main gradient)
#   - stability cue: alpha scaled by n_boot
#   - cleaner labels: PC labels on top, minimal bottom margin
#
# Outputs:
#   - missing_variance_rayleigh_<YEAR>_redesign.pdf      (variable-width)
#   - missing_variance_rayleigh_<YEAR>_redesign_equalwidth.pdf (tile heatmap)
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(jsonlite)
})

## ------------------------------------------------------------------
## Settings
## ------------------------------------------------------------------
YEAR           <- get0("YEAR", ifnotfound = 2024L)
RATER_GSS      <- "gss"
CONDITIONAL    <- TRUE
CUMVAR_TARGET  <- 0.90

# PC selection
KEEP_CUMVAR    <- 0.70   # keep PCs until 70% cumulative variance (within K)
KEEP_MIN_N     <- 4L     # always keep at least first 4 PCs
ADD_OTHER_BIN  <- TRUE
OTHER_LABEL    <- "Other PCs"

# Fill
LIMITS_LOG2    <- c(-1, 1)  # clamp to 0.5x .. 2x

# Ordering
ORDER_MODE     <- "pc1"  # "pc1" or "contrast"

# Stability
STABILITY_MODE <- "alpha" # "alpha" or "border"

# Output
OUTFILE_VW     <- file.path(BASE_VIZ_DIR, sprintf("missing_variance_rayleigh_%d_redesign.pdf", YEAR))
OUTFILE_EQ     <- file.path(BASE_VIZ_DIR, sprintf("missing_variance_rayleigh_%d_redesign_equalwidth.pdf", YEAR))
PDF_WIDTH      <- 22
PDF_HEIGHT     <- 10

## ------------------------------------------------------------------
## PC labels (interpretive, from JSON)
## ------------------------------------------------------------------
pc_labels <- tryCatch({
  pc_labels_tmp <- jsonlite::fromJSON(
    txt = file.path('analysis', 'GSS_PC_explain.json'),
    simplifyVector = TRUE
  )
  pc_labels_tmp$axis_label
}, error = function(e) NULL)

## ------------------------------------------------------------------
## Sanity checks
## ------------------------------------------------------------------
required_objects <- c("BASE_OUT_DIR", "BASE_VIZ_DIR", "persona_vars_global")
missing_objects  <- required_objects[!vapply(required_objects, exists, logical(1), inherits = TRUE)]
if (length(missing_objects)) {
  stop("Missing required objects: ", paste(missing_objects, collapse = ", "))
}

required_fns <- c("load_corr_for_rater", "available_raters",
                   "drop_persona_vars", "cond_belief_corr_on_persona")
missing_fns <- required_fns[!vapply(required_fns, exists, logical(1), mode = "function", inherits = TRUE)]
if (length(missing_fns)) {
  stop("Missing required functions: ", paste(missing_fns, collapse = ", "))
}

## ------------------------------------------------------------------
## Reuse helpers from v3.missing_dimensions.R
## (sourced via master.R, or copy them here if running standalone)
## ------------------------------------------------------------------
# These should already exist: .safe_varnames, .apply_persona,
# get_llm_union_belief_vars, get_overlapping_belief_vars,
# build_gss_basis, rayleigh_ratios_for_rater, summarize_ratios

for (fn in c(".safe_varnames", ".apply_persona", "get_llm_union_belief_vars",
             "get_overlapping_belief_vars", "build_gss_basis",
             "rayleigh_ratios_for_rater", "summarize_ratios")) {
  if (!exists(fn, mode = "function")) {
    stop("Helper function '", fn, "' not found. Source v3.missing_dimensions.R first.")
  }
}

## ------------------------------------------------------------------
## Step 1: Data pipeline (same as v3)
## ------------------------------------------------------------------
message("\n=== Rayleigh Redesign: Building data ===")

raters_all <- unique(c(available_raters(BASE_OUT_DIR, YEAR), RATER_GSS))
raters_llm <- setdiff(raters_all, RATER_GSS)

llm_union_vars <- get_llm_union_belief_vars(
  raters_llm   = raters_llm,
  base_out_dir = BASE_OUT_DIR,
  year         = YEAR,
  persona_vars = persona_vars_global,
  conditional  = CONDITIONAL
)

overlap_vars <- get_overlapping_belief_vars(
  raters         = raters_all,
  llm_union_vars = llm_union_vars,
  base_out_dir   = BASE_OUT_DIR,
  year           = YEAR,
  persona_vars   = persona_vars_global,
  conditional    = CONDITIONAL
)

gss_basis <- build_gss_basis(
  overlap_vars   = overlap_vars,
  base_out_dir   = BASE_OUT_DIR,
  year           = YEAR,
  persona_vars   = persona_vars_global,
  conditional    = CONDITIONAL,
  cumvar_target  = CUMVAR_TARGET
)

message("GSS basis: ", length(gss_basis$varnames), " vars, K = ", gss_basis$K,
        " PCs (", round(100 * gss_basis$cumvar[gss_basis$K], 1), "% var)")

# Compute Rayleigh ratios for all raters
summary_list <- list()
for (r in raters_all) {
  ratio_dt <- rayleigh_ratios_for_rater(
    rater = r, gss_basis = gss_basis, base_out_dir = BASE_OUT_DIR,
    year = YEAR, persona_vars = persona_vars_global,
    conditional = CONDITIONAL, K_use = gss_basis$K
  )
  if (!is.null(ratio_dt)) summary_list[[r]] <- summarize_ratios(ratio_dt)
}

summary_dt <- rbindlist(summary_list, use.names = TRUE, fill = TRUE)
if (!nrow(summary_dt)) stop("No summary results.")

message("Summary computed for ", uniqueN(summary_dt$rater), " raters")

## ------------------------------------------------------------------
## Step 2: PC selection + "Other PCs" bin
## ------------------------------------------------------------------
K <- gss_basis$K
lambdaK <- gss_basis$lambda[1:K]
shareK  <- lambdaK / sum(lambdaK)
cumsK   <- cumsum(shareK)

# Keep PCs until cumulative share reaches KEEP_CUMVAR, but at least KEEP_MIN_N
keep_idx <- which(cumsK <= KEEP_CUMVAR)
keep_idx <- unique(c(seq_len(min(KEEP_MIN_N, K)), keep_idx))
keep_idx <- sort(keep_idx)

other_idx <- setdiff(seq_len(K), keep_idx)

message("Keeping ", length(keep_idx), " PCs; pooling ", length(other_idx), " into 'Other'")

# PC label helper — always prefix with "PCk:" to guarantee uniqueness
get_pc_label <- function(k) {
  if (is.character(pc_labels) && length(pc_labels) >= k &&
      !is.na(pc_labels[[k]]) && nzchar(pc_labels[[k]])) {
    return(paste0("PC", k, ": ", pc_labels[[k]]))
  }
  paste0("PC", k)
}

# Build geometry table for kept PCs
pc_geom <- data.table(
  pc_index  = keep_idx,
  share_raw = shareK[keep_idx],
  pc_label  = vapply(keep_idx, function(k) get_pc_label(k), character(1))
)

# Add "Other PCs" row
if (ADD_OTHER_BIN && length(other_idx) > 0) {
  other_share <- sum(shareK[other_idx])
  pc_geom <- rbind(pc_geom, data.table(
    pc_index  = 999L,
    share_raw = other_share,
    pc_label  = OTHER_LABEL
  ))
}

# Normalize widths to sum to 1
pc_geom[, width := share_raw / sum(share_raw)]
pc_geom[, x_max := cumsum(width)]
pc_geom[, x_min := c(0, head(x_max, -1))]
pc_geom[, x_mid := (x_min + x_max) / 2]
pc_geom[, pct_label := sprintf("%.1f%%", 100 * share_raw)]

## ------------------------------------------------------------------
## Step 3: Aggregate "Other PCs" in summary_dt
## ------------------------------------------------------------------
# Kept PCs: pass through
summary_kept <- summary_dt[pc_index %in% keep_idx]

# Other PCs: share-weighted mean in log2 space
if (ADD_OTHER_BIN && length(other_idx) > 0) {
  summary_other <- summary_dt[pc_index %in% other_idx]

  if (nrow(summary_other) > 0) {
    # Merge in shares for weighting
    other_shares <- data.table(pc_index = other_idx, w = shareK[other_idx])
    summary_other <- merge(summary_other, other_shares, by = "pc_index")

    summary_other_agg <- summary_other[, .(
      pc_index   = 999L,
      ratio_med  = 2^weighted.mean(log2(pmax(ratio_med, 1e-6)), w = w, na.rm = TRUE),
      ratio_low  = 2^weighted.mean(log2(pmax(ratio_low, 1e-6)), w = w, na.rm = TRUE),
      ratio_high = 2^weighted.mean(log2(pmax(ratio_high, 1e-6)), w = w, na.rm = TRUE),
      n_boot     = as.integer(median(n_boot, na.rm = TRUE))
    ), by = rater]

    summary_plot <- rbind(summary_kept, summary_other_agg, use.names = TRUE, fill = TRUE)
  } else {
    summary_plot <- summary_kept
  }
} else {
  summary_plot <- summary_kept
}

## ------------------------------------------------------------------
## Step 4: Rater ordering
## ------------------------------------------------------------------
if (ORDER_MODE == "pc1") {
  # Sort by log2(ratio) for PC1 descending; GSS first
  rater_score <- summary_plot[pc_index == 1L, .(score = log2(ratio_med)), by = rater]
} else {
  # "contrast": log2(PC1) - mean(log2(PC2..PCk_kept))
  kept_not1 <- setdiff(keep_idx, 1L)
  score_pc1 <- summary_plot[pc_index == 1L, .(s1 = log2(ratio_med)), by = rater]
  score_rest <- summary_plot[pc_index %in% kept_not1,
                             .(sr = mean(log2(pmax(ratio_med, 1e-6)), na.rm = TRUE)), by = rater]
  rater_score <- merge(score_pc1, score_rest, by = "rater", all = TRUE)
  rater_score[, score := s1 - fifelse(is.na(sr), 0, sr)]
}

rater_order <- c(
  RATER_GSS,
  rater_score[rater != RATER_GSS][order(-score)]$rater
)
rater_order <- unique(rater_order)

## ------------------------------------------------------------------
## Step 5: Build the variable-width redesigned plot
## ------------------------------------------------------------------
message("\nBuilding redesigned plot...")

plot_dt <- merge(
  summary_plot[pc_index %in% pc_geom$pc_index],
  pc_geom,
  by = "pc_index",
  all.x = TRUE
)

plot_dt[, rater := factor(rater, levels = rater_order)]
plot_dt[, rater_num := as.integer(rater)]

# Fill on log2 scale (null = 0 = ratio of 1)
plot_dt[, fill_val := log2(pmax(ratio_med, 1e-6))]
plot_dt[, fill_val := pmin(pmax(fill_val, LIMITS_LOG2[1]), LIMITS_LOG2[2])]

# Stability cue: alpha scaled by n_boot
boot_range <- range(plot_dt$n_boot, na.rm = TRUE)
if (boot_range[1] == boot_range[2]) {
  plot_dt[, alpha_val := 1.0]
} else {
  plot_dt[, alpha_val := scales::rescale(n_boot, to = c(0.35, 1.0),
                                         from = boot_range)]
}

# Significance border: CI excludes null (log2 space: low > 0 or high < 0)
plot_dt[, sig := (log2(pmax(ratio_low, 1e-6)) > 0) | (log2(pmax(ratio_high, 1e-6)) < 0)]

plot_df <- as.data.frame(plot_dt)

max_y <- max(plot_df$rater_num, na.rm = TRUE)
min_y <- min(plot_df$rater_num, na.rm = TRUE)
pc_top_df   <- as.data.frame(unique(pc_geom[, .(x_mid, pct_label)]))
pc_label_df <- as.data.frame(pc_geom[, .(x_mid, pc_label)])

p_vw <- ggplot(plot_df, aes(
  xmin = x_min, xmax = x_max,
  ymin = rater_num - 0.45, ymax = rater_num + 0.45
)) +
  # Base rectangles with alpha stability cue
  geom_rect(aes(fill = fill_val, alpha = alpha_val)) +

  # Significance border overlay
  geom_rect(
    data = plot_df[plot_df$sig == TRUE, ],
    aes(fill = fill_val),
    color = "black", linewidth = 0.4, alpha = 1
  ) +

  # PC percentage labels at top
  geom_text(
    data = pc_top_df,
    aes(x = x_mid, y = max_y + 0.9, label = pct_label),
    inherit.aes = FALSE, size = 3.2, color = "gray40"
  ) +

  # PC name labels at top
  geom_text(
    data = pc_label_df,
    aes(x = x_mid, y = max_y + 1.4, label = pc_label),
    inherit.aes = FALSE, size = 3.8, fontface = "bold"
  ) +

  # Subtle PC separators
  geom_vline(
    xintercept = pc_geom$x_min[-1],
    color = "gray85", linewidth = 0.25
  ) +

  # Diverging fill anchored at 0 (= ratio of 1)
  scale_fill_gradient2(
    low      = "#2166AC",
    mid      = "white",
    high     = "#B2182B",
    midpoint = 0,
    limits   = LIMITS_LOG2,
    oob      = scales::squish,
    breaks   = c(-1, -0.5, 0, 0.5, 1),
    labels   = c("0.5\u00d7", "0.71\u00d7", "1\u00d7", "1.41\u00d7", "2\u00d7"),
    name     = expression(frac(Var(LLM), Var(GSS)) ~ "(log"[2] * " scale)")
  ) +

  scale_alpha_identity() +

  scale_y_continuous(
    breaks = seq_along(rater_order),
    labels = rater_order,
    name   = NULL,
    expand = c(0, 0),
    limits = c(min_y - 0.8, max_y + 2.0)
  ) +

  scale_x_continuous(
    breaks = NULL, name = NULL,
    expand = c(0, 0), limits = c(0, 1)
  ) +

  coord_cartesian(clip = "off") +

  theme_minimal(base_size = 11) +
  theme(
    axis.text.y    = element_text(size = 10),
    axis.text.x    = element_blank(),
    axis.ticks.x   = element_blank(),
    panel.grid     = element_blank(),
    plot.margin    = margin(t = 40, r = 10, b = 10, l = 10),
    axis.title.x   = element_blank(),
    legend.position = "right"
  )

## ------------------------------------------------------------------
## Step 6: Equal-width tile heatmap variant
## ------------------------------------------------------------------
message("Building equal-width tile variant...")

# Create discrete PC labels in order
pc_levels <- pc_geom[order(pc_index)]$pc_label

tile_dt <- copy(summary_plot[pc_index %in% pc_geom$pc_index])
tile_dt <- merge(tile_dt, pc_geom[, .(pc_index, pc_label)], by = "pc_index")
tile_dt[, pc_label := factor(pc_label, levels = pc_levels)]
tile_dt[, rater := factor(rater, levels = rater_order)]
tile_dt[, fill_val := log2(pmax(ratio_med, 1e-6))]
tile_dt[, fill_val := pmin(pmax(fill_val, LIMITS_LOG2[1]), LIMITS_LOG2[2])]
tile_dt[, sig := (log2(pmax(ratio_low, 1e-6)) > 0) | (log2(pmax(ratio_high, 1e-6)) < 0)]

if (boot_range[1] == boot_range[2]) {
  tile_dt[, alpha_val := 1.0]
} else {
  tile_dt[, alpha_val := scales::rescale(n_boot, to = c(0.35, 1.0), from = boot_range)]
}

tile_df <- as.data.frame(tile_dt)

p_eq <- ggplot(tile_df, aes(x = pc_label, y = rater)) +
  geom_tile(aes(fill = fill_val, alpha = alpha_val), color = "gray90", linewidth = 0.3) +

  # Significance border
  geom_tile(
    data = tile_df[tile_df$sig == TRUE, ],
    aes(fill = fill_val),
    color = "black", linewidth = 0.5, alpha = 1
  ) +

  scale_fill_gradient2(
    low      = "#2166AC",
    mid      = "white",
    high     = "#B2182B",
    midpoint = 0,
    limits   = LIMITS_LOG2,
    oob      = scales::squish,
    breaks   = c(-1, -0.5, 0, 0.5, 1),
    labels   = c("0.5\u00d7", "0.71\u00d7", "1\u00d7", "1.41\u00d7", "2\u00d7"),
    name     = expression(frac(Var(LLM), Var(GSS)) ~ "(log"[2] * " scale)")
  ) +

  scale_alpha_identity() +

  labs(x = NULL, y = NULL) +

  theme_minimal(base_size = 11) +
  theme(
    axis.text.x    = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
    axis.text.y    = element_text(size = 10),
    panel.grid     = element_blank(),
    plot.margin    = margin(t = 10, r = 10, b = 10, l = 10),
    legend.position = "right"
  )

## ------------------------------------------------------------------
## Step 7: Save
## ------------------------------------------------------------------
message("\nSaving variable-width plot: ", OUTFILE_VW)
ggsave(OUTFILE_VW, plot = p_vw, width = PDF_WIDTH, height = PDF_HEIGHT,
       device = "pdf", limitsize = FALSE)

message("Saving equal-width plot: ", OUTFILE_EQ)
ggsave(OUTFILE_EQ, plot = p_eq, width = 20, height = PDF_HEIGHT,
       device = "pdf", limitsize = FALSE)

message("\n=== Rayleigh Redesign Complete ===\n")
