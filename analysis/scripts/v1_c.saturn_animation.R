################################################################################
# Saturn Plot Animation — Cycling Through Quantiles (Q95 → Q5)
#
# Creates animated GIF showing how constraint contours evolve across the
# distribution of pairwise correlations, from tightest (Q95) to weakest (Q5).
#
# Animation shows:
# - Each frame represents one quantile level (Q95, Q90, Q85, ..., Q10, Q5)
# - Models with significant Δ vs GSS are highlighted in color
# - Gray spaghetti for non-highlighted models
# - Bold black GSS contour for comparison
# - Reference circle (ρ=0) baseline
#
# Requires: gganimate, gifski packages for animation
#
# Assumes:
# - BASE_OUT_DIR, BASE_VIZ_DIR, YEAR exist (if running main execution)
# - available_raters() and load_corr_for_rater() exist (from utilities)
################################################################################

suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' required")
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Package 'data.table' required")
  if (!requireNamespace("gganimate", quietly = TRUE)) stop("Package 'gganimate' required")
  if (!requireNamespace("gifski", quietly = TRUE)) stop("Package 'gifski' required for GIF rendering")

  library(ggplot2)
  library(data.table)
  library(gganimate)
})

# ---- Helper functions (reuse from v1_b) ----

compute_mac <- function(R) {
  if (is.null(R)) return(NA_real_)
  p <- ncol(R)
  if (p < 2) return(NA_real_)
  upper_tri <- R[upper.tri(R, diag = FALSE)]
  mean(abs(upper_tri), na.rm = TRUE)
}

compute_mac_bootstraps <- function(corr_list) {
  vapply(corr_list, compute_mac, numeric(1L))
}

extract_all_abs_cors <- function(R) {
  if (is.null(R)) return(numeric(0))
  p <- ncol(R)
  if (p < 2) return(numeric(0))
  abs(R[upper.tri(R, diag = FALSE)])
}

compute_correlation_quantiles <- function(corr_list, probs) {
  qmat <- matrix(NA_real_, nrow = length(corr_list), ncol = length(probs))
  colnames(qmat) <- paste0("q", probs * 100)

  for (i in seq_along(corr_list)) {
    v <- extract_all_abs_cors(corr_list[[i]])
    if (length(v) > 0) qmat[i, ] <- quantile(v, probs = probs, na.rm = TRUE)
  }

  data.frame(
    quantile   = paste0("Q", probs * 100),
    prob       = probs,
    median_rho = apply(qmat, 2, median,   na.rm = TRUE),
    q025       = apply(qmat, 2, quantile, probs = 0.025, na.rm = TRUE),
    q975       = apply(qmat, 2, quantile, probs = 0.975, na.rm = TRUE),
    row.names  = NULL
  )
}

ellipse_from_rho <- function(rho, prob = 0.5, n = 361L) {
  rho <- pmax(pmin(rho, 0.9999), -0.9999)
  Sigma <- matrix(c(1, rho, rho, 1), 2, 2)
  r <- sqrt(qchisq(prob, df = 2))
  theta <- seq(0, 2 * pi, length.out = n)
  U <- rbind(cos(theta), sin(theta))
  e <- eigen(Sigma, symmetric = TRUE)
  A <- e$vectors %*% diag(sqrt(pmax(e$values, 0))) %*% t(e$vectors)
  XY <- t(A %*% (r * U))
  colnames(XY) <- c("x", "y")
  XY
}

reference_circle <- function(prob = 0.5, n = 361L) {
  r0 <- sqrt(qchisq(prob, df = 2))
  theta <- seq(0, 2 * pi, length.out = n)
  data.frame(x = r0 * cos(theta), y = r0 * sin(theta))
}

# ---- Animation builder ----

create_saturn_animation <- function(
    base_out_dir      = BASE_OUT_DIR,
    base_viz_dir      = BASE_VIZ_DIR,
    year              = YEAR,
    prob              = 0.5,
    quantile_seq      = seq(0.95, 0.05, by = -0.01),  # Q95 → Q5

    highlight_gss     = TRUE,

    # Option C parameters
    delta_min         = 0.10,
    top_n_total       = 10L,  # Keep top-N models overall (across all quantiles)

    # Styling
    spaghetti_alpha   = 0.12,
    spaghetti_width   = 0.35,
    gss_width         = 1.2,
    highlight_width   = 0.8,

    # Animation settings
    fps               = 10,
    duration          = 10,  # seconds
    loop              = TRUE,

    output_file       = NULL
) {

  message("\n=== Saturn Animation (Q95 → Q5) ===")
  message("Quantile range: ", min(quantile_seq), " to ", max(quantile_seq))
  message("Frames: ", length(quantile_seq), " | FPS: ", fps, " | Duration: ~", duration, "s")

  # Get all raters
  raters <- available_raters(base_out_dir = base_out_dir, year = year)
  message("Found ", length(raters), " raters")

  if (length(raters) == 0) {
    stop("No raters found in ", base_out_dir)
  }

  # Load correlation matrices for all raters
  all_corr_data <- list()
  mac_summary <- data.table(rater = character(), median_mac = numeric())

  for (rater in raters) {
    message("Loading: ", rater)
    corr_list <- tryCatch(
      load_corr_for_rater(rater = rater, base_out_dir = base_out_dir, year = year, strict = FALSE),
      error = function(e) NULL
    )

    if (!is.null(corr_list) && length(corr_list) > 0) {
      all_corr_data[[rater]] <- corr_list
      mac_vec <- compute_mac_bootstraps(corr_list)
      mac_summary <- rbindlist(list(
        mac_summary,
        data.table(rater = rater, median_mac = median(mac_vec, na.rm = TRUE))
      ))
    }
  }

  if (length(all_corr_data) == 0) {
    stop("No correlation data loaded successfully")
  }

  # Compute quantiles for ALL quantile levels in sequence
  message("\nComputing quantiles for ", length(quantile_seq), " levels...")

  all_quantile_data <- list()

  for (rater in names(all_corr_data)) {
    quant_stats <- compute_correlation_quantiles(all_corr_data[[rater]], probs = quantile_seq)
    quant_stats$rater <- rater
    all_quantile_data[[rater]] <- quant_stats
  }

  quantile_dt <- rbindlist(all_quantile_data)

  # Identify GSS
  quantile_dt[, is_gss := tolower(rater) == "gss"]

  # For each quantile, compute delta vs GSS and identify highlights
  gss_rho <- quantile_dt[is_gss == TRUE, .(quantile, gss_rho = median_rho)]

  model_rho <- quantile_dt[is_gss == FALSE][gss_rho, on = "quantile"]
  model_rho[, delta := median_rho - gss_rho]

  # Identify top-N models by maximum delta across all quantiles
  top_models <- model_rho[, .(max_delta = max(delta, na.rm = TRUE)), by = rater]
  top_models <- top_models[order(-max_delta)]
  top_models <- head(top_models, top_n_total)

  message("\nTop ", top_n_total, " models by max delta:")
  print(top_models)

  # Mark highlights: top-N models, and only where delta >= delta_min
  quantile_dt[, is_highlight := FALSE]
  if (nrow(top_models) > 0) {
    for (i in 1:nrow(top_models)) {
      m <- top_models$rater[i]
      quantile_dt[rater == m & !is_gss, is_highlight := TRUE]
    }

    # Further filter: only highlight where delta >= threshold
    model_rho_highlight <- model_rho[rater %in% top_models$rater]
    weak_pairs <- model_rho_highlight[delta < delta_min, .(rater, quantile)]
    if (nrow(weak_pairs) > 0) {
      quantile_dt[weak_pairs, is_highlight := FALSE, on = .(rater, quantile)]
    }
  }

  # Generate ellipse data for all (rater, quantile) pairs
  message("\nGenerating ellipse paths...")

  ellipse_list <- list()

  for (i in 1:nrow(quantile_dt)) {
    rater_i <- quantile_dt$rater[i]
    quant_i <- quantile_dt$quantile[i]
    rho_i   <- quantile_dt$median_rho[i]
    is_gss_i <- quantile_dt$is_gss[i]
    is_hi_i  <- quantile_dt$is_highlight[i]

    xy <- ellipse_from_rho(rho = rho_i, prob = prob, n = 361L)

    ellipse_list[[paste(rater_i, quant_i, sep = "_")]] <- data.frame(
      x          = xy[, "x"],
      y          = xy[, "y"],
      rater      = rater_i,
      quantile   = quant_i,
      rho        = rho_i,
      is_gss     = is_gss_i,
      is_highlight = is_hi_i,
      stringsAsFactors = FALSE
    )
  }

  ellipse_dt <- rbindlist(ellipse_list)
  ellipse_dt[, group_id := paste(rater, quantile, sep = "_")]

  # Order raters by MAC
  mac_summary <- mac_summary[order(-median_mac)]
  rater_order <- mac_summary$rater
  ellipse_dt[, rater := factor(rater, levels = rater_order)]

  # Order quantiles for animation (Q95 → Q5)
  quantile_levels <- paste0("Q", quantile_seq * 100)
  ellipse_dt[, quantile := factor(quantile, levels = quantile_levels)]

  # Create layer subsets
  gss_dt  <- ellipse_dt[is_gss == TRUE]
  hi_dt   <- ellipse_dt[is_gss == FALSE & is_highlight == TRUE]
  rest_dt <- ellipse_dt[is_gss == FALSE & is_highlight == FALSE]

  message("Total ellipses: ", length(unique(ellipse_dt$group_id)))
  message("GSS ellipses: ", length(unique(gss_dt$group_id)))
  message("Highlighted ellipses: ", length(unique(hi_dt$group_id)))
  message("Background ellipses: ", length(unique(rest_dt$group_id)))

  # Reference circle
  ref <- reference_circle(prob = prob, n = 361L)

  # Color palette for highlighted models
  hi_models <- unique(as.character(hi_dt$rater))
  hi_models <- hi_models[order(match(hi_models, rater_order))]

  if (length(hi_models) > 0) {
    hi_palette <- setNames(
      grDevices::hcl.colors(length(hi_models), palette = "Dark 3"),
      hi_models
    )
  } else {
    hi_palette <- NULL
  }

  # Plot limits
  r0 <- sqrt(qchisq(prob, df = 2))
  lim <- 1.25 * r0

  # Create animated plot
  message("\nCreating animation...")

  p <- ggplot() +
    geom_hline(yintercept = 0, color = "gray90", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "gray90", linewidth = 0.3) +
    geom_path(
      data = ref, aes(x, y),
      color = "gray80", linewidth = 0.4, linetype = "dashed", lineend = "round"
    ) +

    # Background spaghetti
    geom_path(
      data = rest_dt,
      aes(x, y, group = group_id),
      color = "gray50",
      alpha = spaghetti_alpha,
      linewidth = spaghetti_width,
      lineend = "round"
    ) +

    # Highlighted models
    {if (nrow(hi_dt) > 0) geom_path(
      data = hi_dt,
      aes(x, y, group = group_id, color = rater),
      alpha = 0.65,
      linewidth = highlight_width,
      lineend = "round"
    )} +

    # GSS on top
    {if (nrow(gss_dt) > 0) geom_path(
      data = gss_dt,
      aes(x, y, group = group_id),
      color = if (highlight_gss) "black" else "gray20",
      alpha = if (highlight_gss) 1.0 else 0.4,
      linewidth = if (highlight_gss) gss_width else 0.6,
      lineend = "round"
    )} +

    coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim), expand = FALSE) +
    labs(
      title = "",
      x = "Standardized Dimension 1",
      y = "Standardized Dimension 2"
    ) +
    theme_classic(base_size = 16) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(face = "bold"),
      panel.grid = element_blank(),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.margin = margin(10, 10, 10, 10)
    )

  # Add color scale if highlights exist
  if (!is.null(hi_palette)) {
    p <- p + scale_color_manual(
      values = hi_palette,
      name = "Highlighted models"
    ) +
      guides(
        color = guide_legend(
          ncol = 3,
          override.aes = list(alpha = 1, linewidth = 1.5)
        )
      )
  }

  # Add animation transitions
  p_anim <- p +
    transition_states(
      quantile,
      transition_length = 2,
      state_length = 1
    ) +
    ease_aes('cubic-in-out')

  # Render animation
  if (is.null(output_file)) {
    mvn_dir <- file.path(base_viz_dir, sprintf("mvn_%d", year))
    dir.create(mvn_dir, recursive = TRUE, showWarnings = FALSE)
    output_file <- file.path(mvn_dir, "saturn_animation.gif")
  }

  message("\nRendering GIF: ", output_file)
  message("This may take a few minutes...")

  nframes <- length(quantile_seq)

  anim <- animate(
    p_anim,
    nframes = nframes,
    fps = fps,
    duration = duration,
    width = 1100,
    height = 1100,
    renderer = gifski_renderer(output_file, loop = loop)
  )

  message("\n=== Animation Complete ===")
  message("Saved to: ", output_file)
  message("Frames: ", nframes, " | FPS: ", fps, " | Loop: ", loop)

  invisible(anim)
}


################################################################################
# MAIN EXECUTION
################################################################################

if (exists("BASE_OUT_DIR") && exists("BASE_VIZ_DIR") && exists("YEAR")) {
  message("\n====================================================================")
  message("Saturn Animation: Cycling Through Quantiles (Q95 → Q5)")
  message("====================================================================\n")

  saturn_anim <- create_saturn_animation(
    base_out_dir   = BASE_OUT_DIR,
    base_viz_dir   = BASE_VIZ_DIR,
    year           = YEAR,
    prob           = 0.5,
    quantile_seq   = seq(0.95, 0.05, by = -0.01),  # 19 frames

    delta_min      = 0.10,
    top_n_total    = 10L,  # Keep top 10 most constrained models

    fps            = 10,
    duration       = 10,
    loop           = TRUE
  )

} else {
  message("BASE_OUT_DIR, BASE_VIZ_DIR, and YEAR must be defined to run.")
}
