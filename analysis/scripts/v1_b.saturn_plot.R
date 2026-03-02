################################################################################
# Saturn Plots (Faceted) — Option C Highlighting
# Highlight models that are meaningfully more constrained than GSS at "relevant"
# quantiles, using Δ = rho_model(quantile) - rho_GSS(quantile).
#
# FACETED LAYOUT: Each quantile (Q25, Q50, Q75, Q90) shown in separate panel
# - Horizontal layout (nrow = 1) for easy left-to-right comparison
# - Each panel shows constraint at that specific quantile level
# - All non-highlighted non-GSS models: thin grey spaghetti
# - Highlighted model/quantile paths: colored (small legend)
# - GSS: bold black on top (optional)
# - Reference circle: rho = 0 at the same probability mass (prob)
#
# Assumes:
# - BASE_OUT_DIR, BASE_VIZ_DIR, YEAR exist (if running main execution)
# - available_raters() and load_corr_for_rater() exist (from your sourced utilities)
################################################################################

suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required")
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Package 'data.table' is required")
  library(ggplot2)
  library(data.table)
})

# ---- Metrics ----

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

compute_correlation_quantiles <- function(corr_list, probs = c(0.25, 0.5, 0.75, 0.9)) {
  qmat <- matrix(NA_real_, nrow = length(corr_list), ncol = length(probs))
  colnames(qmat) <- paste0("q", probs * 100)
  
  for (i in seq_along(corr_list)) {
    v <- extract_all_abs_cors(corr_list[[i]])
    if (length(v) > 0) qmat[i, ] <- quantile(v, probs = probs, na.rm = TRUE)
  }
  
  data.frame(
    quantile   = paste0("Q", round(probs * 100)),
    prob       = probs,
    median_rho = apply(qmat, 2, median,   na.rm = TRUE),
    q025       = apply(qmat, 2, quantile, probs = 0.025, na.rm = TRUE),
    q975       = apply(qmat, 2, quantile, probs = 0.975, na.rm = TRUE),
    row.names  = NULL
  )
}

# ---- Geometry ----

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

# ---- Plot builder (Option C highlighting) ----

create_saturn_plot <- function(
    base_out_dir      = BASE_OUT_DIR,
    base_viz_dir      = BASE_VIZ_DIR,
    year              = YEAR,
    prob              = 0.5,
    quantiles         = c(0.1, 0.3, 0.5, 0.7, 0.9),  # Only quantiles that will be highlighted

    highlight_gss     = TRUE,

    # Option C parameters:
    relevant_q        = c("Q90", "Q75", "Q50"),  # which quantiles to consider "relevant"
    delta_min         = 0.10,            # highlight if rho_model - rho_GSS >= delta_min
    top_n_per_q       = NULL,            # optional: keep only top-N by delta per quantile
    
    # Styling:
    spaghetti_alpha   = 0.12,
    spaghetti_width   = 0.35,
    gss_width         = 1.2,
    highlight_width   = 0.8,
    
    save_pdf          = TRUE,
    output_file       = NULL,
    save_summaries    = TRUE
) {
  
  message("\n=== Saturn Plot (Improved — Option C: Δ vs GSS) ===")
  message("prob: ", prob, " | quantiles: ", paste(quantiles, collapse = ", "))
  message("relevant_q: ", paste(relevant_q, collapse = ", "),
          " | delta_min: ", delta_min,
          if (!is.null(top_n_per_q)) paste0(" | top_n_per_q: ", top_n_per_q) else "")
  
  raters <- available_raters(base_out_dir = base_out_dir, year = year)
  if (!length(raters)) stop("No raters found in ", base_out_dir, " for year ", year)
  
  mac_summary <- data.table(
    rater      = character(),
    median_mac = numeric(),
    mean_mac   = numeric(),
    sd_mac     = numeric(),
    q025_mac   = numeric(),
    q975_mac   = numeric()
  )
  
  quantile_summary <- list()
  ellipse_data <- list()
  
  for (rater in raters) {
    message("Processing: ", rater)
    
    corr_list <- tryCatch(
      load_corr_for_rater(
        rater        = rater,
        base_out_dir = base_out_dir,
        year         = year,
        strict       = FALSE
      ),
      error = function(e) {
        warning("Could not load ", rater, ": ", e$message)
        NULL
      }
    )
    
    if (is.null(corr_list) || !length(corr_list)) next
    
    mac_vec <- compute_mac_bootstraps(corr_list)
    mac_vec <- mac_vec[!is.na(mac_vec)]
    if (!length(mac_vec)) next
    
    med_mac  <- median(mac_vec)
    mean_mac <- mean(mac_vec)
    sd_mac   <- sd(mac_vec)
    q025     <- as.numeric(quantile(mac_vec, 0.025))
    q975     <- as.numeric(quantile(mac_vec, 0.975))
    
    mac_summary <- rbindlist(list(
      mac_summary,
      data.table(
        rater      = rater,
        median_mac = med_mac,
        mean_mac   = mean_mac,
        sd_mac     = sd_mac,
        q025_mac   = q025,
        q975_mac   = q975
      )
    ))
    
    quant_stats <- compute_correlation_quantiles(corr_list, probs = quantiles)
    quant_stats$rater <- rater
    quantile_summary[[rater]] <- quant_stats
    
    # Ellipses for each quantile (based on median rho across bootstraps)
    for (q in seq_len(nrow(quant_stats))) {
      rho_q <- quant_stats$median_rho[q]
      qlab  <- quant_stats$quantile[q]
      
      xy <- ellipse_from_rho(rho = rho_q, prob = prob, n = 361L)
      ellipse_data[[paste(rater, qlab, sep = "_")]] <- data.frame(
        x        = xy[, "x"],
        y        = xy[, "y"],
        rater    = rater,
        quantile = qlab,
        prob_q   = quant_stats$prob[q],
        rho      = rho_q,
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (!length(ellipse_data)) {
    stop("No ellipse data created. Check that bootstrap correlation matrices exist and load correctly.")
  }
  
  ellipse_dt <- rbindlist(ellipse_data, use.names = TRUE, fill = TRUE)
  ellipse_dt[, rater := as.character(rater)]
  ellipse_dt[, is_gss := tolower(rater) == "gss"]
  
  # Order raters by constraint (median MAC)
  mac_summary <- mac_summary[order(-median_mac)]
  rater_order <- mac_summary$rater
  ellipse_dt[, rater := factor(rater, levels = rater_order)]
  
  # Quantile ordering (highest first)
  quantile_levels <- paste0("Q", round(sort(quantiles * 100, decreasing = TRUE)))
  ellipse_dt[, quantile := factor(quantile, levels = quantile_levels)]
  
  # Path grouping
  ellipse_dt[, group_id := paste(rater, quantile, sep = "_")]
  
  # ---------------------------------------------------------------------------
  # OPTION C: Determine which (model, quantile) pairs to highlight based on Δ vs GSS
  # ---------------------------------------------------------------------------
  
  # One row per (rater, quantile): rho is constant within each ellipse path
  rq <- ellipse_dt[, .(rho = rho[1], is_gss = is_gss[1]), by = .(rater, quantile)]
  
  # Ensure GSS exists
  if (!any(rq$is_gss)) {
    warning("No GSS rater found. Option C highlighting (Δ vs GSS) cannot run; no highlights will be applied.")
    ellipse_dt[, is_hi := FALSE]
  } else {
    # GSS rho per quantile
    gss_rho <- rq[is_gss == TRUE, .(quantile, gss_rho = rho)]
    
    # Join to non-GSS rows and compute delta
    rq2 <- rq[is_gss == FALSE][gss_rho, on = "quantile"]
    rq2[, delta := rho - gss_rho]
    
    # Candidate highlights by threshold within relevant quantiles
    hi_pairs <- rq2[quantile %in% relevant_q & !is.na(delta) & delta >= delta_min,
                    .(rater, quantile, delta)]
    
    # Optional: keep only top-N deltas per quantile
    if (!is.null(top_n_per_q) && nrow(hi_pairs) > 0) {
      top_n_per_q <- as.integer(top_n_per_q)
      hi_pairs <- hi_pairs[order(-delta)]
      hi_pairs <- hi_pairs[, .SD[1:min(.N, top_n_per_q)], by = quantile]
    }
    
    # Mark highlights back onto full ellipse_dt (safe, no recycling issues)
    ellipse_dt[, is_hi := FALSE]
    if (nrow(hi_pairs) > 0) {
      ellipse_dt[hi_pairs[, .(rater, quantile)], is_hi := TRUE, on = .(rater, quantile)]
    }
    
    message("Highlighted pairs (count): ", sum(ellipse_dt$is_hi))
    if (nrow(hi_pairs) > 0) {
      message("Example highlights:")
      print(hi_pairs[order(-delta)][1:min(.N, 12)])
    } else {
      message("No (model, quantile) pairs met delta_min in relevant_q; no highlights.")
    }
  }
  
  # Layer subsets
  gss_dt  <- ellipse_dt[ is_gss == TRUE]
  hi_dt   <- ellipse_dt[ is_gss == FALSE & is_hi == TRUE]
  rest_dt <- ellipse_dt[ is_gss == FALSE & is_hi == FALSE]
  
  # Reference circle for rho=0 at same prob
  ref <- reference_circle(prob = prob, n = 361L)
  
  # Linetypes by quantile (small legend)
  linetypes <- c("Q90"="solid", "Q75"="solid", "Q50"="dashed", "Q25"="dotted")
  if (!all(levels(ellipse_dt$quantile) %in% names(linetypes))) {
    lt <- rep(c("solid", "dashed", "dotted", "dotdash"), length.out = length(levels(ellipse_dt$quantile)))
    names(lt) <- levels(ellipse_dt$quantile)
    linetypes <- lt
  }
  
  # Colors only for highlighted raters (keeps legend small)
  hi_levels <- unique(as.character(hi_dt$rater))
  hi_levels <- hi_levels[order(match(hi_levels, rater_order))]
  hi_palette <- NULL
  if (length(hi_levels)) {
    hi_cols <- grDevices::hcl.colors(length(hi_levels), palette = "Dark 3")
    hi_palette <- setNames(hi_cols, hi_levels)
  }
  
  # Plot limits
  r0 <- sqrt(qchisq(prob, df = 2))
  lim <- 1.25 * r0
  
  # --- FACETED PLOT: Each quantile in its own panel ---

  p <- ggplot() +
    geom_hline(yintercept = 0, color = "gray90", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "gray90", linewidth = 0.3) +
    geom_path(
      data = ref, aes(x, y),
      color = "gray80", linewidth = 0.4, linetype = "dashed", lineend = "round"
    ) +

    # Non-highlighted LLM spaghetti (neutral) - no linetype aesthetic needed per facet
    geom_path(
      data = rest_dt,
      aes(x, y, group = group_id),
      color = "gray50",
      alpha = spaghetti_alpha,
      linewidth = spaghetti_width,
      lineend = "round"
    ) +

    # Highlighted (model, quantile) paths - no linetype needed
    { if (nrow(hi_dt)) geom_path(
      data = hi_dt,
      aes(x, y, group = group_id, color = rater),
      alpha = 0.65,
      linewidth = highlight_width,
      lineend = "round"
    )
    } +

    # GSS on top
    { if (nrow(gss_dt)) geom_path(
      data = gss_dt,
      aes(x, y, group = group_id),
      color = if (highlight_gss) "black" else "gray20",
      alpha = if (highlight_gss) 1.0 else 0.4,
      linewidth = if (highlight_gss) gss_width else 0.6,
      lineend = "round"
    )
    } +

    # Facet by quantile (each panel shows one constraint level)
    facet_wrap(~ quantile, nrow = 3, labeller = labeller(quantile = function(x) {
      paste0("Constraint: ", x)
    })) +

    coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim), expand = FALSE) +
    labs(x = "Standardized Dimension 1", y = "Standardized Dimension 2") +
    theme_classic(base_size = 14) +
    theme(
      # Legend at bottom (no need for linetype legend with facets)
      legend.position   = "bottom",
      legend.direction  = "horizontal",
      legend.box.just   = "center",

      legend.title      = element_text(face = "bold", size = 18),
      legend.text       = element_text(size = 16),
      panel.grid        = element_blank(),

      # Facet strip styling
      strip.background  = element_rect(fill = "gray95", color = "gray70", linewidth = 0.5),
      strip.text        = element_text(face = "bold", size = 12),

      # Make legend block fit on the page
      legend.margin     = margin(t = 6, r = 6, b = 6, l = 6),
      legend.box.margin = margin(t = 6),
      plot.margin       = margin(10, 10, 18, 10),

      # Dense model legends: tighter keys
      legend.key.height = unit(0.9, "lines"),
      legend.key.width  = unit(1.6, "lines"),

      # Panel spacing
      panel.spacing     = unit(1, "lines")
    )

  # Guides / scales:
  # - Only show color legend for highlighted models (no linetype legend needed with facets)
  if (length(hi_levels)) {
    p <- p +
      scale_color_manual(values = hi_palette, name = "Highlighted models") +
      guides(
        color = guide_legend(
          ncol = 3,          # More columns for horizontal layout across facets
          byrow = TRUE,
          title.position = "top",
          title.hjust = 0.5,
          override.aes = list(alpha = 1, linewidth = 1.2)
        )
      )
  } else {
    p <- p +
      guides(color = "none")
  }
  
  # Save plot
  if (save_pdf) {
    if (is.null(output_file)) {
      mvn_dir <- file.path(base_viz_dir, sprintf("mvn_%d", year))
      output_file <- file.path(mvn_dir, "saturn_plot.pdf")
    }
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    message("Saving: ", output_file)
    
    # 2 cols × 3 rows layout for 6 quantile facets
    ggsave(filename = output_file, plot = p, width = 12, height = 20, device = "pdf")
  }
  
  # Save summaries
  if (save_summaries) {
    mac_file <- file.path(base_out_dir, sprintf("mac_summary_%d.csv", year))
    fwrite(mac_summary, mac_file)
    message("Saved MAC summary: ", mac_file)
    
    quantile_dt <- rbindlist(quantile_summary, use.names = TRUE, fill = TRUE)
    quant_file <- file.path(base_out_dir, sprintf("correlation_quantiles_%d.csv", year))
    fwrite(quantile_dt, quant_file)
    message("Saved correlation quantiles: ", quant_file)
  }
  
  invisible(p)
}


# ---- Facet Grid Overview (Q95 → Q5) ----

create_saturn_facet_grid <- function(
    base_out_dir      = BASE_OUT_DIR,
    base_viz_dir      = BASE_VIZ_DIR,
    year              = YEAR,
    prob              = 0.5,
    quantile_seq      = seq(0.95, 0.05, by = -0.05),  # Q95 → Q5

    highlight_gss     = TRUE,
    llm_color         = "steelblue",  # All LLMs same color
    llm_alpha         = 0.3,
    llm_width         = 0.4,
    gss_width         = 1.2,

    save_pdf          = TRUE,
    output_file       = NULL
) {

  message("\n=== Saturn Facet Grid (Q95 → Q5) - Ultra Minimal ===")
  message("Quantile range: ", min(quantile_seq), " to ", max(quantile_seq))
  message("Number of facets: ", length(quantile_seq))

  # Get all raters
  raters <- available_raters(base_out_dir = base_out_dir, year = year)
  message("Found ", length(raters), " raters")

  if (length(raters) == 0) {
    stop("No raters found in ", base_out_dir)
  }

  # Load correlation matrices and compute quantiles
  all_quantile_data <- list()

  for (rater in raters) {
    message("Processing: ", rater)

    corr_list <- tryCatch(
      load_corr_for_rater(rater = rater, base_out_dir = base_out_dir, year = year, strict = FALSE),
      error = function(e) NULL
    )

    if (!is.null(corr_list) && length(corr_list) > 0) {
      quant_stats <- compute_correlation_quantiles(corr_list, probs = quantile_seq)
      quant_stats$rater <- rater
      all_quantile_data[[rater]] <- quant_stats
    }
  }

  if (length(all_quantile_data) == 0) {
    stop("No correlation data loaded successfully")
  }

  quantile_dt <- rbindlist(all_quantile_data)
  quantile_dt[, is_gss := tolower(rater) == "gss"]

  # Generate ellipse data for all (rater, quantile) pairs
  message("\nGenerating ellipse paths...")

  ellipse_list <- list()

  for (i in 1:nrow(quantile_dt)) {
    rater_i <- quantile_dt$rater[i]
    quant_i <- quantile_dt$quantile[i]
    rho_i   <- quantile_dt$median_rho[i]
    is_gss_i <- quantile_dt$is_gss[i]

    xy <- ellipse_from_rho(rho = rho_i, prob = prob, n = 361L)

    ellipse_list[[paste(rater_i, quant_i, sep = "_")]] <- data.frame(
      x          = xy[, "x"],
      y          = xy[, "y"],
      rater      = rater_i,
      quantile   = quant_i,
      rho        = rho_i,
      is_gss     = is_gss_i,
      stringsAsFactors = FALSE
    )
  }

  ellipse_dt <- rbindlist(ellipse_list)
  ellipse_dt[, group_id := paste(rater, quantile, sep = "_")]

  # Order quantiles for facets (Q95 → Q5)
  quantile_levels <- paste0("Q", round(quantile_seq * 100))
  ellipse_dt[, quantile := factor(quantile, levels = quantile_levels)]

  # Separate GSS and LLMs
  gss_dt <- ellipse_dt[is_gss == TRUE]
  llm_dt <- ellipse_dt[is_gss == FALSE]

  # Plot limits
  r0 <- sqrt(qchisq(prob, df = 2))
  lim <- 1.25 * r0

  # Create ultra-minimal faceted plot
  message("\nCreating ultra-minimal facet strip...")

  p <- ggplot() +
    # All LLMs in same color
    geom_path(
      data = llm_dt,
      aes(x, y, group = group_id),
      color = llm_color,
      alpha = llm_alpha,
      linewidth = llm_width,
      lineend = "round"
    ) +

    # GSS in black
    {if (nrow(gss_dt) > 0) geom_path(
      data = gss_dt,
      aes(x, y, group = group_id),
      color = "black",
      alpha = 1.0,
      linewidth = gss_width,
      lineend = "round"
    )} +

    # All facets in one row, no labels
    facet_wrap(~ quantile, nrow = 1) +

    coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim), expand = FALSE) +

    # Ultra-minimal theme - nothing but the contours
    theme_void() +
    theme(
      # No facet labels
      strip.text = element_blank(),
      strip.background = element_blank(),

      # No legend
      legend.position = "none",

      # Minimal margins
      plot.margin = margin(2, 2, 2, 2),
      panel.spacing = unit(0.1, "lines")
    )

  # Save plot
  if (save_pdf) {
    if (is.null(output_file)) {
      mvn_dir <- file.path(base_viz_dir, sprintf("mvn_%d", year))
      dir.create(mvn_dir, recursive = TRUE, showWarnings = FALSE)
      output_file <- file.path(mvn_dir, "saturn_facet_grid.pdf")
    }

    # Calculate dimensions for horizontal strip (all in one row)
    n_facets <- length(quantile_seq)

    # Width: ~2 inches per facet for horizontal strip
    width <- n_facets * 2
    # Height: fixed at 3 inches for single row
    height <- 3

    message("\nSaving facet grid to: ", output_file)
    message("Dimensions: ", width, " x ", height, " inches (", n_facets, " facets in one row)")

    ggsave(
      filename = output_file,
      plot = p,
      width = width,
      height = height,
      device = "pdf"
    )
  }

  message("\n=== Facet Grid Complete ===\n")

  invisible(p)
}


################################################################################
# Belief Pair Facet Grid — All pairwise correlations as individual panels
################################################################################

#' Extract pairwise correlations from a correlation matrix
extract_pairwise_cors <- function(R) {
  if (is.null(R)) return(NULL)
  R <- as.matrix(R)
  vars <- colnames(R)
  if (is.null(vars) || length(vars) < 2) return(NULL)

  pairs <- which(upper.tri(R), arr.ind = TRUE)

  data.table(
    var1 = vars[pairs[, 1]],
    var2 = vars[pairs[, 2]],
    rho  = R[pairs]
  )
}

#' Compute median pairwise correlations across bootstrap samples
compute_median_pairwise_cors <- function(corr_list) {
  all_pairs <- rbindlist(lapply(corr_list, extract_pairwise_cors), fill = TRUE)
  if (nrow(all_pairs) == 0) return(NULL)

  median_pairs <- all_pairs[, .(median_rho = median(rho, na.rm = TRUE)),
                            by = .(var1, var2)]

  median_pairs[, pair_label := mapply(function(a, b) {
    paste(sort(c(a, b)), collapse = " × ")
  }, var1, var2)]

  median_pairs
}

#' Create belief pair facet grid
#'
#' Each panel shows one belief pair (var1 × var2) with LLMs in faint grey and GSS in black.
#' @param pairs_per_page Number of pairs per page (NULL = all on one page)
#' @param page_ncol Number of columns per page (default 8)
#' @param page_nrow Number of rows per page (default 6, so 48 pairs per page)
#' @param compact_mode If TRUE, use compact layout for 3 landscape pages (25x18 grid)
#' @param cell_size Size multiplier per cell (default 1.4, use smaller for compact)
create_belief_pair_grid <- function(
    base_out_dir      = BASE_OUT_DIR,
    base_viz_dir      = BASE_VIZ_DIR,
    year              = YEAR,
    prob              = 0.5,
    vars_subset       = NULL,
    llm_alpha         = 0.15,
    llm_width         = 0.3,
    llm_color         = "gray50",
    gss_width         = 0.8,
    gss_color         = "black",
    ncol              = NULL,
    page_ncol         = 8,
    page_nrow         = 6,
    pairs_per_page    = NULL,
    compact_mode      = FALSE,
    cell_size         = 1.4,
    save_pdf          = TRUE,
    output_file       = NULL
) {

  message("\n=== Belief Pair Facet Grid ===")
  message("Generating ellipses for all pairwise correlations")

  raters <- available_raters(base_out_dir = base_out_dir, year = year)
  if (!length(raters)) stop("No raters found in ", base_out_dir, " for year ", year)

  message("Found ", length(raters), " raters")

  all_rater_pairs <- list()

  for (rater in raters) {
    message("Processing: ", rater)

    corr_list <- tryCatch(
      load_corr_for_rater(rater = rater, base_out_dir = base_out_dir,
                          year = year, strict = FALSE),
      error = function(e) { warning("Could not load ", rater, ": ", e$message); NULL }
    )

    if (is.null(corr_list) || !length(corr_list)) next

    median_pairs <- compute_median_pairwise_cors(corr_list)
    if (is.null(median_pairs) || nrow(median_pairs) == 0) next

    median_pairs[, rater := rater]
    median_pairs[, is_gss := tolower(rater) == "gss"]

    all_rater_pairs[[rater]] <- median_pairs
  }

  if (!length(all_rater_pairs)) stop("No pairwise correlation data loaded")

  pairs_dt <- rbindlist(all_rater_pairs, use.names = TRUE, fill = TRUE)

  if (!is.null(vars_subset)) {
    pairs_dt <- pairs_dt[var1 %in% vars_subset & var2 %in% vars_subset]
    message("Filtered to ", length(unique(pairs_dt$pair_label)), " pairs from ",
            length(vars_subset), " variables")
  }

  n_pairs <- length(unique(pairs_dt$pair_label))
  message("Total pairs: ", n_pairs)

  if (n_pairs == 0) stop("No pairs remaining after filtering")

  message("\nGenerating ellipse paths...")

  ellipse_list <- list()

  for (i in 1:nrow(pairs_dt)) {
    rater_i  <- pairs_dt$rater[i]
    pair_i   <- pairs_dt$pair_label[i]
    rho_i    <- pairs_dt$median_rho[i]
    is_gss_i <- pairs_dt$is_gss[i]

    if (is.na(rho_i)) next

    xy <- ellipse_from_rho(rho = rho_i, prob = prob, n = 181L)

    ellipse_list[[paste(rater_i, pair_i, sep = "___")]] <- data.frame(
      x = xy[, "x"], y = xy[, "y"],
      rater = rater_i, pair_label = pair_i,
      rho = rho_i, is_gss = is_gss_i,
      stringsAsFactors = FALSE
    )
  }

  ellipse_dt <- rbindlist(ellipse_list)
  ellipse_dt[, group_id := paste(rater, pair_label, sep = "___")]

  # Order pairs by median LLM correlation: high positive → zero → high negative
  # (thin diagonal-right → circular → thin diagonal-left)
  llm_pair_rho <- ellipse_dt[is_gss == FALSE, .(median_llm_rho = median(rho, na.rm = TRUE)),
                              by = pair_label]
  llm_pair_rho <- llm_pair_rho[order(-median_llm_rho)]  # descending: positive first
  pair_levels <- llm_pair_rho$pair_label

  # Filter to only pairs that have LLM data (exclude GSS-only pairs)
  ellipse_dt <- ellipse_dt[pair_label %in% pair_levels]

  ellipse_dt[, pair_label := factor(pair_label, levels = pair_levels)]

  gss_dt <- ellipse_dt[is_gss == TRUE]
  llm_dt <- ellipse_dt[is_gss == FALSE]

  n_pairs <- length(unique(ellipse_dt$pair_label))  # update after filtering
  message("Pairs with LLM data: ", n_pairs)
  message("Ordering: high positive rho → zero → high negative rho")

  r0 <- sqrt(qchisq(prob, df = 2))
  lim <- 1.3 * r0

  # Compact mode: dense grid with readable labels (~10 wide landscape pages)

  if (compact_mode) {
    page_ncol <- 20
    page_nrow <- 9
    cell_size <- 0.95
    message("\nCompact mode: targeting ~10 landscape pages (20 cols × 9 rows = 180 pairs/page)")
  }

  # Determine pagination
  if (is.null(pairs_per_page)) {
    pairs_per_page <- page_ncol * page_nrow
  }

  # Split pairs into pages
  n_pages <- ceiling(n_pairs / pairs_per_page)

  if (n_pages == 1) {
    # Single page - use original behavior
    if (is.null(ncol)) ncol <- ceiling(sqrt(n_pairs * 1.5))
    message("\nCreating single-page grid with ", ncol, " columns...")
  } else {
    ncol <- page_ncol
    message("\nCreating ", n_pages, " pages with ", page_ncol, " cols x ", page_nrow, " rows each...")
  }

  # Build list of plots (one per page)
  plots <- list()

  for (page in seq_len(n_pages)) {
    start_idx <- (page - 1) * pairs_per_page + 1
    end_idx <- min(page * pairs_per_page, n_pairs)
    page_pairs <- pair_levels[start_idx:end_idx]

    page_llm_dt <- llm_dt[pair_label %in% page_pairs]
    page_gss_dt <- gss_dt[pair_label %in% page_pairs]

    # Re-factor to only include this page's pairs (keeps facet order)
    page_llm_dt[, pair_label := factor(pair_label, levels = page_pairs)]
    page_gss_dt[, pair_label := factor(pair_label, levels = page_pairs)]

    p <- ggplot() +
      geom_path(data = page_llm_dt, aes(x, y, group = group_id),
                color = llm_color, alpha = llm_alpha,
                linewidth = if (compact_mode) 0.2 else llm_width, lineend = "round") +
      {if (nrow(page_gss_dt) > 0) geom_path(data = page_gss_dt, aes(x, y, group = group_id),
                color = gss_color, alpha = 1.0,
                linewidth = if (compact_mode) 0.5 else gss_width, lineend = "round")} +
      facet_wrap(~ pair_label, ncol = ncol) +
      coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim), expand = FALSE) +
      theme_void() +
      theme(
        strip.text = element_text(size = if (compact_mode) 4 else 7, margin = margin(1, 0, 1, 0)),
        strip.background = element_blank(),
        legend.position = "none",
        plot.margin = if (compact_mode) margin(5, 5, 5, 5) else margin(10, 10, 10, 10),
        panel.spacing = unit(if (compact_mode) 0.15 else 0.3, "lines")
      )

    plots[[page]] <- p
  }

  if (save_pdf) {
    mvn_dir <- file.path(base_viz_dir, sprintf("mvn_%d", year))
    dir.create(mvn_dir, recursive = TRUE, showWarnings = FALSE)

    # Calculate dimensions per page (landscape for compact mode)
    width  <- page_ncol * cell_size + 0.5
    height <- page_nrow * cell_size + 0.5

    for (page in seq_len(n_pages)) {
      if (n_pages == 1) {
        out_file <- file.path(mvn_dir, "belief_pair_grid.pdf")
      } else if (compact_mode) {
        out_file <- file.path(mvn_dir, sprintf("belief_pair_grid_compact_%d.pdf", page))
      } else {
        out_file <- file.path(mvn_dir, sprintf("belief_pair_grid%d.pdf", page))
      }

      start_idx <- (page - 1) * pairs_per_page + 1
      end_idx <- min(page * pairs_per_page, n_pairs)
      n_page_pairs <- end_idx - start_idx + 1

      message("Saving page ", page, "/", n_pages, ": ", out_file,
              " (pairs ", start_idx, "-", end_idx, ")")

      ggsave(filename = out_file, plot = plots[[page]], width = width, height = height,
             device = "pdf", limitsize = FALSE)
    }

    orientation <- if (width > height) "landscape" else "portrait"
    message("\nDimensions per page: ", round(width, 1), " x ", round(height, 1), " inches (", orientation, ")")
  }

  message("\n=== Belief Pair Grid Complete ===\n")
  invisible(plots)
}


################################################################################
# MAIN EXECUTION
################################################################################

if (exists("BASE_OUT_DIR") && exists("BASE_VIZ_DIR") && exists("YEAR")) {
  message("\n====================================================================")
  message("Saturn Plot (Faceted): Option C Highlighting (Δ vs GSS at quantiles)")
  message("Each quantile shown in separate panel for easy comparison")
  message("====================================================================\n")
  
  saturn_plot <- create_saturn_plot(
    base_out_dir   = BASE_OUT_DIR,
    base_viz_dir   = BASE_VIZ_DIR,
    year           = YEAR,
    prob           = 0.5,
    quantiles      = c(0.15, 0.3, 0.5, 0.65, 0.8, 0.95),  # Six quantiles

    # Option C controls
    relevant_q     = c("Q95", "Q80", "Q65", "Q50", "Q30", "Q15"),
    delta_min      = 0.10,    # increase to highlight fewer, decrease to highlight more
    top_n_per_q    = 6L,      # set NULL to keep all above threshold
    
    highlight_gss  = TRUE,
    save_pdf       = TRUE
  )

  # Optional: Facet grid overview (Q95 → Q5) - ultra-minimal horizontal strip
  # Uncomment to generate wide horizontal strip showing constraint evolution:
  saturn_facet_grid <- create_saturn_facet_grid(
    base_out_dir   = BASE_OUT_DIR,
    base_viz_dir   = BASE_VIZ_DIR,
    year           = YEAR,
    quantile_seq   = seq(0.95, 0.05, by = -0.05),  # 19 quantiles, all in one row
    llm_color      = "steelblue",
    llm_alpha      = 0.3,
    save_pdf       = TRUE
  )

  # Belief pair grid - COMPACT VERSION for appendix (~10 landscape pages)
  # With readable labels at 15×12 = 180 pairs per page
  belief_pair_grid_compact <- create_belief_pair_grid(
    base_out_dir   = BASE_OUT_DIR,
    base_viz_dir   = BASE_VIZ_DIR,
    year           = YEAR,
    prob           = 0.5,
    vars_subset    = NULL,
    llm_alpha      = 0.15,
    llm_color      = "gray50",
    gss_color      = "black",
    compact_mode   = TRUE,      # 3 landscape pages (25 cols × 18 rows)
    save_pdf       = TRUE
  )

} else {
  message("BASE_OUT_DIR, BASE_VIZ_DIR, and YEAR must be defined to run.")
}