################################################################################
# Saturn Plots (Improved) — Option C Highlighting
# Highlight models that are meaningfully more constrained than GSS at “relevant”
# quantiles, using Δ = rho_model(quantile) - rho_GSS(quantile).
#
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
    quantile   = paste0("Q", probs * 100),
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
    quantiles         = c(0.25, 0.5, 0.75, 0.9),
    
    highlight_gss     = TRUE,
    
    # Option C parameters:
    relevant_q        = c("Q90", "Q75"),  # which quantiles to consider “relevant”
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
  quantile_levels <- paste0("Q", sort(quantiles * 100, decreasing = TRUE))
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
  
  # --- PLOT (legend fix: stack legends vertically + wrap model legend into columns) ---
  
  p <- ggplot() +
    geom_hline(yintercept = 0, color = "gray90", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "gray90", linewidth = 0.3) +
    geom_path(
      data = ref, aes(x, y),
      color = "gray80", linewidth = 0.4, linetype = "dashed", lineend = "round"
    ) +
    
    # Non-highlighted LLM spaghetti (neutral)
    geom_path(
      data = rest_dt,
      aes(x, y, group = group_id, linetype = quantile),
      color = "gray50",
      alpha = spaghetti_alpha,
      linewidth = spaghetti_width,
      lineend = "round"
    ) +
    
    # Highlighted (model, quantile) paths
    { if (nrow(hi_dt)) geom_path(
      data = hi_dt,
      aes(x, y, group = group_id, linetype = quantile, color = rater),
      alpha = 0.65,
      linewidth = highlight_width,
      lineend = "round"
    )
    } +
    
    # GSS on top
    { if (nrow(gss_dt)) geom_path(
      data = gss_dt,
      aes(x, y, group = group_id, linetype = quantile),
      color = if (highlight_gss) "black" else "gray20",
      alpha = if (highlight_gss) 1.0 else 0.4,
      linewidth = if (highlight_gss) gss_width else 0.6,
      lineend = "round"
    )
    } +
    
    scale_linetype_manual(values = linetypes, name = "Constraint quantile") +
    coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim), expand = FALSE) +
    labs(x = "Standardized Dimension 1", y = "Standardized Dimension 2") +
    theme_classic(base_size = 16) +
    theme(
      # Put legends at bottom, but STACK them vertically so one can't push the other off
      legend.position   = "bottom",
      legend.box        = "vertical",
      legend.direction  = "horizontal",
      legend.box.just   = "center",
      
      legend.title      = element_text(face = "bold"),
      panel.grid        = element_blank(),
      
      # Make legend block fit on the page (prevents clipping)
      legend.margin     = margin(t = 6, r = 6, b = 6, l = 6),
      legend.box.margin = margin(t = 6),
      plot.margin       = margin(10, 10, 18, 10),
      
      # Dense model legends: tighter keys
      legend.key.height = unit(0.9, "lines"),
      legend.key.width  = unit(1.6, "lines")
    )
  
  # Guides / scales:
  # - If highlights exist: show a wrapped color legend (multiple columns) + quantile legend
  # - If no highlights: hide color legend, keep quantile legend
  if (length(hi_levels)) {
    p <- p +
      scale_color_manual(values = hi_palette, name = "Highlighted models") +
      guides(
        color = guide_legend(
          order = 1,
          ncol = 3,          # <-- change this: fewer columns = taller; more columns = wider
          byrow = TRUE,
          title.position = "top",
          title.hjust = 0.5,
          override.aes = list(alpha = 1, linewidth = 1.2)
        ),
        linetype = guide_legend(
          order = 2,
          nrow = 1,
          title.position = "top",
          title.hjust = 0.5,
          override.aes = list(alpha = 1, linewidth = 1.0)
        )
      )
  } else {
    p <- p +
      guides(
        color = "none",
        linetype = guide_legend(
          order = 1,
          nrow = 1,
          title.position = "top",
          title.hjust = 0.5,
          override.aes = list(alpha = 1, linewidth = 1.0)
        )
      )
  }
  
  # Save plot
  if (save_pdf) {
    if (is.null(output_file)) {
      mvn_dir <- file.path(base_viz_dir, sprintf("mvn_%d", year))
      output_file <- file.path(mvn_dir, "saturn_plot.pdf")
    }
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    message("Saving: ", output_file)
    
    # Slightly taller device gives legends more room (prevents clipping)
    ggsave(filename = output_file, plot = p, width = 12, height = 14, device = "pdf")
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


################################################################################
# MAIN EXECUTION
################################################################################

if (exists("BASE_OUT_DIR") && exists("BASE_VIZ_DIR") && exists("YEAR")) {
  message("\n====================================================================")
  message("Saturn Plot (Improved): Option C Highlighting (Δ vs GSS at quantiles)")
  message("====================================================================\n")
  
  saturn_plot <- create_saturn_plot(
    base_out_dir   = BASE_OUT_DIR,
    base_viz_dir   = BASE_VIZ_DIR,
    year           = YEAR,
    prob           = 0.5,
    quantiles      = c(0.25, 0.5, 0.75, 0.9),
    
    # Option C controls
    relevant_q     = c("Q90", "Q75", "Q50"),
    delta_min      = 0.10,    # increase to highlight fewer, decrease to highlight more
    top_n_per_q    = 6L,      # set NULL to keep all above threshold
    
    highlight_gss  = TRUE,
    save_pdf       = TRUE
  )
  
} else {
  message("BASE_OUT_DIR, BASE_VIZ_DIR, and YEAR must be defined to run.")
}
