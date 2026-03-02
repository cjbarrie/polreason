################################################################################
# Constraint diagnostics: PC1 variance share + effective dependence (De)
#
# Overview
# --------
# This script computes and visualizes two scalar “constraint” diagnostics from
# bootstrap polychoric correlation matrices for multiple raters (GSS + LLMs),
# split by education group (overall / low edu / high edu):
#
#   (A) PC1 variance explained:
#       share = λ1 / sum(λj)
#
#   (B) Effective dependence (Peña & Rodríguez, 2003):
#       De = 1 - |R|^(1/k)
#       (computed stably via the geometric mean of eigenvalues)
#
# Key design choices
# ------------------
#   - LLM-eligibility restriction for GSS:
#       To avoid anchoring the analysis on GSS-only items, GSS correlation
#       matrices are restricted (per education group) to the union of variables
#       that appear in at least one LLM bootstrap draw. This restriction is
#       performed by restrict_corr_list_to_vars(), which lives in utils.
#
#   - Bootstrap-first estimation:
#       All diagnostics are computed per bootstrap draw and then summarized
#       (means by default) to obtain point estimates, rather than computing
#       diagnostics on an averaged correlation matrix.
#
# Workflow
# --------
#   1) Discover raters on disk for YEAR via available_raters().
#   2) For each education group, compute the union of variables appearing in at
#      least one LLM draw (LLM var-union).
#   3) For each rater × edu_group × bootstrap draw:
#        - load the correlation matrix list (per subgroup file),
#        - if rater == "gss", restrict to the LLM var-union for that subgroup,
#        - compute PC1 variance share, effective dependence, and matrix dimension.
#   4) Summarize bootstrap distributions to point estimates per rater × edu_group.
#   5) Save bootstrap-level and point-estimate tables to BASE_OUT_DIR.
#   6) Produce:
#        - horizontal half-violin PDFs for each statistic (by edu group),
#        - a ggplot scatter of point estimates: De vs PC1 share (by edu group).
#
# Prerequisites (provided elsewhere in the project)
# ------------------------------------------------
# This script assumes the following objects/functions exist in the workspace:
#   - BASE_OUT_DIR : directory containing "<rater>-<YEAR>" subdirectories
#   - BASE_VIZ_DIR : directory to write PDF outputs
#   - available_raters(base_out_dir, year, exclude=...) -> character vector
#   - load_corr_for_rater(rater, base_out_dir, year, suffix, align_to_mode_varset,
#                         strict) -> list of correlation matrices (or NULL)
#   - restrict_corr_list_to_vars(corr_list, vars_keep) -> list(corr_list, bootstrap_id)
#   - pc1_var_explained(R)                # returns NA_real_ if invalid
#   - effective_dependence(R)             # returns NA_real_ if invalid
#
#
# Packages used
# -------------
#   - data.table
#   - ggplot2
#   - scales
#   - grid
#
# Output
# ------
#   - BASE_OUT_DIR/constraint_bootstrap_stats_<YEAR>.csv
#   - BASE_OUT_DIR/constraint_point_stats_<YEAR>.csv
#   - BASE_VIZ_DIR/constraint_violins_<YEAR>/
#       * constraint_pc1_var_explained_by_edu_<YEAR>.pdf
#       * constraint_effective_dependence_by_edu_<YEAR>.pdf
#       * constraint_pc1_vs_effective_dependence_by_edu.pdf
################################################################################


#' Plot horizontal half-violin distributions by rater and education group
#'
#' Produces a PDF with horizontal “half-violin” plots for a scalar statistic
#' (e.g., PC1 variance share or effective dependence), split by education group:
#' \code{overall}, \code{loedu}, \code{hiedu}.
#'
#' Each rater appears on a separate row (Y-axis). For each rater, the function
#' overlays half-violins for low-edu, high-edu, and overall groups.
#'
#' If a \code{bootstrap_id} column is present, the function also computes the
#' bootstrap distribution of deltas (\code{hiedu - loedu}) and annotates each
#' rater with median and a 95% bootstrap interval.
#'
#' @param dt A data.frame/data.table containing at least \code{rater},
#'   \code{edu_group}, and \code{value_col}. Optionally \code{bootstrap_id}.
#' @param value_col Character scalar giving the column name to plot.
#' @param xlab X-axis label.
#' @param main Plot title.
#' @param out_file Output PDF path.
#' @param gss_overall,gss_loedu,gss_hiedu Fill colors for GSS (overall/low/high).
#' @param llm_overall,llm_loedu,llm_hiedu Fill colors for non-GSS raters.
#'
#' @return Invisibly \code{NULL}. Called for side-effect (writes PDF).
#'
#' @keywords internal
plot_half_violins_edu <- function(
    dt,
    value_col,
    xlab,
    main,
    out_file,
    lolli_pad_frac   = -0.01,
    lolli_width_frac = 0.14,
    gss_overall  = "#FDB863",
    gss_loedu    = "#FEE0B6",
    gss_hiedu    = "#E08214",
    llm_overall  = "#BDBDBD",
    llm_loedu    = "#A6CEE3",
    llm_hiedu    = "#1F78B4"
) {
  stopifnot(is.data.frame(dt))
  
  if (!("rater" %in% names(dt))) stop("plot_half_violins_edu() requires a 'rater' column.")
  if (!("edu_group" %in% names(dt))) stop("plot_half_violins_edu() requires an 'edu_group' column.")
  if (!(value_col %in% names(dt))) stop("Column ", sQuote(value_col), " not found in data.")
  
  if (!inherits(dt, "data.table")) data.table::setDT(dt)
  
  dt[, edu_group := as.character(edu_group)]
  dt[is.na(edu_group) | edu_group == "", edu_group := "overall"]
  dt <- dt[edu_group %in% c("overall", "loedu", "hiedu")]
  
  dt <- dt[is.finite(get(value_col))]
  if (!nrow(dt)) {
    warning("No finite values for ", value_col, "; skipping plot.")
    return(invisible(NULL))
  }
  dt[, (value_col) := as.numeric(get(value_col))]
  
  delta_stats <- NULL
  if ("bootstrap_id" %in% names(dt)) {
    delta_stats <- dt[
      edu_group %in% c("loedu", "hiedu") & !is.na(bootstrap_id),
      {
        wide <- data.table::dcast(.SD, bootstrap_id ~ edu_group, value.var = value_col)
        if (!all(c("loedu", "hiedu") %in% names(wide))) {
          return(data.table::data.table(delta_median = NA_real_, delta_lwr = NA_real_, delta_upr = NA_real_))
        }
        deltas <- wide$hiedu - wide$loedu
        deltas <- deltas[is.finite(deltas)]
        if (!length(deltas)) {
          return(data.table::data.table(delta_median = NA_real_, delta_lwr = NA_real_, delta_upr = NA_real_))
        }
        data.table::data.table(
          delta_median = stats::median(deltas, na.rm = TRUE),
          delta_lwr    = stats::quantile(deltas, 0.025, na.rm = TRUE),
          delta_upr    = stats::quantile(deltas, 0.975, na.rm = TRUE)
        )
      },
      by = rater
    ]
  }
  
  med_dt <- dt[edu_group == "overall", .(
    median_value = as.numeric(stats::median(get(value_col), na.rm = TRUE))
  ), by = rater]
  
  if (nrow(med_dt) < length(unique(dt$rater))) {
    missing_raters <- setdiff(unique(dt$rater), med_dt$rater)
    if (length(missing_raters)) {
      med_missing <- dt[rater %in% missing_raters, .(
        median_value = as.numeric(stats::median(get(value_col), na.rm = TRUE))
      ), by = rater]
      med_dt <- data.table::rbindlist(list(med_dt, med_missing), use.names = TRUE, fill = TRUE)
    }
  }
  
  med_dt <- med_dt[order(median_value)]
  raters_order <- med_dt$rater
  
  dt[, rater := factor(rater, levels = raters_order)]
  n_raters <- length(raters_order)
  
  pdf_width  <- 10
  pdf_height <- max(3, 0.3 * n_raters + 1)
  
  grDevices::pdf(out_file, width = pdf_width, height = pdf_height)
  on.exit(grDevices::dev.off(), add = TRUE)
  
  x_range <- range(dt[[value_col]], finite = TRUE)
  if (!is.finite(diff(x_range)) || diff(x_range) == 0) {
    x_range <- x_range + c(-0.5, 0.5)
  } else {
    pad     <- 0.05 * diff(x_range)
    x_range <- x_range + c(-pad, pad)
  }
  
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  
  par(mar = c(4, 18, 4, 8))
  
  plot(
    NA,
    xlim = x_range,
    ylim = c(0.5, n_raters + 0.9),
    xlab = xlab,
    ylab = "",
    yaxt = "n",
    main = main,
    bty  = "n"
  )
  
  usr <- par("usr")
  segments(usr[1], usr[3], usr[2], usr[3])  # bottom
  segments(usr[1], usr[4], usr[2], usr[4])  # top
  segments(usr[1], usr[3], usr[1], usr[4])  # left
  
  axis(2, at = seq_len(n_raters), labels = raters_order, las = 2)
  
  grid_ticks <- pretty(x_range)
  abline(v = grid_ticks, col = "grey92", lwd = 0.45)
  
  max_width <- 0.4
  par(xpd = NA)
  
  # Lollipop layout (to the right of the plot area)
  usr <- par("usr")
  x_right <- usr[2]
  lolli_pad    <- lolli_pad_frac * diff(usr[1:2])
  lolli_width  <- lolli_width_frac * diff(usr[1:2])
  lolli_center <- x_right + lolli_pad + (lolli_width / 2)
  
  max_abs_delta <- NA_real_
  if (!is.null(delta_stats) && nrow(delta_stats)) {
    max_abs_delta <- max(abs(delta_stats$delta_median), na.rm = TRUE)
  }
  if (!is.finite(max_abs_delta) || max_abs_delta == 0) {
    max_abs_delta <- 1
  }
  
  # Delta axis (center line only; no tick labels)
  segments(
    x0 = lolli_center,
    x1 = lolli_center,
    y0 = 0.5,
    y1 = n_raters + 0.5,
    col = grDevices::adjustcolor("grey70", alpha.f = 0.8),
    lwd = 0.8
  )
  
  draw_half <- function(vals, base_y, fill_col, median_col = fill_col) {
    vals <- vals[is.finite(vals)]
    if (length(vals) < 2L || length(unique(vals)) < 2L) return(invisible(NULL))
    
    dens <- stats::density(vals)
    scaled_y <- dens$y / max(dens$y, na.rm = TRUE) * max_width
    
    x_poly <- c(dens$x, rev(dens$x))
    y_poly <- c(rep(base_y, length(dens$x)), base_y + rev(scaled_y))
    
    polygon(x_poly, y_poly, col = grDevices::adjustcolor(fill_col, alpha.f = 0.6), border = NA)
    lines(dens$x, base_y + scaled_y, col = grDevices::adjustcolor("black", alpha.f = 0.5))
    
    med_val <- stats::median(vals, na.rm = TRUE)
    segments(
      x0 = med_val, x1 = med_val,
      y0 = base_y, y1 = base_y + max_width * 1.1,
      lwd = 1.8, col = median_col
    )
  }
  
  if (!is.null(delta_stats) && nrow(delta_stats)) {
    delta_stats <- delta_stats[rater %in% raters_order]
  }
  
  for (i in seq_along(raters_order)) {
    r    <- raters_order[i]
    dt_r <- dt[rater == r]
    
    is_gss <- identical(as.character(r), "gss")
    
    vals_overall <- dt_r[edu_group == "overall", get(value_col)]
    vals_lo      <- dt_r[edu_group == "loedu",   get(value_col)]
    vals_hi      <- dt_r[edu_group == "hiedu",   get(value_col)]
    
    col_overall <- if (is_gss) gss_overall else llm_overall
    col_loedu   <- if (is_gss) gss_loedu   else llm_loedu
    col_hiedu   <- if (is_gss) gss_hiedu   else llm_hiedu
    
    if (length(vals_lo)) {
      draw_half(vals_lo, i, col_loedu, median_col = if (is_gss) col_loedu else "skyblue")
    }
    if (length(vals_hi)) {
      draw_half(vals_hi, i, col_hiedu)
    }
    if (length(vals_overall)) {
      draw_half(vals_overall, i, col_overall)
    }
    
    if (!is.null(delta_stats) && nrow(delta_stats)) {
      drow <- delta_stats[rater == r]
      if (nrow(drow) == 1L && is.finite(drow$delta_median)) {
        delta_val <- drow$delta_median
        x_end <- lolli_center + (delta_val / max_abs_delta) * (lolli_width / 2)
        
        sig_delta <- is.finite(drow$delta_lwr) && is.finite(drow$delta_upr) &&
          (drow$delta_lwr > 0 || drow$delta_upr < 0)
        delta_col <- if (sig_delta) "black" else "grey60"
        
        if (is.finite(drow$delta_lwr) && is.finite(drow$delta_upr)) {
          x_lwr <- lolli_center + (drow$delta_lwr / max_abs_delta) * (lolli_width / 2)
          x_upr <- lolli_center + (drow$delta_upr / max_abs_delta) * (lolli_width / 2)
          segments(x0 = x_lwr, x1 = x_upr, y0 = i, y1 = i, lwd = 0.9, col = delta_col)
          segments(x0 = x_lwr, x1 = x_lwr, y0 = i - 0.08, y1 = i + 0.08, lwd = 0.9, col = delta_col)
          segments(x0 = x_upr, x1 = x_upr, y0 = i - 0.08, y1 = i + 0.08, lwd = 0.9, col = delta_col)
        }
        
        segments(x0 = lolli_center, x1 = x_end, y0 = i, y1 = i, lwd = 1.2, col = delta_col)
        points(x = x_end, y = i, pch = 16, cex = 0.8, col = delta_col)
        text(
          x = x_end,
          y = i + 0.25,
          labels = sprintf("%.3f", delta_val),
          cex = 0.75,
          col = delta_col
        )
      }
    }
  }
  
  legend(
    "topleft",
    inset  = 0.02,
    bty    = "n",
    legend = c("Low edu", "Overall", "High edu"),
    fill   = c(llm_loedu, llm_overall, llm_hiedu)
  )
  
  invisible(NULL)
}


# ---- plotting helper: scatter of De vs PC1 by edu ------------------------- #

#' Scatter plot of effective dependence vs PC1 share, colored by education
#'
#' Produces a scatter plot (via ggplot2) of point-estimates per rater and
#' education group, plotting:
#'   - x = PC1 variance share
#'   - y = effective dependence (De)
#'
#' Raters are mapped to distinct single-character “glyph” shapes. Education
#' groups are mapped to colors; GSS is overplotted using an orange palette and
#' omitted from the color legend.
#'
#' @param dt A data.frame/data.table with columns \code{rater}, \code{edu_group},
#'   \code{pc1_var_explained}, \code{effective_dependence}.
#' @param out_file Optional output path. If \code{NULL}, prints to current device.
#' @param main Plot title.
#'
#' @return Invisibly, the ggplot object.
#'
#' @keywords internal
plot_scatter_De_vs_PC1_gg <- function(
    dt,
    out_file = NULL,
    main     = "Constraints on Generated Belief Systems"
) {
  stopifnot(is.data.frame(dt))
  
  req_cols <- c("rater", "edu_group", "pc1_var_explained", "effective_dependence")
  missing  <- setdiff(req_cols, names(dt))
  if (length(missing)) {
    stop("plot_scatter_De_vs_PC1_gg() missing columns: ", paste(missing, collapse = ", "))
  }
  
  if (!inherits(dt, "data.table")) data.table::setDT(dt)
  
  dt <- dt[is.finite(pc1_var_explained) & is.finite(effective_dependence)]
  if (!nrow(dt)) {
    warning("No finite rows to plot in plot_scatter_De_vs_PC1_gg().")
    return(invisible(NULL))
  }
  
  dt[is.na(edu_group) | edu_group == "", edu_group := "overall"]
  dt <- dt[edu_group %in% c("overall", "loedu", "hiedu")]
  dt[, edu_group := factor(edu_group, levels = c("overall", "loedu", "hiedu"))]
  
  rater_lvls <- sort(unique(dt$rater))
  if ("gss" %in% rater_lvls) rater_lvls <- c("gss", setdiff(rater_lvls, "gss"))
  dt[, rater := factor(rater, levels = rater_lvls)]
  
  make_shape_values <- function(levels_vec) {
    glyphs <- c(
      LETTERS[!LETTERS %in% c("O", "I")],
      letters[!letters %in% c("o", "i", "l")],
      as.character(0:9),
      c("+", "x", "*", "@", "#", "%", "&", "?", "!", "=", ":", ";", "~", "^")
    )
    glyphs <- unique(glyphs)
    
    if (length(levels_vec) > length(glyphs)) {
      warning("There are ", length(levels_vec), " raters but only ", length(glyphs),
              " unique glyph-shapes; shapes will repeat.")
      glyphs <- rep(glyphs, length.out = length(levels_vec))
    } else {
      glyphs <- glyphs[seq_len(length(levels_vec))]
    }
    
    names(glyphs) <- levels_vec
    glyphs
  }
  
  shape_vals <- make_shape_values(rater_lvls)
  
  edu_cols_llm <- c(overall = "dodgerblue", loedu = "grey70", hiedu = "dodgerblue4")
  
  non_gss <- dt
  gss_dt  <- dt[rater == "gss"]
  
  rater_labels_expr <- lapply(rater_lvls, function(x) {
    if (identical(x, "gss")) bquote(bold(.(toupper(x)))) else bquote(.(x))
  })
  rater_labels_expr <- do.call(expression, rater_labels_expr)
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data  = non_gss,
      mapping = ggplot2::aes(
        x      = pc1_var_explained,
        y      = effective_dependence,
        colour = edu_group,
        shape  = rater
      ),
      size  = 2.1,
      alpha = 0.85
    )
  
  if (nrow(gss_dt)) {
    p <- p +
      ggplot2::geom_point(
        data = gss_dt[edu_group == "overall"],
        ggplot2::aes(x = pc1_var_explained, y = effective_dependence, shape = rater),
        colour = "orange", size = 2.4, alpha = 0.95, show.legend = FALSE
      ) +
      ggplot2::geom_point(
        data = gss_dt[edu_group == "loedu"],
        ggplot2::aes(x = pc1_var_explained, y = effective_dependence, shape = rater),
        colour = "lightgoldenrod", size = 2.4, alpha = 0.95, show.legend = FALSE
      ) +
      ggplot2::geom_point(
        data = gss_dt[edu_group == "hiedu"],
        ggplot2::aes(x = pc1_var_explained, y = effective_dependence, shape = rater),
        colour = "orangered", size = 2.4, alpha = 0.95, show.legend = FALSE
      )
  }
  
  p <- p +
    ggplot2::scale_x_continuous(
      name   = "PC1 share of variance",
      limits = c(0, 1),
      expand = c(0, 0),
      breaks = scales::pretty_breaks(5),
      labels = scales::label_number(accuracy = 0.01)
    ) +
    ggplot2::scale_y_continuous(
      name   = "Effective dependence De (1 - |R|^(1/p))",
      limits = c(0, 1),
      expand = c(0, 0),
      breaks = scales::pretty_breaks(5),
      labels = scales::label_number(accuracy = 0.01)
    ) +
    ggplot2::scale_colour_manual(
      name   = "Education",
      values = edu_cols_llm,
      breaks = c("overall", "loedu", "hiedu"),
      labels = c("Overall", "Low edu", "High edu")
    ) +
    ggplot2::scale_shape_manual(
      name   = "Sample",
      values = shape_vals,
      breaks = rater_lvls,
      labels = rater_labels_expr
    ) +
    ggplot2::ggtitle(main) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      legend.position      = "right",
      legend.box           = "vertical",
      legend.title         = ggplot2::element_text(size = 9),
      legend.text          = ggplot2::element_text(size = 7),
      legend.spacing.y     = grid::unit(0.1, "lines"),
      legend.spacing.x     = grid::unit(0.2, "lines"),
      legend.key.height    = grid::unit(0.35 * 2, "lines"),
      legend.key.width     = grid::unit(0.8, "lines"),
      legend.margin        = ggplot2::margin(2, 2, 2, 2),
      legend.box.margin    = ggplot2::margin(2, 2, 2, 2),
      plot.title.position  = "plot",
      plot.title           = ggplot2::element_text(face = "bold"),
      panel.grid.minor     = ggplot2::element_blank()
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(order = 1, override.aes = list(alpha = 1, size = 2.2)),
      shape  = ggplot2::guide_legend(order = 2, override.aes = list(alpha = 1, size = 2.2))
    )
  
  if (!is.null(out_file)) {
    ggplot2::ggsave(filename = out_file, plot = p, width = 8, height = 3.5, units = "in")
  } else {
    print(p)
  }
  
  invisible(p)
}


# ---- RUN ------------------------------------------------------------------- #

# Pick up YEAR from global environment (default to 2024 if missing)
year <- get0("YEAR", ifnotfound = 2024L)

# ---- output directory for plots -------------------------------------------- #

viz_dir <- file.path(BASE_VIZ_DIR, sprintf("constraint_violins_%s", year))
dir.create(viz_dir, recursive = TRUE, showWarnings = FALSE)

# ---- print monitor for varset ---------------------------------------------- #
PRINT_VARSETS <- get0("PRINT_VARSETS", ifnotfound = FALSE)

# Discover which raters have directories of the form "<rater>-<year>"
raters <- available_raters(base_out_dir = BASE_OUT_DIR, year = year)
if (!length(raters)) {
  warning("No rater directories found in ", BASE_OUT_DIR, " for year ", year)
} else {
  message("Found raters: ", paste(raters, collapse = ", "))
}

# Map education labels -> filename suffixes used in polychor_bootstrap files
edu_suffixes <- c(overall = "", loedu = "loedu", hiedu = "hiedu")

# ---- compute union of variables present in at least 1 LLM (by edu group) --- #

llm_raters <- setdiff(raters, "gss")

llm_vars_by_group <- setNames(vector("list", length(edu_suffixes)), names(edu_suffixes))
for (grp in names(llm_vars_by_group)) llm_vars_by_group[[grp]] <- character(0)

if (length(llm_raters)) {
  for (r in llm_raters) {
    for (grp in names(edu_suffixes)) {
      suffix <- edu_suffixes[[grp]]
      
      # NOTE: we use strict = FALSE here so missing subgroup files simply skip
      corr_list <- load_corr_for_rater(
        rater                = r,
        base_out_dir         = BASE_OUT_DIR,
        year                 = year,
        suffix               = suffix,
        align_to_mode_varset = FALSE,
        strict               = FALSE
      )
      if (is.null(corr_list)) next
      
      vars_here <- unique(unlist(lapply(corr_list, colnames), use.names = FALSE))
      vars_here <- vars_here[!is.na(vars_here) & nzchar(vars_here)]
      
      llm_vars_by_group[[grp]] <- union(llm_vars_by_group[[grp]], vars_here)
    }
  }
}


for (grp in names(llm_vars_by_group)) {
  vars <- sort(unique(llm_vars_by_group[[grp]]))
  message("LLM var-union for ", grp, ": ", length(vars), " vars")
  if (isTRUE(PRINT_VARSETS)) {
    message("Vars used for ", grp, ":\n", paste0("  - ", vars, collapse = "\n"))
  }
}

# ---- compute bootstrap-level constraint stats ----------------------------- #

constraint_list <- vector("list", length(raters))
names(constraint_list) <- raters

for (r in raters) {
  per_rater <- vector("list", length(edu_suffixes))
  names(per_rater) <- names(edu_suffixes)
  
  for (grp in names(edu_suffixes)) {
    suffix <- edu_suffixes[[grp]]
    
    corr_list <- load_corr_for_rater(
      rater                = r,
      base_out_dir         = BASE_OUT_DIR,
      year                 = year,
      suffix               = suffix,
      align_to_mode_varset = FALSE,
      strict               = FALSE
    )
    if (is.null(corr_list)) next
    
    bootstrap_ids <- seq_along(corr_list)
    
    # For GSS only: restrict vars to those present in ≥1 LLM for this edu group
    if (identical(r, "gss")) {
      vars_keep <- llm_vars_by_group[[grp]]
      restricted    <- restrict_corr_list_to_vars(corr_list, vars_keep)
      corr_list     <- restricted$corr_list
      bootstrap_ids <- restricted$bootstrap_id
      if (is.null(corr_list) || !length(corr_list)) next
    }
    
    pc1_vals <- vapply(corr_list, pc1_var_explained,    numeric(1L))
    De_vals  <- vapply(corr_list, effective_dependence, numeric(1L))
    
    n_vars_vec <- vapply(
      corr_list,
      function(R) {
        R <- as.matrix(R)
        if (!is.numeric(R))               return(NA_integer_)
        if (nrow(R) < 2L || ncol(R) < 2L) return(NA_integer_)
        if (nrow(R) != ncol(R))           return(NA_integer_)
        ncol(R)
      },
      integer(1L)
    )
    
    per_rater[[grp]] <- data.table::data.table(
      rater                = r,
      edu_group            = grp,
      bootstrap_id         = bootstrap_ids,
      pc1_var_explained    = pc1_vals,
      effective_dependence = De_vals,
      n_vars               = n_vars_vec
    )
  }
  
  per_rater <- per_rater[!vapply(per_rater, is.null, logical(1L))]
  if (length(per_rater)) {
    constraint_list[[r]] <- data.table::rbindlist(per_rater, use.names = TRUE, fill = TRUE)
  }
}

non_empty_constraint <- constraint_list[!vapply(constraint_list, is.null, logical(1L))]
if (length(non_empty_constraint)) {
  constraint_dt <- data.table::rbindlist(non_empty_constraint, use.names = TRUE, fill = TRUE)
} else {
  constraint_dt <- data.table::data.table()
}


# ---- point-estimate summaries per rater × edu_group ----------------------- #

if (nrow(constraint_dt)) {
  constraint_point_dt <- constraint_dt[
    ,
    .(
      pc1_var_explained    = mean(pc1_var_explained,    na.rm = TRUE),
      effective_dependence = mean(effective_dependence, na.rm = TRUE),
      n_bootstraps         = .N,
      n_vars               = as.integer(stats::median(n_vars, na.rm = TRUE))
    ),
    by = .(rater, edu_group)
  ]
} else {
  constraint_point_dt <- data.table::data.table()
}


# ---- save tables ---------------------------------------------------------- #

if (nrow(constraint_dt)) {
  boot_csv <- file.path(BASE_OUT_DIR, sprintf("constraint_bootstrap_stats_%s.csv", year))
  data.table::fwrite(constraint_dt, boot_csv)
}

if (nrow(constraint_point_dt)) {
  point_csv <- file.path(BASE_OUT_DIR, sprintf("constraint_point_stats_%s.csv", year))
  data.table::fwrite(constraint_point_dt, point_csv)
}



# ---- violin plots --------------------------------------------------------- #

if (nrow(constraint_dt)) {
  plot_half_violins_edu(
    dt        = constraint_dt,
    value_col = "pc1_var_explained",
    xlab      = "Proportion of variance (attitudes + persona vars) explained by PC1",
    main      = "",
    out_file  = file.path(viz_dir, sprintf("constraint_pc1_var_explained_by_edu_%s.pdf", year))
  )
  
  plot_half_violins_edu(
    dt        = constraint_dt,
    value_col = "effective_dependence",
    xlab      = "Effective dependence De (1 - |R|^(1/p))",
    main      = "",
    lolli_pad_frac = 0.01,
    out_file  = file.path(viz_dir, sprintf("constraint_effective_dependence_by_edu_%s.pdf", year))
  )
}


# ---- scatter plot (point estimates) -------------------------------------- #

if (nrow(constraint_point_dt)) {
  scatter_file <- file.path(viz_dir, "constraint_pc1_vs_effective_dependence_by_edu.pdf")
  plot_scatter_De_vs_PC1_gg(dt = constraint_point_dt, out_file = scatter_file)
}

invisible(NULL)
