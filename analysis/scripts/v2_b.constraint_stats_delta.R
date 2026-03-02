################################################################################
# Persona conditioning diagnostics: beliefs vs beliefs|persona + variance split
#
# Overview
# --------
# This script quantifies how much “persona variables” constrain the correlation
# structure among beliefs, using bootstrap correlation matrices that include both
# belief variables and persona variables.
#
# Crucially, this script is written to:
#   - KEEP each bootstrap draw “as-is” (no alignment to a modal variable set),
#   - STILL enforce an “LLM-eligibility” restriction for non-LLM raters (e.g. GSS):
#       *A variable must appear in at least one LLM bootstrap draw* (within the
#       same education group) to be considered eligible for non-LLM raters.
#   - ALWAYS include persona variables in the eligible set, so conditioning remains
#     possible even under LLM-eligibility restriction.
#
# The analysis has two parts:
#
# Part A: paired constraint diagnostics (overall group; bootstrap distributions)
# ---------------------------------------------------------------------------
# For each rater and each bootstrap draw, we compute two scalar statistics on:
#   (1) R_BB : the belief-only correlation matrix (unconditional beliefs)
#   (2) R_B|P: the belief correlation matrix conditional on persona variables
#
# Conditional belief correlation is computed via the Schur complement:
#   Σ_{B|P} = R_BB - R_BP R_PP^{-1} R_PB
# then re-standardised into a correlation matrix:
#   R_{B|P} = D^{-1/2} Σ_{B|P} D^{-1/2}
#
# We produce paired half-violin plots per rater comparing:
#   - PC1 share of variance among beliefs: pc1_var_explained()
#   - Effective dependence among beliefs:  effective_dependence()
# for unconditional beliefs vs beliefs|persona.
#
# Part B: cumulative explained-variance decomposition (overall/loedu/hiedu)
# ------------------------------------------------------------------------
# For each rater × education group and each bootstrap draw, we decompose belief
# variance along belief PCs into persona-driven vs non-persona contributions:
#
#   - Compute persona-induced covariance among beliefs:
#       Σ_persona = R_BP R_PP^{-1} R_PB
#   - Eigendecompose belief correlation: R_BB = U Λ U'
#   - Persona variance along PC j:      u_j' Σ_persona u_j
#   - Non-persona variance along PC j:  λ_j - persona_j  
#   - Convert to shares of total belief variance tr(R_BB), then cumulate by PC.
#
# We produce two faceted PDFs (one facet per rater):
#   1) Cumulative persona vs non-persona shares (solid vs dashed), with education
#      overlaid by colour.
#   2) Cumulative total explained variance by belief PCs, with vertical markers
#      and labels for the PC index needed to reach 90% variance (per edu series).
#
# LLM-eligibility restriction (non-LLM raters only)
# ------------------------------------------------
# To avoid anchoring on survey-only items, non-LLM raters (default: "gss") are
# restricted per education group to:
#   eligible_vars(edu) = union_{LLM rater, draw} vars_in_draw(edu)
#                        UNION persona_vars
# This is applied by subsetting each non-LLM draw to intersect(colnames(R), eligible_vars).
# Draws are *not* otherwise modified or aligned, and are kept unless the subset
# leaves fewer than 2 variables (in which case the draw is dropped).
#
# Prerequisites (provided elsewhere in the project)
# ------------------------------------------------
# Required globals:
#   - BASE_OUT_DIR
#   - BASE_VIZ_DIR
#   - PERSONA_VARS_CANONICAL : character vector of persona variable names
#
# Required utils (assumed already sourced):
#   - available_raters(base_out_dir, year, exclude=...)
#   - load_corr_for_rater(rater, base_out_dir, year, suffix, align_to_mode_varset, strict)
#   - restrict_corr_list_to_vars(corr_list, vars_keep)
#   - pc1_var_explained(R)
#   - effective_dependence(R)
#
# Packages used
# -------------
#   - data.table
#   - ggplot2
#
# Output
# ------
# Writes PDFs into:
#   file.path(BASE_VIZ_DIR, sprintf("constraint_violins_%s", YEAR))
#
# Part A (overall only):
#   - pc1_beliefs_conditional_<YEAR>.pdf
#   - effective_dependence_beliefs_conditional_<YEAR>.pdf
#
# Part B (overall/loedu/hiedu):
#   - nonpersona_cumulative_pc_share_<YEAR>.pdf
#   - total_cumulative_pc_share_<YEAR>.pdf
################################################################################


#' Drop persona variables from a correlation matrix
#'
#' Given a correlation matrix \code{R} whose rows/columns correspond to a mixture
#' of belief variables and persona variables, this helper removes all persona
#' variables and returns the beliefs-only correlation matrix \eqn{R_{BB}}.
#'
#' This function is used for “unconditional beliefs” diagnostics: it describes
#' the correlation geometry among belief variables alone, without conditioning.
#'
#' @param R A square numeric correlation matrix. Must have column names.
#'   Row/column order is assumed to be consistent.
#' @param persona_vars Character vector of persona variable names to remove.
#'
#' @return A square numeric matrix containing only belief variables
#'   (\eqn{R_{BB}}), or \code{NULL} if:
#' \itemize{
#'   \item \code{R} is \code{NULL} or non-numeric,
#'   \item \code{R} has no column names (cannot match persona vars),
#'   \item fewer than 2 belief variables remain after dropping persona variables.
#' }
#'
#' @details
#' This function does not attempt to validate that \code{R} is a proper
#' correlation matrix (e.g., symmetric/PSD). It assumes upstream computation
#' produced a valid correlation matrix.
#'
#' @keywords internal
drop_persona_vars <- function(R, persona_vars = PERSONA_VARS_CANONICAL) {
  if (is.null(R)) return(NULL)
  R <- as.matrix(R)
  if (!is.numeric(R)) return(NULL)
  
  cn <- colnames(R)
  if (is.null(cn)) stop("drop_persona_vars(): R must have column names.")
  if (is.null(persona_vars) || !length(persona_vars)) {
    stop("drop_persona_vars(): persona_vars not found or empty.")
  }
  
  keep <- !(cn %in% persona_vars)
  if (sum(keep) < 2L) return(NULL)
  R[keep, keep, drop = FALSE]
}


#' Conditional belief correlation given persona variables
#'
#' Computes the beliefs-only conditional correlation matrix \eqn{R_{B|P}} from a
#' correlation matrix over \emph{(belief + persona)} variables.
#'
#' Let \eqn{B} denote belief variables (all non-persona variables) and \eqn{P}
#' denote persona variables (those listed in \code{persona_vars}). Consider the
#' corresponding block partition of \eqn{R}:
#' \deqn{
#' R =
#' \begin{bmatrix}
#' R_{BB} & R_{BP} \\
#' R_{PB} & R_{PP}
#' \end{bmatrix}.
#' }
#'
#' The conditional covariance of beliefs given persona is the Schur complement:
#' \deqn{\Sigma_{B|P} = R_{BB} - R_{BP}R_{PP}^{-1}R_{PB}.}
#' Since \eqn{\Sigma_{B|P}} is a covariance matrix (not correlation), we
#' re-standardise it into a correlation matrix using its diagonal:
#' \deqn{
#' R_{B|P} = D^{-1/2} \Sigma_{B|P} D^{-1/2},
#' \quad D = \mathrm{diag}(\Sigma_{B|P}).
#' }
#'
#' @param R A square numeric correlation matrix over beliefs + persona variables.
#'   Must have column names.
#' @param persona_vars Character vector of persona variable names to condition on.
#'
#' @return A square numeric beliefs-only correlation matrix \eqn{R_{B|P}}, or
#'   \code{NULL} if:
#' \itemize{
#'   \item none of \code{persona_vars} are found in \code{colnames(R)},
#'   \item fewer than 2 belief variables remain after excluding persona vars,
#'   \item \eqn{R_{PP}} inversion fails (singular/ill-conditioned),
#'   \item conditional variances are non-finite or non-positive (cannot standardise).
#' }
#'
#' @details
#' Numerical safeguards:
#' \itemize{
#'   \item Output is symmetrised via \code{(R + t(R))/2} to damp tiny asymmetries.
#'   \item Diagonal is forced to 1 after standardisation.
#' }
#'
#' Conceptual note:
#' Conditioning removes the linear dependence among beliefs that can be explained
#' by persona variables. Comparing diagnostics on \eqn{R_{BB}} vs \eqn{R_{B|P}}
#' helps quantify how much persona structure “drives” apparent constraint among beliefs.
#'
#' @keywords internal
cond_belief_corr_on_persona <- function(R, persona_vars = PERSONA_VARS_CANONICAL) {
  if (is.null(R)) return(NULL)
  R <- as.matrix(R)
  if (!is.numeric(R)) stop("cond_belief_corr_on_persona(): R must be numeric.")
  
  cn <- colnames(R)
  if (is.null(cn)) stop("cond_belief_corr_on_persona(): R must have column names.")
  if (is.null(persona_vars) || !length(persona_vars)) {
    stop("cond_belief_corr_on_persona(): persona_vars not found or empty.")
  }
  
  idxP <- match(persona_vars, cn)
  idxP <- idxP[!is.na(idxP)]
  if (!length(idxP)) return(NULL)
  
  idxB <- setdiff(seq_len(ncol(R)), idxP)
  if (length(idxB) < 2L) return(NULL)
  
  R_BB <- R[idxB, idxB, drop = FALSE]
  R_BP <- R[idxB, idxP, drop = FALSE]
  R_PP <- R[idxP, idxP, drop = FALSE]
  
  inv_R_PP <- try(solve(R_PP), silent = TRUE)
  if (inherits(inv_R_PP, "try-error")) return(NULL)
  
  Sigma_B_given_P <- R_BB - R_BP %*% inv_R_PP %*% t(R_BP)
  
  v <- diag(Sigma_B_given_P)
  if (any(!is.finite(v)) || any(v <= 0)) return(NULL)
  
  D_inv_sqrt  <- diag(1 / sqrt(v))
  R_B_given_P <- D_inv_sqrt %*% Sigma_B_given_P %*% D_inv_sqrt
  
  R_B_given_P <- (R_B_given_P + t(R_B_given_P)) / 2
  diag(R_B_given_P) <- 1
  R_B_given_P
}


#' Plot paired half-violin distributions for two conditions per rater
#'
#' Writes a PDF that compares two bootstrap distributions of a scalar statistic,
#' for each rater. Each rater is displayed on its own horizontal row and
#' contains:
#' \itemize{
#'   \item a half-violin for the “full” distribution (\code{full_col}),
#'   \item a half-violin for the “reduced” distribution (\code{reduced_col}),
#'   \item a median marker for each distribution,
#'   \item a right-margin annotation showing
#'     \eqn{\Delta = full - reduced} with a 95% bootstrap interval.
#' }
#'
#' Raters are ordered by decreasing median \eqn{\Delta} so that raters with the
#' largest persona-conditioning effect appear at the top.
#'
#' @param dt Data frame/data.table with one row per bootstrap draw per rater,
#'   containing at least:
#'   \itemize{
#'     \item \code{rater} (character)
#'     \item \code{full_col} (numeric)
#'     \item \code{reduced_col} (numeric)
#'   }
#' @param full_col Column name for the “full/unconditional” statistic.
#' @param reduced_col Column name for the “reduced/conditional” statistic.
#' @param xlab X-axis label.
#' @param main Title displayed at top of plot.
#' @param out_file Output PDF path.
#' @param full_label,reduced_label Legend labels for the two conditions.
#' @param full_color,reduced_color Fill/median colours for the two conditions.
#' @param legend_pos Legend position (passed to \code{legend()}).
#'
#' @return Invisibly \code{NULL}. Called for its side effect (writes a PDF).
#'
#' @details
#' Density estimation uses \code{stats::density()} with default bandwidth
#' selection. Distributions with fewer than 2 distinct finite values are skipped.
#'
#' The delta interval is computed from the bootstrap distribution of
#' \eqn{\Delta_i = full_i - reduced_i} within each rater:
#' \itemize{
#'   \item median: \code{median(delta)}
#'   \item interval: \code{quantile(delta, c(.025, .975))}
#' }
#'
#' @keywords internal
plot_half_violins_two_conditions <- function(
    dt,
    full_col,
    reduced_col,
    xlab,
    main,
    out_file,
    full_label,
    reduced_label,
    full_color    = "#A6CEE3",
    reduced_color = "#FDB863",
    legend_pos    = "topright",
    lolli_pad_frac   = -0.01,
    lolli_width_frac = 0.14
) {
  stopifnot(is.data.frame(dt))
  if (!("rater" %in% names(dt))) stop("plot_half_violins_two_conditions() requires a 'rater' column.")
  if (!(full_col %in% names(dt))) stop("Column ", sQuote(full_col), " not found.")
  if (!(reduced_col %in% names(dt))) stop("Column ", sQuote(reduced_col), " not found.")
  
  if (!inherits(dt, "data.table")) data.table::setDT(dt)
  if (!nrow(dt)) {
    warning("Empty data in plot_half_violins_two_conditions(); skipping plot.")
    return(invisible(NULL))
  }
  
  delta_dt <- dt[, {
    full_vals    <- as.numeric(get(full_col))
    reduced_vals <- as.numeric(get(reduced_col))
    delta <- full_vals - reduced_vals
    delta <- delta[is.finite(delta)]
    if (!length(delta)) {
      .(delta_median = NA_real_, delta_lower = NA_real_, delta_upper = NA_real_)
    } else {
      .(
        delta_median = stats::median(delta, na.rm = TRUE),
        delta_lower  = stats::quantile(delta, probs = 0.025, na.rm = TRUE),
        delta_upper  = stats::quantile(delta, probs = 0.975, na.rm = TRUE)
      )
    }
  }, by = rater]
  
  delta_dt <- delta_dt[order(-delta_median, na.last = TRUE)]
  raters_order <- delta_dt$rater
  
  dt[, rater := factor(rater, levels = raters_order)]
  n_raters <- length(raters_order)
  
  all_vals <- unlist(list(as.numeric(dt[[full_col]]), as.numeric(dt[[reduced_col]])), use.names = FALSE)
  all_vals <- all_vals[is.finite(all_vals)]
  if (!length(all_vals)) {
    warning("No finite values in ", full_col, " or ", reduced_col, "; skipping plot.")
    return(invisible(NULL))
  }
  
  x_range <- range(all_vals, finite = TRUE)
  if (!is.finite(diff(x_range)) || diff(x_range) == 0) {
    x_range <- x_range + c(-0.5, 0.5)
  } else {
    pad     <- 0.05 * diff(x_range)
    x_range <- x_range + c(-pad, pad)
  }
  
  pdf_width  <- 11
  pdf_height <- max(3, 0.3 * n_raters + 1)
  
  grDevices::pdf(out_file, width = pdf_width, height = pdf_height)
  on.exit(grDevices::dev.off(), add = TRUE)
  
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  
  par(mar = c(4, 18, 4, 8), xpd = NA)
  
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

  # Draw vertical gridlines ONLY inside the plot panel (no bleed into margins)
  grid_x <- pretty(x_range)
  
  op_xpd <- par(xpd = FALSE)  # clip to plot region
  usr <- par("usr")
  segments(x0 = grid_x, y0 = usr[3], x1 = grid_x, y1 = usr[4],
           col = "grey92", lwd = 0.45)
  par(op_xpd)
  
  # Restore margin drawing for right-side delta labels
  par(xpd = NA)
  
  
  max_width <- 0.4
  usr <- par("usr")
  
  # Lollipop layout (to the right of the plot area)
  x_right <- usr[2]
  lolli_pad    <- lolli_pad_frac * diff(usr[1:2])
  lolli_width  <- lolli_width_frac * diff(usr[1:2])
  lolli_center <- x_right + lolli_pad + (lolli_width / 2)
  
  max_abs_delta <- max(abs(delta_dt$delta_median), na.rm = TRUE)
  if (!is.finite(max_abs_delta) || max_abs_delta == 0) max_abs_delta <- 1
  
  draw_half <- function(vals, base_y, fill_col) {
    vals <- vals[is.finite(vals)]
    if (length(vals) < 2L || length(unique(vals)) < 2L) return()
    dens <- stats::density(vals)
    scaled_y <- dens$y / max(dens$y, na.rm = TRUE) * max_width
    x_poly   <- c(dens$x, rev(dens$x))
    y_poly   <- c(rep(base_y, length(dens$x)), base_y + rev(scaled_y))
    polygon(x_poly, y_poly, col = grDevices::adjustcolor(fill_col, alpha.f = 0.7), border = NA)
    lines(dens$x, base_y + scaled_y, col = grDevices::adjustcolor("black", alpha.f = 0.5))
  }
  
  for (i in seq_along(raters_order)) {
    r <- raters_order[i]
    dt_r <- dt[rater == r]
    
    full_vals    <- as.numeric(dt_r[[full_col]])
    reduced_vals <- as.numeric(dt_r[[reduced_col]])
    
    draw_half(reduced_vals, i, reduced_color)
    if (any(is.finite(reduced_vals))) {
      med_reduced <- stats::median(reduced_vals, na.rm = TRUE)
      segments(med_reduced, i, med_reduced, i + max_width * 1.05, lwd = 1.8, col = reduced_color)
    }
    
    draw_half(full_vals, i, full_color)
    if (any(is.finite(full_vals))) {
      med_full <- stats::median(full_vals, na.rm = TRUE)
      segments(med_full, i, med_full, i + max_width * 1.05, lwd = 1.8, col = full_color)
    }
    
    drow <- delta_dt[rater == r]
    if (nrow(drow) == 1L && is.finite(drow$delta_median)) {
      delta_val <- drow$delta_median
      x_end <- lolli_center + (delta_val / max_abs_delta) * (lolli_width / 2)
      
      sig_delta <- is.finite(drow$delta_lower) && is.finite(drow$delta_upper) &&
        (drow$delta_lower > 0 || drow$delta_upper < 0)
      delta_col <- if (sig_delta) "black" else "grey60"
      
      if (is.finite(drow$delta_lower) && is.finite(drow$delta_upper)) {
        x_lwr <- lolli_center + (drow$delta_lower / max_abs_delta) * (lolli_width / 2)
        x_upr <- lolli_center + (drow$delta_upper / max_abs_delta) * (lolli_width / 2)
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
  
  legend(
    legend_pos,
    inset  = 0.02,
    bty    = "n",
    legend = c(full_label, reduced_label),
    fill   = c(
      grDevices::adjustcolor(full_color,    alpha.f = 0.7),
      grDevices::adjustcolor(reduced_color, alpha.f = 0.7)
    )
  )
  
  invisible(NULL)
}


#' Bootstrap summary of persona vs non-persona cumulative explained variance by PC
#'
#' For a given rater and education group, this function loads bootstrap
#' correlation matrices over (belief + persona) variables and computes, per draw,
#' how much explained belief variance is attributable to persona vs non-persona
#' sources along belief principal components (PCs).
#'
#' The computation per bootstrap draw is:
#' \enumerate{
#'   \item Partition variables into beliefs (B) and persona (P).
#'   \item Form blocks \eqn{R_{BB}}, \eqn{R_{BP}}, \eqn{R_{PP}}.
#'   \item Compute persona-induced covariance among beliefs:
#'         \eqn{\Sigma_{persona} = R_{BP} R_{PP}^{-1} R_{PB}}.
#'   \item Eigendecompose belief correlation \eqn{R_{BB} = U \Lambda U'}.
#'   \item For each PC j:
#'         \itemize{
#'           \item total variance: \eqn{\lambda_j}
#'           \item persona variance: \eqn{u_j' \Sigma_{persona} u_j}
#'           \item non-persona variance: \eqn{\max(\lambda_j - persona_j, 0)}
#'         }
#'   \item Convert to shares of total belief variance \eqn{tr(R_{BB})}, then
#'         take cumulative sums across PCs.
#' }
#'
#' This script intentionally does *not* align variable sets across bootstraps.
#' Each draw is used with whatever variable set it contains, after applying an
#' optional restriction set \code{vars_keep}. Because draws may have different
#' numbers of PCs, cumulative vectors are padded with \code{NA} for aggregation.
#'
#' @param rater Character rater id.
#' @param base_out_dir Base directory containing "<rater>-<year>" subdirectories.
#' @param year Integer year.
#' @param persona_vars Character vector of persona variable names.
#' @param suffix One of \code{""}, \code{"loedu"}, \code{"hiedu"} selecting which
#'   bootstrap file to load for the rater.
#' @param edu_group Label stored in output (\code{"overall"}, \code{"loedu"}, \code{"hiedu"}).
#' @param max_pcs Maximum number of PCs to retain (default \code{Inf}).
#' @param vars_keep Optional character vector restricting matrices before analysis.
#'   Typically used for non-LLM raters to enforce “var appears in ≥1 LLM” eligibility.
#'
#' @return A \code{data.table} with one row per \code{pc_index} and columns:
#' \itemize{
#'   \item \code{rater}, \code{edu_group}, \code{pc_index}
#'   \item \code{cum_non_*}: cumulative non-persona share (median/low/high)
#'   \item \code{cum_persona_*}: cumulative persona share (median/low/high)
#'   \item \code{cum_total_*}: cumulative total share (median/low/high)
#'   \item \code{share_total_*}: per-PC total share (median/low/high), computed
#'     from successive differences of cumulative totals
#' }
#'
#' @details
#' Aggregation across bootstraps:
#' \itemize{
#'   \item Compute cumulative vectors per draw.
#'   \item Let \code{max_len} be the maximum vector length across draws (capped by
#'     \code{max_pcs} if finite).
#'   \item Pad shorter vectors with \code{NA} up to \code{max_len}.
#'   \item Take column-wise medians and 2.5%/97.5% quantiles with \code{na.rm=TRUE}.
#' }
#'
#' Numerical considerations:
#' \itemize{
#'   \item If \eqn{R_{PP}} inversion fails, the draw is skipped.
#' }
#'
#' @keywords internal
nonpersona_cumvar_for_rater <- function(
    rater,
    base_out_dir = BASE_OUT_DIR,
    year         = get0("YEAR", ifnotfound = 2024L),
    persona_vars = PERSONA_VARS_CANONICAL,
    suffix       = "",
    edu_group    = "overall",
    max_pcs      = Inf,
    vars_keep    = NULL
) {
  corr_list <- load_corr_for_rater(
    rater                = rater,
    base_out_dir         = base_out_dir,
    year                 = year,
    suffix               = suffix,
    align_to_mode_varset = FALSE,  # KEEP each draw as-is
    strict               = FALSE
  )
  if (is.null(corr_list) || !length(corr_list)) return(NULL)
  
  if (!is.null(vars_keep) && length(vars_keep)) {
    restricted <- restrict_corr_list_to_vars(corr_list, vars_keep = vars_keep)
    corr_list <- restricted$corr_list
    if (is.null(corr_list) || !length(corr_list)) return(NULL)
  }
  
  cum_non_list <- vector("list", length(corr_list))
  cum_tot_list <- vector("list", length(corr_list))
  
  for (i in seq_along(corr_list)) {
    R_full <- corr_list[[i]]
    if (is.null(R_full)) { cum_non_list[[i]] <- NA_real_; cum_tot_list[[i]] <- NA_real_; next }
    
    R_full <- as.matrix(R_full)
    if (!is.numeric(R_full)) { cum_non_list[[i]] <- NA_real_; cum_tot_list[[i]] <- NA_real_; next }
    
    cn <- colnames(R_full)
    if (is.null(cn)) { cum_non_list[[i]] <- NA_real_; cum_tot_list[[i]] <- NA_real_; next }
    
    idxP <- match(persona_vars, cn)
    idxP <- idxP[!is.na(idxP)]
    if (!length(idxP)) { cum_non_list[[i]] <- NA_real_; cum_tot_list[[i]] <- NA_real_; next }
    
    idxB <- setdiff(seq_len(ncol(R_full)), idxP)
    if (length(idxB) < 2L) { cum_non_list[[i]] <- NA_real_; cum_tot_list[[i]] <- NA_real_; next }
    
    R_BB <- R_full[idxB, idxB, drop = FALSE]
    R_BP <- R_full[idxB, idxP, drop = FALSE]
    R_PP <- R_full[idxP, idxP, drop = FALSE]
    
    inv_R_PP <- try(solve(R_PP), silent = TRUE)
    if (inherits(inv_R_PP, "try-error")) { cum_non_list[[i]] <- NA_real_; cum_tot_list[[i]] <- NA_real_; next }
    
    Sigma_persona <- R_BP %*% inv_R_PP %*% t(R_BP)
    
    ev_R <- try(eigen(R_BB, symmetric = TRUE), silent = TRUE)
    if (inherits(ev_R, "try-error")) { cum_non_list[[i]] <- NA_real_; cum_tot_list[[i]] <- NA_real_; next }
    
    lambda <- Re(ev_R$values)
    U      <- ev_R$vectors
    ok_eig <- is.finite(lambda)
    if (!any(ok_eig)) { cum_non_list[[i]] <- NA_real_; cum_tot_list[[i]] <- NA_real_; next }
    
    lambda <- lambda[ok_eig]
    U      <- U[, ok_eig, drop = FALSE]
    
    total_var <- sum(diag(R_BB))
    if (!is.finite(total_var) || total_var <= 0) { cum_non_list[[i]] <- NA_real_; cum_tot_list[[i]] <- NA_real_; next }
    
    Sigma_pc       <- crossprod(U, Sigma_persona %*% U)
    var_persona_pc <- Re(diag(Sigma_pc))
    
    var_nonpersona_pc   <- lambda - var_persona_pc
    
    share_nonpersona_pc <- var_nonpersona_pc / total_var
    share_total_pc      <- lambda / total_var
    
    cum_non <- cumsum(share_nonpersona_pc)
    cum_tot <- cumsum(share_total_pc)
    
    if (is.finite(max_pcs)) {
      m <- min(length(cum_non), max_pcs)
      cum_non <- cum_non[seq_len(m)]
      cum_tot <- cum_tot[seq_len(m)]
    }
    
    cum_non_list[[i]] <- cum_non
    cum_tot_list[[i]] <- cum_tot
  }
  
  ok <- vapply(cum_non_list, function(x) is.numeric(x) && any(is.finite(x)), logical(1L))
  if (!any(ok)) return(NULL)
  cum_non_list <- cum_non_list[ok]
  cum_tot_list <- cum_tot_list[ok]
  
  len_vec  <- vapply(cum_non_list, length, integer(1L))
  max_len0 <- max(len_vec)
  max_len  <- if (is.finite(max_pcs)) min(max_len0, max_pcs) else max_len0
  
  n_draws <- length(cum_non_list)
  cum_non_mat <- matrix(NA_real_, nrow = n_draws, ncol = max_len)
  cum_tot_mat <- matrix(NA_real_, nrow = n_draws, ncol = max_len)
  
  for (j in seq_len(n_draws)) {
    lj <- min(length(cum_non_list[[j]]), max_len)
    if (lj > 0L) {
      cum_non_mat[j, 1:lj] <- cum_non_list[[j]][1:lj]
      cum_tot_mat[j, 1:lj] <- cum_tot_list[[j]][1:lj]
    }
  }
  
  col_med  <- function(m) apply(m, 2L, stats::median,   na.rm = TRUE)
  col_low  <- function(m) apply(m, 2L, stats::quantile, probs = 0.025, na.rm = TRUE)
  col_high <- function(m) apply(m, 2L, stats::quantile, probs = 0.975, na.rm = TRUE)
  
  med_non   <- col_med(cum_non_mat); low_non <- col_low(cum_non_mat); high_non <- col_high(cum_non_mat)
  med_tot   <- col_med(cum_tot_mat); low_tot <- col_low(cum_tot_mat); high_tot <- col_high(cum_tot_mat)
  
  cum_persona_mat <- cum_tot_mat - cum_non_mat
  med_per   <- col_med(cum_persona_mat); low_per <- col_low(cum_persona_mat); high_per <- col_high(cum_persona_mat)
  
  diff_mat <- cum_non_mat - cum_persona_mat
  med_diff <- col_med(diff_mat); low_diff <- col_low(diff_mat); high_diff <- col_high(diff_mat)
  
  share_tot_mat <- cbind(
    cum_tot_mat[, 1, drop = FALSE],
    cum_tot_mat[, -1, drop = FALSE] - cum_tot_mat[, -max_len, drop = FALSE]
  )
  share_tot_med  <- col_med(share_tot_mat)
  share_tot_low  <- col_low(share_tot_mat)
  share_tot_high <- col_high(share_tot_mat)
  
  data.table::data.table(
    rater            = rater,
    edu_group        = edu_group,
    pc_index         = seq_len(max_len),
    cum_non_med      = med_non,
    cum_non_low      = low_non,
    cum_non_high     = high_non,
    cum_total_med    = med_tot,
    cum_total_low    = low_tot,
    cum_total_high   = high_tot,
    cum_persona_med  = med_per,
    cum_persona_low  = low_per,
    cum_persona_high = high_per,
    diff_med         = med_diff,
    diff_low         = low_diff,
    diff_high        = high_diff,
    share_total_med  = share_tot_med,
    share_total_low  = share_tot_low,
    share_total_high = share_tot_high
  )
}


# ---- RUN ------------------------------------------------------------------- #

# Pick up YEAR from global environment (default to 2024 if missing)
year <- get0("YEAR", ifnotfound = 2024L)

# ---- output directory for plots -------------------------------------------- #
viz_dir <- file.path(BASE_VIZ_DIR, sprintf("constraint_violins_%s", year))
dir.create(viz_dir, recursive = TRUE, showWarnings = FALSE)

# ---- make sure that persona vars canon is defined -------------------------- #
persona_vars_global <- get0("PERSONA_VARS_CANONICAL", ifnotfound = NULL)
stopifnot(!is.null(persona_vars_global), length(persona_vars_global) > 0)

# Discover which raters have directories of the form "<rater>-<year>"
raters <- available_raters(base_out_dir = BASE_OUT_DIR, year = year)
if (!length(raters)) stop("No rater directories found in ", BASE_OUT_DIR, " for year ", year)

non_llm_raters <- c("gss")
llm_raters     <- setdiff(raters, non_llm_raters)

edu_suffixes <- c(overall = "", loedu = "loedu", hiedu = "hiedu")

# Compute var-union across LLMs per education group WITHOUT aligning draws.
llm_vars_by_group <- setNames(vector("list", length(edu_suffixes)), names(edu_suffixes))
for (grp in names(llm_vars_by_group)) llm_vars_by_group[[grp]] <- character(0)

if (length(llm_raters)) {
  for (r in llm_raters) {
    for (grp in names(edu_suffixes)) {
      suffix <- edu_suffixes[[grp]]
      
      corr_list <- load_corr_for_rater(
        rater                = r,
        base_out_dir         = BASE_OUT_DIR,
        year                 = year,
        suffix               = suffix,
        align_to_mode_varset = FALSE,  # keep each draw as-is
        strict               = FALSE
      )
      if (is.null(corr_list) || !length(corr_list)) next
      
      vars_here <- unique(unlist(lapply(corr_list, colnames), use.names = FALSE))
      vars_here <- vars_here[!is.na(vars_here) & nzchar(vars_here)]
      llm_vars_by_group[[grp]] <- union(llm_vars_by_group[[grp]], vars_here)
    }
  }
}



# Always include persona vars so conditioning remains possible under restriction.
llm_vars_by_group_raw <- llm_vars_by_group
PRINT_VARSETS <- get0("PRINT_VARSETS", ifnotfound = FALSE)

for (grp in names(llm_vars_by_group)) {
  raw   <- sort(unique(llm_vars_by_group_raw[[grp]]))
  final <- sort(unique(union(raw, persona_vars_global)))
  added <- sort(setdiff(final, raw))
  
  llm_vars_by_group[[grp]] <- final
  
  message("LLM var-union (observed) for ", grp, ": ", length(raw), " vars")
  message("LLM-eligible vars for ", grp, ": ", length(final), " (incl persona)")
  message("Vars added from persona canon for ", grp, ": ", length(added))
  
  if (isTRUE(PRINT_VARSETS)) {
    message("Observed union vars for ", grp, ":\n", paste0("  - ", raw, collapse = "\n"))
    if (length(added)) {
      message("Added from persona canon for ", grp, ":\n", paste0("  - ", added, collapse = "\n"))
    } else {
      message("Added from persona canon for ", grp, ": (none)")
    }
    message("Final eligible vars for ", grp, ":\n", paste0("  - ", final, collapse = "\n"))
  }
}


## =============================================================================
## Part A: beliefs vs beliefs|persona (overall only; keep each draw as-is)
## =============================================================================

constraint_list <- vector("list", length(raters))
names(constraint_list) <- raters

for (r in raters) {
  corr_list <- load_corr_for_rater(
    rater                = r,
    base_out_dir         = BASE_OUT_DIR,
    year                 = year,
    suffix               = "",
    align_to_mode_varset = FALSE,  # IMPORTANT: no modal alignment
    strict               = FALSE
  )
  if (is.null(corr_list) || !length(corr_list)) next
  
  bootstrap_ids <- seq_along(corr_list)
  
  # Restrict non-LLM raters to LLM-eligible vars (overall)
  if (r %in% non_llm_raters) {
    restricted <- restrict_corr_list_to_vars(corr_list, vars_keep = llm_vars_by_group[["overall"]])
    corr_list     <- restricted$corr_list
    bootstrap_ids <- restricted$bootstrap_id
    if (is.null(corr_list) || !length(corr_list)) next
  }
  
  n_draws <- length(corr_list)
  
  pc1_beliefs      <- rep(NA_real_, n_draws)
  pc1_beliefs_cond <- rep(NA_real_, n_draws)
  De_beliefs       <- rep(NA_real_, n_draws)
  De_beliefs_cond  <- rep(NA_real_, n_draws)
  
  for (i in seq_along(corr_list)) {
    R_full <- corr_list[[i]]
    
    R_beliefs <- drop_persona_vars(R_full, persona_vars = persona_vars_global)
    if (!is.null(R_beliefs)) {
      pc1_beliefs[i] <- pc1_var_explained(R_beliefs)
      De_beliefs[i]  <- effective_dependence(R_beliefs)
    }
    
    R_cond <- cond_belief_corr_on_persona(R_full, persona_vars = persona_vars_global)
    if (!is.null(R_cond)) {
      pc1_beliefs_cond[i] <- pc1_var_explained(R_cond)
      De_beliefs_cond[i]  <- effective_dependence(R_cond)
    }
  }
  
  constraint_list[[r]] <- data.table::data.table(
    rater            = r,
    bootstrap_id     = bootstrap_ids,
    pc1_beliefs      = pc1_beliefs,
    pc1_beliefs_cond = pc1_beliefs_cond,
    De_beliefs       = De_beliefs,
    De_beliefs_cond  = De_beliefs_cond
  )
}

constraint_dt <- data.table::rbindlist(
  constraint_list[!vapply(constraint_list, is.null, logical(1L))],
  use.names = TRUE, fill = TRUE
)

if (nrow(constraint_dt) && !isTRUE(get0("SKIP_V2B_PLOTS", ifnotfound = FALSE))) {
  plot_half_violins_two_conditions(
    dt            = constraint_dt,
    full_col      = "pc1_beliefs",
    reduced_col   = "pc1_beliefs_cond",
    xlab          = "PC1 share of variance among beliefs",
    main          = "",
    out_file      = file.path(viz_dir, sprintf("pc1_beliefs_conditional_%s.pdf", year)),
    full_label    = "Beliefs (unconditional)",
    reduced_label = "Beliefs | persona",
    legend_pos    = "topright"
  )
  
  plot_half_violins_two_conditions(
    dt            = constraint_dt,
    full_col      = "De_beliefs",
    reduced_col   = "De_beliefs_cond",
    xlab          = "Effective dependence De among beliefs",
    main          = "Effective dependence among beliefs: unconditional vs conditional on persona",
    out_file      = file.path(viz_dir, sprintf("effective_dependence_beliefs_conditional_%s.pdf", year)),
    full_label    = "Beliefs (unconditional)",
    reduced_label = "Beliefs | persona",
    legend_pos    = "topleft"
  )
} else {
  warning("No constraint statistics computed in Part A; skipping Part A plots.")
}


## =============================================================================
## Part B: cumulative variance decomposition by PC (edu overlay; keep draws as-is)
## =============================================================================

nonpersona_cumvar_list <- list()

for (r in raters) {
  for (grp in names(edu_suffixes)) {
    suf <- edu_suffixes[[grp]]
    vars_keep <- if (r %in% non_llm_raters) llm_vars_by_group[[grp]] else NULL
    
    res <- nonpersona_cumvar_for_rater(
      rater        = r,
      base_out_dir = BASE_OUT_DIR,
      year         = year,
      persona_vars = persona_vars_global,
      suffix       = suf,
      edu_group    = grp,
      max_pcs      = Inf,
      vars_keep    = vars_keep
    )
    if (!is.null(res) && nrow(res)) {
      nonpersona_cumvar_list[[paste(r, grp, sep = "_")]] <- res
    }
  }
}

nonpersona_cumvar_dt <- data.table::rbindlist(nonpersona_cumvar_list, use.names = TRUE, fill = TRUE)
if (!nrow(nonpersona_cumvar_dt)) stop("No non-persona cumulative-variance statistics computed; nothing to plot.")

# Identify rater type and order raters by PC1 total share (overall if available)
rater_gss <- "gss"
nonpersona_cumvar_dt[, is_gss     := as.character(rater) == rater_gss]
nonpersona_cumvar_dt[, rater_type := ifelse(is_gss, "GSS", "LLM")]

# Order by PC1 share of total variance (prefer overall; fall back to any PC1)
order_dt <- nonpersona_cumvar_dt[pc_index == 1 & edu_group == "overall",
                                 .(pc1_share = share_total_med),
                                 by = rater]
if (!nrow(order_dt)) {
  order_dt <- nonpersona_cumvar_dt[pc_index == 1,
                                   .(pc1_share = share_total_med),
                                   by = rater]
}
order_dt <- order_dt[order(-pc1_share)]
ordered_raters <- as.character(order_dt$rater)

# Ensure all raters appear (avoid NA facets if some raters lack PC1 rows)
all_raters <- unique(as.character(nonpersona_cumvar_dt$rater))
ordered_raters <- c(ordered_raters, setdiff(all_raters, ordered_raters))

nonpersona_cumvar_dt[, rater     := factor(as.character(rater), levels = ordered_raters)]
nonpersona_cumvar_dt[, edu_group := factor(edu_group, levels = c("overall", "loedu", "hiedu"))]

# Colour palette consistent with earlier constraint plots
edu_palette_df <- data.table::data.table(
  rater_type = rep(c("LLM", "GSS"), each = 3L),
  edu_group  = rep(c("overall", "loedu", "hiedu"), times = 2L)
)

edu_palette_df[
  rater_type == "LLM" & edu_group == "overall", colour := "dodgerblue"
][rater_type == "LLM" & edu_group == "loedu",   colour := "grey70"
][rater_type == "LLM" & edu_group == "hiedu",   colour := "dodgerblue4"
][rater_type == "GSS" & edu_group == "overall", colour := "orange"
][rater_type == "GSS" & edu_group == "loedu",   colour := "lightgoldenrod"
][rater_type == "GSS" & edu_group == "hiedu",   colour := "orangered"
]

edu_palette_df[, edu_series := paste(rater_type, edu_group, sep = ": ")]
edu_colours <- edu_palette_df$colour
names(edu_colours) <- edu_palette_df$edu_series

nonpersona_cumvar_dt <- merge(
  nonpersona_cumvar_dt,
  edu_palette_df[, .(rater_type, edu_group, edu_series)],
  by = c("rater_type", "edu_group"),
  all.x = TRUE
)

# Re-impose facet ordering (merge may drop factor levels)
nonpersona_cumvar_dt[, rater := factor(as.character(rater), levels = ordered_raters)]

## ---- Plot 1: cumulative persona vs non-persona shares -------------------- ##

df_comp <- data.table::rbindlist(list(
  nonpersona_cumvar_dt[, .(
    rater, edu_group, rater_type, edu_series, pc_index,
    cum_low = cum_non_low, cum_med = cum_non_med, cum_high = cum_non_high,
    component = "Non-persona"
  )],
  nonpersona_cumvar_dt[, .(
    rater, edu_group, rater_type, edu_series, pc_index,
    cum_low = cum_persona_low, cum_med = cum_persona_med, cum_high = cum_persona_high,
    component = "Persona"
  )]
), use.names = TRUE)

df_comp[, component := factor(component, levels = c("Non-persona", "Persona"))]

pc_breaks  <- sort(unique(df_comp$pc_index))
pc_breaks5 <- pc_breaks[pc_breaks %% 5L == 0L]

# Re-impose facet ordering for plotted tables (defensive against coercions)
df_comp[, rater := factor(as.character(rater), levels = ordered_raters)]

# --- Background GSS median lines in every LLM facet ---------------------------

llm_raters_only <- setdiff(levels(nonpersona_cumvar_dt$rater), "gss")

gss_base_comp <- df_comp[rater == "gss", .(
  edu_series, pc_index, cum_med, component
)]

# cartesian product: replicate GSS curves into every LLM facet
llm_facets <- data.table::data.table(rater = llm_raters_only)
llm_facets[, key := 1L]
gss_base_comp[, key := 1L]

gss_bg_comp <- merge(
  llm_facets,
  gss_base_comp,
  by = "key",
  allow.cartesian = TRUE
)[, key := NULL]

# Re-impose facet ordering for background data (defensive)
gss_bg_comp[, rater := factor(as.character(rater), levels = ordered_raters)]

if (!isTRUE(get0("SKIP_V2B_PLOTS", ifnotfound = FALSE))) {
  out_file_non_cum <- file.path(viz_dir, sprintf("nonpersona_cumulative_pc_share_%s.pdf", year))
  grDevices::pdf(out_file_non_cum, width = 17.5, height = 15)  # landscape
  
  p <- 
  ggplot2::ggplot(df_comp, ggplot2::aes(x = pc_index, group = interaction(edu_series, component))) +
    ggplot2::geom_line(
      data = gss_bg_comp,
      ggplot2::aes(
        x = pc_index,
        y = cum_med,
        group = interaction(edu_series, component),
        colour = edu_series,
        linetype = component
      ),
      linewidth = 0.9,
      alpha = 0.12,
      show.legend = FALSE
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = cum_low, ymax = cum_high, fill = edu_series),
      alpha = 0.20, colour = NA
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = cum_med, colour = edu_series, linetype = component),
      linewidth = 0.7,
      alpha = 0.70
    ) +
    ggplot2::scale_colour_manual(name = "Rater × education", values = edu_colours) +
    ggplot2::scale_fill_manual(name = "Rater × education", values = edu_colours) +
    ggplot2::scale_linetype_manual(name = "Component", values = c("Non-persona" = "solid", "Persona" = "dashed")) +
    ggplot2::scale_x_continuous(breaks = pc_breaks5, minor_breaks = NULL) +
    ggplot2::coord_cartesian(ylim = c(0, max(df_comp$cum_high, na.rm = TRUE))) +
    ggplot2::facet_wrap(~ rater, ncol = 5) +
  ggplot2::labs(
    x     = "Principal component index",
    y     = "Cumulative share of total belief variance",
    title = ""
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    strip.text       = ggplot2::element_text(face = "bold"),
    legend.position  = "bottom"
  )
  
  print(p)
  grDevices::dev.off()
}

## ---- Plot 1b: nonpersona cumulative (all models in one panel) ------------ ##

if (nrow(df_comp)) {
  out_file_non_cum_all <- file.path(viz_dir, sprintf("nonpersona_cumulative_pc_share_allinone_%s.pdf", year))
  grDevices::pdf(out_file_non_cum_all, width = 10, height = 6)
  
  df_comp_overall <- df_comp[edu_group == "overall"]
  df_comp_overall[, is_gss := as.character(rater) == "gss"]
  threshold <- 0.9
  pc90_overall <- nonpersona_cumvar_dt[edu_group == "overall", .(
    rater, pc_index, cum_total_med
  )]
  pc90_overall <- pc90_overall[order(pc_index),
                               {
                                 idx <- which(cum_total_med >= threshold)
                                 if (length(idx)) .SD[idx[1L]] else NULL
                               },
                               by = rater]
  pc90_overall <- pc90_overall[order(pc_index)]
  top3 <- tail(pc90_overall$rater, 3)
  bot3 <- head(pc90_overall$rater, 3)
  keep_llm <- setdiff(unique(c(top3, bot3)), "gss")
  
  df_comp_overall[, highlight := as.character(rater) %in% keep_llm]
  df_comp_overall[, color_group := ifelse(highlight, as.character(rater), "other")]
  
  highlight_cols <- setNames(c("dodgerblue4", "dodgerblue2", "dodgerblue"), keep_llm[seq_len(min(3, length(keep_llm)))])
  if (length(keep_llm) > 3) {
    extra <- keep_llm[4:length(keep_llm)]
    highlight_cols <- c(highlight_cols, setNames(c("orangered3", "orangered2", "orange")[seq_len(length(extra))], extra))
  }
  highlight_cols <- c(highlight_cols, other = "grey70", gss = "black")
  
  p_all <- ggplot2::ggplot(
    df_comp_overall,
    ggplot2::aes(x = pc_index, y = cum_med, group = interaction(rater, component), colour = color_group)
  ) +
    ggplot2::geom_line(
      ggplot2::aes(linetype = component),
      linewidth = 0.7,
      alpha = 0.8
    ) +
    ggplot2::scale_colour_manual(values = highlight_cols, breaks = c(keep_llm, "gss", "other"), name = NULL) +
    ggplot2::scale_linetype_manual(name = "Component", values = c("Non-persona" = "solid", "Persona" = "dashed")) +
    ggplot2::scale_x_continuous(breaks = pc_breaks5, minor_breaks = NULL) +
    ggplot2::coord_cartesian(ylim = c(0, max(df_comp_overall$cum_high, na.rm = TRUE))) +
    ggplot2::labs(
      x = "Principal component index",
      y = "Cumulative share of total belief variance",
      title = ""
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "right"
    )
  
  print(p_all)
  grDevices::dev.off()
}


## ---- Plot 2: cumulative total share + PC count to reach 90% -------------- ##

df_total <- nonpersona_cumvar_dt[, .(
  rater, edu_group, rater_type, edu_series, pc_index,
  cum_total_low, cum_total_med, cum_total_high
)]

threshold <- 0.9

pc90_dt <- df_total[order(pc_index),
                    {
                      idx <- which(cum_total_med >= threshold)
                      if (length(idx)) .SD[idx[1L]] else NULL
                    },
                    by = .(rater, edu_series)
]

y_min_tot <- 0
y_max_tot <- max(df_total$cum_total_high, na.rm = TRUE)
y_range   <- y_max_tot - y_min_tot
y_base    <- y_min_tot + 0.04 * y_range

pc90_dt[, edu_rank := as.integer(factor(edu_series, levels = names(edu_colours)))]
mid_rank <- mean(range(pc90_dt$edu_rank, na.rm = TRUE))
pc90_dt[, y_lab := y_base + (edu_rank - mid_rank) * 0.03 * y_range]
pc90_dt[, y_lab := y_lab + 0.05 * y_range]

y_breaks <- sort(unique(c(threshold, pretty(c(y_min_tot, y_max_tot)))))

pc_breaks  <- sort(unique(df_total$pc_index))
pc_breaks5 <- pc_breaks[pc_breaks %% 5L == 0L]

# Re-impose facet ordering for plotted tables (defensive against coercions)
df_total[, rater := factor(as.character(rater), levels = ordered_raters)]
pc90_dt[,  rater := factor(as.character(rater), levels = ordered_raters)]

# --- Background GSS median (total) lines in every LLM facet -------------------

gss_base_total <- df_total[rater == "gss", .(
  edu_series, pc_index, cum_total_med
)]

llm_facets <- data.table::data.table(rater = llm_raters_only)
llm_facets[, key := 1L]
gss_base_total[, key := 1L]

gss_bg_total <- merge(
  llm_facets,
  gss_base_total,
  by = "key",
  allow.cartesian = TRUE
)[, key := NULL]

# Re-impose facet ordering for background data (defensive)
gss_bg_total[, rater := factor(as.character(rater), levels = ordered_raters)]

out_file_tot_cum <- file.path(viz_dir, sprintf("total_cumulative_pc_share_%s.pdf", year))
grDevices::pdf(out_file_tot_cum, width = 17.5, height = 15)

p <- 
ggplot2::ggplot(df_total, ggplot2::aes(x = pc_index)) +
  ggplot2::geom_line(
    data = gss_bg_total,
    ggplot2::aes(
      x = pc_index,
      y = cum_total_med,
      group = edu_series,
      colour = edu_series
    ),
    linewidth = 0.9,
    alpha = 0.12,
    show.legend = FALSE
  ) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = cum_total_low, ymax = cum_total_high, fill = edu_series),
    alpha = 0.20, colour = NA
  ) +
  ggplot2::geom_line(
    ggplot2::aes(y = cum_total_med, colour = edu_series),
    linewidth = 0.7,
    alpha = 0.70
  ) +
  ggplot2::geom_hline(yintercept = threshold, linetype = "dotted") +
  ggplot2::geom_vline(
    data = pc90_dt,
    ggplot2::aes(xintercept = pc_index, colour = edu_series),
    linetype = "dotted", linewidth = 0.6, show.legend = FALSE
  ) +
  ggplot2::geom_text(
    data = pc90_dt,
    ggplot2::aes(x = pc_index, y = y_lab, label = pc_index, colour = edu_series),
    fontface = "bold", size = 3.8, vjust = 0.5, show.legend = FALSE
  ) +
  ggplot2::scale_colour_manual(name = "Rater × education", values = edu_colours) +
  ggplot2::scale_fill_manual(name = "Rater × education", values = edu_colours) +
  ggplot2::scale_x_continuous(breaks = pc_breaks5, minor_breaks = NULL) +
  ggplot2::scale_y_continuous(breaks = y_breaks, minor_breaks = NULL) +
  ggplot2::coord_cartesian(ylim = c(y_min_tot, y_max_tot)) +
  ggplot2::facet_wrap(~ rater, ncol = 5) +
  ggplot2::labs(
    x     = "Principal component index",
    y     = "Cumulative share of total belief variance",
    title = ""
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    strip.text       = ggplot2::element_text(face = "bold"),
    plot.title       = ggplot2::element_text(face = "bold"),
    legend.position  = "bottom"
  )
print(p)
grDevices::dev.off()

invisible(NULL)
