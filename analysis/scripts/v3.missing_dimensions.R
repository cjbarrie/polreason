################################################################################
# Rayleigh-ratio variance diagnostics: LLM vs GSS along a shared GSS PCA basis
#
# Overview
# --------
# This script compares how much variance different raters (LLMs and GSS) exhibit
# along a *shared* belief-system basis derived from GSS correlation structure.
#
# We:
#   1) Identify which raters have data on disk for YEAR (via available_raters()).
#   2) Build an LLM-anchored “eligible belief universe”:
#        - union of belief variables appearing in at least one LLM bootstrap draw,
#          after persona handling (drop vs conditional).
#   3) Compute the belief-variable overlap across all raters (including GSS),
#      restricted to the LLM-eligible union, so the basis is not anchored on
#      GSS-only items.
#   4) Build a GSS PCA basis on the overlap set:
#        - load GSS bootstrap correlation matrices,
#        - (optionally) align to a modal variable set (required for averaging),
#        - average across draws,
#        - apply persona handling,
#        - eigendecompose the resulting belief correlation matrix,
#        - retain K PCs to reach CUMVAR_TARGET cumulative variance.
#   5) For each rater and each bootstrap draw that contains the full basis varset,
#      compute Rayleigh quotients along each retained GSS PC direction:
#         v_j(rater, draw) = u_j' R_rater(draw) u_j
#      and compare to the corresponding GSS eigenvalue:
#         ratio_j = v_j(rater, draw) / lambda_j(GSS)
#   6) Summarize ratios across bootstraps (median and 95% interval) and plot a
#      heatmap-like rectangle plot where PC segments have widths proportional
#      to GSS explained-variance shares.
#
# Persona handling
# ---------------
# A single switch controls whether we:
#   - CONDITIONAL = TRUE : compute belief correlations conditional on persona vars
#                          using cond_belief_corr_on_persona(), or
#   - CONDITIONAL = FALSE: simply drop persona vars via drop_persona_vars().
#
# IMPORTANT: Some upstream conditioning utilities may drop dimnames. This script
# repairs dimnames to prevent silent misalignment across raters.
#
# Bootstrap “as-is” policy
# ------------------------
# For union/overlap discovery and Rayleigh ratio computation, we load correlation
# lists with align_to_mode_varset = FALSE to preserve each draw’s native varset.
# We ONLY align to a modal varset when we must average matrices (GSS basis step).
#
# Prerequisites (provided elsewhere in the project)
# ------------------------------------------------
# Required globals:
#   - BASE_OUT_DIR
#   - BASE_VIZ_DIR
#   - persona_vars_global : character vector of persona variable names
#
# Required utils (assumed already sourced)
# ---------------------------------------
#   - available_raters(base_out_dir, year, exclude=...)
#   - load_corr_for_rater(rater, base_out_dir, year, suffix, align_to_mode_varset, strict)
#   - drop_persona_vars(R_full, persona_vars=...)
#   - cond_belief_corr_on_persona(R_full, persona_vars=...)
#
# Packages used
# -------------
#   - data.table
#   - ggplot2
#   - jsonlite
#   - scales
#
# Output
# ------
# Writes a PDF to:
#   OUTFILE_PDF = file.path(BASE_VIZ_DIR, sprintf("missing_variance_rayleigh_%s.pdf", YEAR))
################################################################################



## ------------------------------------------------------------------
## User-facing settings
## ------------------------------------------------------------------
YEAR              <- get0("YEAR", ifnotfound = 2024L)
RATER_GSS         <- "gss"

# Persona handling:
#   TRUE  => partial correlations among belief vars conditional on persona vars
#   FALSE => just drop persona vars
CONDITIONAL       <- TRUE

# How many GSS PCs to retain for computation (cumulative explained variance)
CUMVAR_TARGET     <- 0.90

# For plotting: only show PCs that explain more than this fraction (within first K)
MIN_PC_SHARE      <- 0.01

# Output
OUTFILE_PDF       <- file.path(BASE_VIZ_DIR, sprintf("missing_variance_rayleigh_%s.pdf", YEAR))
PDF_WIDTH         <- 22.5
PDF_HEIGHT        <- 10

## ------------------------------------------------------------------
## Optional: human-readable PC labels (edit as appropriate)
## ------------------------------------------------------------------
# These labels are interpretive; they do not affect the numerical analysis.
pc_labels_tmp <- jsonlite::fromJSON(txt = 'analysis/GSS_PC_explain.json', simplifyVector = TRUE)
pc_labels <- pc_labels_tmp$axis_label

## ------------------------------------------------------------------
## Sanity checks for prerequisites
## ------------------------------------------------------------------
required_objects <- c("BASE_OUT_DIR", "BASE_VIZ_DIR", "persona_vars_global")
missing_objects  <- required_objects[!vapply(required_objects, exists, logical(1), inherits = TRUE)]
if (length(missing_objects)) {
  stop("Missing required objects in the workspace: ", paste(missing_objects, collapse = ", "))
}

required_fns <- c(
  "load_corr_for_rater",
  "available_raters",
  "drop_persona_vars",
  "cond_belief_corr_on_persona"
)
missing_fns <- required_fns[!vapply(required_fns, exists, logical(1), mode = "function", inherits = TRUE)]
if (length(missing_fns)) {
  stop("Missing required functions in the workspace: ", paste(missing_fns, collapse = ", "))
}

## ------------------------------------------------------------------
## Helper utilities (script-local)
## ------------------------------------------------------------------

#' Extract variable names from a matrix robustly.
#'
#' @description
#' Returns a unique, non-empty set of variable names using column names when
#' available, otherwise row names. Filters out `NA` and empty strings.
#'
#' @param M A matrix-like object, typically a correlation matrix.
#'
#' @return A character vector of unique variable names (possibly length 0).
#'
#' @keywords internal
#' @noRd
.safe_varnames <- function(M) {
  v <- colnames(M)
  if (is.null(v) || !length(v)) v <- rownames(M)
  v <- v[!is.na(v) & nzchar(v)]
  unique(v)
}

#' Apply persona handling (drop or condition) and repair dimnames if needed.
#'
#' @description
#' Produces a belief-by-belief correlation matrix with correct dimnames.
#'
#' When `conditional = TRUE`, this attempts to condition beliefs on persona vars
#' via `cond_belief_corr_on_persona()` when persona variables are present in the
#' input matrix. If no persona vars are present, it falls back to dropping persona
#' variables (with a warning).
#'
#' When `conditional = FALSE`, persona variables are dropped via `drop_persona_vars()`.
#'
#' Some upstream conditioning utilities may return matrices without dimnames. This
#' function restores dimnames based on the original variable order to prevent
#' silent misalignment across raters.
#'
#' @param R_full Full correlation matrix that may include beliefs and persona vars.
#' @param conditional Logical; whether to condition beliefs on persona vars.
#' @param persona_vars Character vector of persona variable names (canonical set).
#'
#' @return A numeric matrix: belief-by-belief correlation matrix with dimnames,
#'   or `NULL` if input is unusable.
#'
#' @keywords internal
#' @noRd
.apply_persona <- function(R_full, conditional, persona_vars) {
  R_full <- as.matrix(R_full)
  
  # Recover the full variable order from the input matrix.
  full_vars <- .safe_varnames(R_full)
  if (!length(full_vars)) return(NULL)
  
  # Identify which persona vars are present in this matrix.
  persona_present <- intersect(persona_vars, full_vars)
  
  # Belief variables are everything except persona vars, in *original order*.
  belief_vars <- full_vars[!(full_vars %in% persona_present)]
  
  if (!conditional) {
    out <- drop_persona_vars(R_full, persona_vars = persona_vars)
    if (is.null(out)) return(NULL)
    out <- as.matrix(out)
    
    # If names were dropped upstream, restore them when dimensions match.
    if (is.null(colnames(out)) || !length(colnames(out))) {
      if (ncol(out) == length(belief_vars)) {
        dimnames(out) <- list(belief_vars, belief_vars)
      }
    }
    return(out)
  }
  
  # conditional == TRUE
  if (!length(persona_present)) {
    warning("No persona vars present in this matrix; falling back to drop_persona_vars().")
    out <- drop_persona_vars(R_full, persona_vars = persona_vars)
  } else {
    out <- cond_belief_corr_on_persona(R_full, persona_vars = persona_present)
  }
  
  if (is.null(out)) return(NULL)
  out <- as.matrix(out)
  
  # Re-attach names if the conditioning function dropped dimnames.
  if (is.null(colnames(out)) || !length(colnames(out))) {
    if (ncol(out) == length(belief_vars)) {
      dimnames(out) <- list(belief_vars, belief_vars)
    } else if (ncol(out) == length(full_vars)) {
      dimnames(out) <- list(full_vars, full_vars)
    } else {
      stop(
        "Conditioned matrix has no dimnames and unexpected size: ncol(out) = ",
        ncol(out), " (expected ", length(belief_vars), " or ", length(full_vars), ")."
      )
    }
  }
  
  # If conditioning returned a full matrix (belief+persona), drop persona here.
  if (length(persona_present) && any(persona_present %in% colnames(out))) {
    out <- out[belief_vars, belief_vars, drop = FALSE]
  }
  
  out
}

## ------------------------------------------------------------------
## Step 1: Determine the variable universe to analyze
## ------------------------------------------------------------------

#' Compute the union of belief variables observed in LLM raters after persona handling.
#'
#' @description
#' Iterates over all LLM raters and their bootstrap draws. For each draw, applies
#' persona handling (drop vs conditional) and collects the set of belief variable
#' names. The result is the union across all LLM raters and draws.
#'
#' Bootstrap draws are kept “as-is”: correlation lists are loaded with
#' `align_to_mode_varset = FALSE` to avoid shrinking the variable universe.
#'
#' @param raters_llm Character vector of LLM rater identifiers (must be non-empty).
#' @param base_out_dir Character; base output directory.
#' @param year Integer; analysis year.
#' @param persona_vars Character vector of persona variable names.
#' @param conditional Logical; whether to condition on persona vars (TRUE) or drop them (FALSE).
#'
#' @return Character vector of belief variable names appearing in at least one LLM draw.
#'
#' @export
get_llm_union_belief_vars <- function(raters_llm, base_out_dir, year, persona_vars, conditional) {
  if (!length(raters_llm)) stop("raters_llm is empty.")
  
  union_vars <- character(0)
  
  for (r in raters_llm) {
    cl <- load_corr_for_rater(
      rater = r,
      base_out_dir = base_out_dir,
      year = year,
      align_to_mode_varset = FALSE,
      strict = FALSE
    )
    if (is.null(cl) || !length(cl)) next
    
    for (R_full in cl) {
      if (is.null(R_full)) next
      R_use <- .apply_persona(R_full, conditional = conditional, persona_vars = persona_vars)
      if (is.null(R_use)) next
      
      vv <- .safe_varnames(R_use)
      if (!length(vv)) next
      union_vars <- union(union_vars, vv)
    }
  }
  
  if (!length(union_vars)) {
    stop("No belief variables found in any LLM after persona handling. ",
         "This typically indicates empty matrices or missing dimnames upstream.")
  }
  
  union_vars
}

#' Compute the overlapping belief-variable set across raters, restricted to the LLM union.
#'
#' @description
#' For each rater, loads its correlation list and uses a representative draw to
#' obtain the set of belief variable names after persona handling. Each rater’s
#' varset is intersected with `llm_union_vars`, then overlap across raters is
#' computed.
#'
#' For stability, this function loads lists with `align_to_mode_varset = TRUE`
#' (if available in your environment) so the representative draw reflects the
#' rater’s modal variable set.
#'
#' @param raters Character vector of rater identifiers (including GSS).
#' @param llm_union_vars Character vector of belief variables eligible from LLMs.
#' @param base_out_dir Character; base output directory.
#' @param year Integer; analysis year.
#' @param persona_vars Character vector of persona variable names.
#' @param conditional Logical; whether to condition on persona vars (TRUE) or drop them (FALSE).
#'
#' @return Character vector of overlapping belief variables across raters after
#'   restricting each rater to `llm_union_vars`.
#'
#' @export
get_overlapping_belief_vars <- function(raters, llm_union_vars,
                                        base_out_dir, year, persona_vars, conditional) {
  if (is.null(llm_union_vars) || !length(llm_union_vars)) {
    stop("llm_union_vars must be non-empty.")
  }
  
  varsets <- list()
  
  for (r in raters) {
    cl <- load_corr_for_rater(
      rater = r,
      base_out_dir = base_out_dir,
      year = year,
      align_to_mode_varset = TRUE,
      strict = FALSE
    )
    if (is.null(cl) || !length(cl)) next
    
    # Representative draw for names after modal-varset alignment.
    R_full <- cl[[1L]]
    if (is.null(R_full)) next
    
    R_use <- .apply_persona(R_full, conditional = conditional, persona_vars = persona_vars)
    if (is.null(R_use)) next
    
    vv <- intersect(.safe_varnames(R_use), llm_union_vars)
    if (length(vv) < 2L) next
    
    varsets[[r]] <- vv
  }
  
  if (!length(varsets)) {
    stop("No usable raters when computing overlapping variable set.")
  }
  
  overlap <- Reduce(intersect, varsets)
  if (!length(overlap)) {
    stop("No overlapping belief variables found across raters after restricting to LLM-union.")
  }
  
  overlap
}

## ------------------------------------------------------------------
## Step 2: Build the GSS basis (eigendecomposition on overlap vars)
## ------------------------------------------------------------------

#' Build a GSS PCA basis on an overlapping, LLM-eligible belief-variable set.
#'
#' @description
#' Loads the GSS bootstrap correlation matrices, aligns them to a modal variable
#' set if possible (required to average), averages across draws, applies persona
#' handling (drop vs conditional), subsets to the overlap vars in GSS order,
#' symmetrizes for numerical stability, then eigendecomposes the resulting
#' correlation matrix.
#'
#' Eigenvalues are clamped at 0 for explained-variance accounting.
#' The number of PCs `K` is chosen as the smallest index reaching `cumvar_target`.
#'
#' @param overlap_vars Character vector of candidate overlap belief vars.
#' @param base_out_dir Character; base output directory.
#' @param year Integer; analysis year.
#' @param persona_vars Character vector of persona variable names.
#' @param conditional Logical; whether to condition on persona vars (TRUE) or drop them (FALSE).
#' @param cumvar_target Numeric in (0,1]; target cumulative variance share for retaining PCs.
#'
#' @return A list with elements:
#'   - `R_use`: correlation matrix used for PCA (belief-only)
#'   - `varnames`: variable names in basis order
#'   - `U`: eigenvectors (columns)
#'   - `lambda`: nonnegative eigenvalues
#'   - `cumvar`: cumulative variance shares
#'   - `K`: number of PCs to retain to reach target
#'   - `year`, `conditional`: metadata
#'
#' @export
build_gss_basis <- function(overlap_vars, base_out_dir, year, persona_vars,
                            conditional, cumvar_target = 0.90) {
  
  cl_gss <- load_corr_for_rater(
    rater = RATER_GSS,
    base_out_dir = base_out_dir,
    year = year,
    align_to_mode_varset = TRUE,  # required for averaging across draws
    strict = TRUE
  )
  
  if (is.null(cl_gss) || !length(cl_gss)) stop("No correlation list found for GSS.")
  
  # Average GSS correlation matrices across bootstrap draws.
  R_full_avg <- Reduce(`+`, cl_gss) / length(cl_gss)
  
  # Apply persona handling (drop vs conditional) with dimname repair.
  R_use_full <- .apply_persona(R_full_avg, conditional = conditional, persona_vars = persona_vars)
  if (is.null(R_use_full)) stop("Persona-handled GSS matrix is NULL.")
  
  # Enforce that overlap vars exist in this matrix, and adopt GSS ordering.
  overlap_vars <- intersect(overlap_vars, .safe_varnames(R_use_full))
  overlap_vars <- colnames(R_use_full)[colnames(R_use_full) %in% overlap_vars]
  
  if (length(overlap_vars) < 2L) stop("Need at least 2 overlapping belief vars to build a basis.")
  
  R_use <- R_use_full[overlap_vars, overlap_vars, drop = FALSE]
  
  # Numerical hygiene: ensure symmetry (important for eigen(..., symmetric=TRUE))
  R_use <- (R_use + t(R_use)) / 2
  
  eig <- eigen(R_use, symmetric = TRUE)
  U_raw <- eig$vectors
  lambda_raw <- eig$values
  
  # Clamp negative eigenvalues for explained-variance accounting.
  lambda <- pmax(Re(lambda_raw), 0)
  
  cumvar <- cumsum(lambda) / sum(lambda)
  K <- which(cumvar >= cumvar_target)[1L]
  if (is.na(K)) K <- length(lambda)
  
  list(
    R_use    = R_use,
    varnames = overlap_vars,
    U        = U_raw,
    lambda   = lambda,
    cumvar   = cumvar,
    K        = K,
    year     = year,
    conditional = conditional
  )
}

## ------------------------------------------------------------------
## Step 3: Rayleigh ratios along the GSS basis for a given rater
## ------------------------------------------------------------------

#' Compute Rayleigh ratios Var(rater)/Var(GSS) along a GSS basis for one rater.
#'
#' @description
#' For each bootstrap correlation matrix for a given rater:
#'   - load draws “as-is” (no mode-varset alignment),
#'   - apply persona handling (drop vs conditional),
#'   - require that all basis variables are present,
#'   - subset and reorder to the basis var order,
#'   - symmetrize,
#'   - compute Rayleigh quotients v_j = u_j' R u_j for retained PCs,
#'   - divide by the corresponding GSS eigenvalues lambda_j.
#'
#' PCs with non-positive GSS eigenvalues are excluded to avoid division-by-zero.
#'
#' @param rater Character; rater identifier.
#' @param gss_basis List returned by [build_gss_basis()].
#' @param base_out_dir Character; base output directory.
#' @param year Integer; analysis year.
#' @param persona_vars Character vector of persona variable names.
#' @param conditional Logical; must match `gss_basis$conditional`.
#' @param K_use Optional integer; number of PCs to consider (defaults to `gss_basis$K`).
#'
#' @return A `data.table` with columns:
#'   `rater`, `bootstrap_id`, `pc_index`, `var_rater`, `var_gss`, `ratio`;
#'   or `NULL` if no usable draws are found.
#'
#' @export
rayleigh_ratios_for_rater <- function(rater, gss_basis, base_out_dir, year,
                                      persona_vars, conditional, K_use = NULL) {
  if (!identical(conditional, gss_basis$conditional)) {
    stop("conditional must match gss_basis$conditional.")
  }
  
  cl <- load_corr_for_rater(
    rater = rater,
    base_out_dir = base_out_dir,
    year = year,
    align_to_mode_varset = FALSE, # preserve bootstrap varsets
    strict = FALSE
  )
  if (is.null(cl) || !length(cl)) return(NULL)
  
  var_ref <- gss_basis$varnames
  U_ref   <- gss_basis$U
  lam_ref <- gss_basis$lambda
  
  if (is.null(K_use)) K_use <- gss_basis$K
  K_use <- min(K_use, ncol(U_ref), length(lam_ref))
  
  # Use only PCs with strictly positive reference eigenvalues.
  pc_ok <- which(lam_ref[seq_len(K_use)] > 0)
  if (!length(pc_ok)) return(NULL)
  
  U_K   <- U_ref[, pc_ok, drop = FALSE]
  lam_K <- lam_ref[pc_ok]
  
  out_list <- vector("list", length(cl))
  keep <- 0L
  
  for (b in seq_along(cl)) {
    R_full <- cl[[b]]
    if (is.null(R_full)) next
    
    R_use <- .apply_persona(R_full, conditional = conditional, persona_vars = persona_vars)
    if (is.null(R_use)) next
    
    # Only retain bootstraps that contain all basis variables.
    if (!all(var_ref %in% colnames(R_use))) next
    
    # Restrict to basis variables in the basis ordering.
    R_use <- R_use[var_ref, var_ref, drop = FALSE]
    R_use <- (R_use + t(R_use)) / 2
    
    # Rayleigh quotient along each basis direction:
    #   v = diag(U' R U)
    # Compute as column-wise inner products:
    #   v_j = sum_i U_{ij} * (R U)_{ij}
    RU <- R_use %*% U_K
    v  <- colSums(U_K * RU)
    
    keep <- keep + 1L
    out_list[[keep]] <- data.table::data.table(
      rater        = rater,
      bootstrap_id = b,
      pc_index     = pc_ok,
      var_rater    = as.numeric(v),
      var_gss      = as.numeric(lam_K),
      ratio        = as.numeric(v / lam_K)
    )
  }
  
  if (!keep) {
    warning("Rater ", rater, ": no bootstrap draw contained the full basis varset; skipping.")
    return(NULL)
  }
  
  data.table::rbindlist(out_list[seq_len(keep)])
}

#' Summarize Rayleigh ratios across bootstrap draws for plotting.
#'
#' @description
#' Aggregates per-(rater, pc_index) statistics: median ratio, 95% interval, and
#' number of bootstrap draws used.
#'
#' @param ratio_dt A `data.table` as returned by [rayleigh_ratios_for_rater()].
#'
#' @return A `data.table` keyed by `rater` and `pc_index` with summary columns:
#'   `ratio_med`, `ratio_low`, `ratio_high`, `n_boot`.
#'
#' @export
summarize_ratios <- function(ratio_dt) {
  ratio_dt[
    ,
    .(
      ratio_med  = median(ratio, na.rm = TRUE),
      ratio_low  = as.numeric(stats::quantile(ratio, 0.025, na.rm = TRUE)),
      ratio_high = as.numeric(stats::quantile(ratio, 0.975, na.rm = TRUE)),
      n_boot     = .N
    ),
    by = .(rater, pc_index)
  ]
}

## ------------------------------------------------------------------
## Step 4: Plot construction helpers
## ------------------------------------------------------------------

#' Build PC geometry for plotting (segment widths, midpoints, and labels).
#'
#' @description
#' Computes plot geometry for each retained PC:
#'   - `x_min`, `x_max`: cumulative-width segment bounds on [0,1]
#'   - `x_mid`: segment midpoint
#'
#' Widths are proportional to each PC's share of variance within the retained
#' set of PCs (`K = gss_basis$K`). PCs with share <= `min_pc_share` are excluded.
#'
#' Optional `pc_labels` can override ordinal labels (1st, 2nd, ...).
#'
#' @param gss_basis List returned by [build_gss_basis()].
#' @param min_pc_share Numeric; minimum within-K variance share to include a PC.
#' @param pc_labels Optional character vector; label for PC k is `pc_labels[[k]]`.
#'
#' @return A `data.table` with one row per plotted PC and columns:
#'   `pc_index`, `x_min`, `x_max`, `x_mid`, `pct_label`, `pc_label`.
#'
#' @export
build_pc_geometry <- function(gss_basis, min_pc_share = 0.01, pc_labels = NULL) {
  K <- gss_basis$K
  lambda <- gss_basis$lambda[seq_len(K)]
  
  # Share of variance among the retained K PCs.
  share_all <- lambda / sum(lambda)
  
  pc_keep <- which(share_all > min_pc_share)
  if (!length(pc_keep)) stop("No PCs exceed min_pc_share; lower the threshold.")
  
  share_keep_raw <- share_all[pc_keep]                     # for % labels
  share_keep_w   <- share_keep_raw / sum(share_keep_raw)   # widths sum to 1
  
  cum_w <- cumsum(share_keep_w)
  x_max <- cum_w
  x_min <- c(0, head(cum_w, -1))
  x_mid <- (x_min + x_max) / 2
  
  ordinal_label <- function(k) {
    k <- as.integer(k)
    if ((k %% 100L) %in% 11L:13L) return(paste0(k, "th"))
    last <- k %% 10L
    suf <- if (last == 1L) "st" else if (last == 2L) "nd" else if (last == 3L) "rd" else "th"
    paste0(k, suf)
  }
  
  get_label <- function(k) {
    if (is.character(pc_labels) &&
        length(pc_labels) >= k &&
        !is.na(pc_labels[[k]]) &&
        nzchar(pc_labels[[k]])) {
      return(pc_labels[[k]])
    }
    ordinal_label(k)
  }
  
  data.table::data.table(
    pc_index  = pc_keep,
    x_min     = x_min,
    x_max     = x_max,
    x_mid     = x_mid,
    pct_label = paste0(round(share_keep_raw * 100, 1), "%"),
    pc_label  = vapply(pc_keep, get_label, character(1))
  )
}

#' Create the Rayleigh-ratio rectangle plot.
#'
#' @description
#' Produces a ggplot2 figure where each rectangle corresponds to a (rater, PC)
#' cell. The x-axis is partitioned into PC segments whose widths reflect each
#' PC's variance share in the retained GSS basis. The fill value is the median
#' Rayleigh ratio per cell.
#'
#' @param summary_dt A `data.table` from [summarize_ratios()].
#' @param pc_geom_dt A `data.table` from [build_pc_geometry()].
#' @param rater_order Character vector giving the desired y-axis order.
#'
#' @return A `ggplot` object.
#'
#' @export
make_rayleigh_plot <- function(summary_dt, pc_geom_dt, rater_order) {
  # Restrict to plotted PCs and merge in geometry.
  plot_dt <- merge(
    summary_dt[pc_index %in% pc_geom_dt$pc_index],
    pc_geom_dt,
    by = "pc_index",
    all.x = TRUE
  )
  
  # Impose rater order as a factor so the y-axis is stable and interpretable.
  plot_dt[, rater := factor(rater, levels = rater_order)]
  plot_dt[, rater_num := as.integer(rater)]
  
  plot_dt[, fill_val := ratio_med]
  plot_df <- as.data.frame(plot_dt)
  
  max_y <- max(plot_df$rater_num, na.rm = TRUE)
  min_y <- min(plot_df$rater_num, na.rm = TRUE)
  
  pc_top_df   <- as.data.frame(unique(pc_geom_dt[, .(x_mid, pct_label)]))
  pc_label_df <- as.data.frame(pc_geom_dt[, .(x_mid, pc_label)])
  
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      xmin = x_min, xmax = x_max,
      ymin = rater_num - 0.5, ymax = rater_num + 0.5,
      fill = fill_val
    )
  ) +
    ggplot2::geom_rect() +
    ggplot2::geom_text(
      data = pc_top_df,
      ggplot2::aes(x = x_mid, y = max_y + 1, label = pct_label),
      inherit.aes = FALSE,
      size = 4
    ) +
    ggplot2::geom_text(
      data = pc_label_df,
      ggplot2::aes(x = x_mid, y = min_y - 0.5, label = pc_label),
      inherit.aes = FALSE,
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size  = 4.5
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq_along(rater_order),
      labels = rater_order,
      name   = "Rater",
      expand = c(0, 0),
      limits = c(min_y - 1.5, max_y + 1.5)
    ) +
    ggplot2::scale_x_continuous(
      breaks = pc_geom_dt$x_mid,
      labels = NULL,
      name   = NULL,
      expand = c(0, 0),
      limits = c(0, 1)
    ) +
    ggplot2::scale_fill_gradientn(
      colours = c("orangered", "white", "dodgerblue", "purple4"),
      values  = c(0, 1/3, 2/3, 1),
      limits  = c(0, 3),
      breaks  = c(0, 1, 2, 3),
      labels  = c("0", "1", "2", "3"),
      oob     = scales::squish,
      name    = "Rayleigh ratio\nVar(LLM) / Var(GSS)"
    ) +
    ggplot2::labs(
      title = "Rayleigh ratios of LLM v. GSS variance over shared GSS belief-system basis.",
      subtitle = if (isTRUE(get0("CONDITIONAL", ifnotfound = CONDITIONAL))) {
        "Basis calculated on belief-systems conditional on persona variables."
      } else {
        "Basis calculated on belief-systems with persona variables dropped."
      }
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y  = ggplot2::element_text(size = 12),
      axis.text.x  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid   = ggplot2::element_blank(),
      plot.margin  = ggplot2::margin(t = 20, r = 10, b = 300, l = 10),
      axis.title.x = ggplot2::element_blank()
    )
}

## ------------------------------------------------------------------
## Main pipeline
## ------------------------------------------------------------------

# Identify the set of raters with available data.
raters_all <- available_raters(BASE_OUT_DIR, YEAR)

# Ensure GSS is included even if available_raters() excludes it for any reason.
raters_all <- unique(c(raters_all, RATER_GSS))

raters_llm <- setdiff(raters_all, RATER_GSS)

# 1) Union of belief variables across LLMs after persona handling.
llm_union_vars <- get_llm_union_belief_vars(
  raters_llm   = raters_llm,
  base_out_dir = BASE_OUT_DIR,
  year         = YEAR,
  persona_vars = persona_vars_global,
  conditional  = CONDITIONAL
)

# 2) Overlap across raters (including GSS), restricted to LLM union.
overlap_vars <- get_overlapping_belief_vars(
  raters         = raters_all,
  llm_union_vars = llm_union_vars,
  base_out_dir   = BASE_OUT_DIR,
  year           = YEAR,
  persona_vars   = persona_vars_global,
  conditional    = CONDITIONAL
)

# 3) Build the GSS basis on the overlap vars.
gss_basis <- build_gss_basis(
  overlap_vars   = overlap_vars,
  base_out_dir   = BASE_OUT_DIR,
  year           = YEAR,
  persona_vars   = persona_vars_global,
  conditional    = CONDITIONAL,
  cumvar_target  = CUMVAR_TARGET
)

message(
  "GSS basis built on ", length(gss_basis$varnames), " belief variables; ",
  "K = ", gss_basis$K, " PCs explain ",
  round(100 * gss_basis$cumvar[gss_basis$K], 1), "% variance."
)

# 4) Compute Rayleigh ratios for each rater, then summarize across bootstraps.
summary_list <- vector("list", length(raters_all))
names(summary_list) <- raters_all

for (i in seq_along(raters_all)) {
  r <- raters_all[[i]]
  
  ratio_dt <- rayleigh_ratios_for_rater(
    rater        = r,
    gss_basis    = gss_basis,
    base_out_dir = BASE_OUT_DIR,
    year         = YEAR,
    persona_vars = persona_vars_global,
    conditional  = CONDITIONAL,
    K_use        = gss_basis$K
  )
  
  if (is.null(ratio_dt)) next
  summary_list[[r]] <- summarize_ratios(ratio_dt)
}

summary_dt <- data.table::rbindlist(summary_list, use.names = TRUE, fill = TRUE)
if (!nrow(summary_dt)) stop("No summary results produced; check rater inputs and overlap vars.")

# 5) Order raters: put GSS first, then sort by closeness to 1 (mean squared deviation).
rater_dist <- summary_dt[, .(dist = mean((ratio_med - 1)^2, na.rm = TRUE)), by = rater]
rater_order <- c(
  RATER_GSS,
  rater_dist[rater != RATER_GSS][order(dist)]$rater
)
rater_order <- unique(rater_order)

# 6) Build PC geometry for the plot (widths + labels).
pc_geom_dt <- build_pc_geometry(
  gss_basis     = gss_basis,
  min_pc_share  = MIN_PC_SHARE,
  pc_labels     = pc_labels
)

# 7) Make plot.
p <- make_rayleigh_plot(summary_dt, pc_geom_dt, rater_order)

# 8) Save to PDF.
grDevices::pdf(OUTFILE_PDF, width = PDF_WIDTH, height = PDF_HEIGHT)
print(p)
grDevices::dev.off()

message("Wrote Rayleigh ratio plot to: ", OUTFILE_PDF)

## ------------------------------------------------------------------
## (Optional) PC interpretation helper
## ------------------------------------------------------------------

#' Inspect variable loadings for a given GSS PC.
#'
#' @description
#' Extracts the eigenvector loadings for PC `pc_j` from a GSS basis object and
#' returns a table sorted by absolute loading magnitude (descending).
#'
#' @param gss_basis List returned by [build_gss_basis()].
#' @param pc_j Integer; which principal component to inspect (column index into `gss_basis$U`).
#'
#' @return A `data.table` with columns `var` and `loading`, sorted by `abs(loading)` desc.
#'
#' @export
interpret_pc <- function(gss_basis, pc_j) {
  loadings <- gss_basis$U[, pc_j]
  data.table::data.table(
    var     = gss_basis$varnames,
    loading = loadings
  )[order(-abs(loading))]
}

for (i in seq_len(ncol(gss_basis$U))) {
  print(paste0("PC", i, ":"))
  print(interpret_pc(gss_basis, pc_j = i))
}

# The following prompt and the associated PCs are passed to GPT PRO 5.2 :
#
# "You are an expert quantitative social scientist interpreting PCA/FA loadings 
# from the General Social Survey (GSS). You will be given one or more principal 
# components (PCs) as a list of {var, loading}. Variable names are GSS codebook 
# variable codes, and a GSS codebook PDF is attached. Your job is to assign each 
# PC a short, human-meaningful label. 
#
# Core requirements 
# 1) USE THE CODEBOOK: For every variable you rely on, look it up in the 
# attached GSS codebook (or the provided mapping) and use its substantive 
# meaning. Do not guess variable meanings. 
# 2) Focus on the signal: identify the strongest loadings by absolute value; 
# ignore tiny loadings unless needed for clarity. 
# 3) Handle sign correctly: PC direction is arbitrary. Treat the PC as an axis. 
# Provide (a) an overall axis label and (b) a label for each pole based on which 
# variables are positive vs negative. 
# 4) Labels must be SHORT: Every label you output (axis label, positive pole 
# label, negative pole label, alternates) must be ≤5 words. 
# 5) Be honest about uncertainty: if too many key variables cannot be found in 
# the codebook, or the component mixes unrelated themes, output an 
# “unclear/mixed” label (still ≤5 words) and set confidence to Low. 
#
# Method  (do this internally; don’t show step-by-step) 
# A) For each PC, select: 
# - Top POS variables: up to 6 with largest positive loadings 
# - Top NEG variables: up to 6 with most negative loadings 
# - Prefer variables with |loading| ≥ 0.10 (if fewer than 4 meet this, just use 
# top 4 by |loading|) 
# B) Look up each selected variable in the GSS codebook and write a 3–10 word 
# “meaning gloss” for it (internally). 
# C) Infer the latent dimension that best summarizes the contrast between the 
# POS and NEG groups. 
# D) Propose 3 candidate axis labels (≤5 words), choose the best one, and also 
# provide ≤5-word labels for each pole.
# 
# Output format (STRICT) Return ONLY valid JSON (no markdown, no commentary) as 
# an array of objects, one per PC, using this schema: 
# [ 
#  { "pc": "PC1", 
#    "axis_label": "...", // ≤5 words 
#    "positive_pole_label": "...", // ≤5 words 
#    "negative_pole_label": "...", // ≤5 words 
#    "alt_axis_labels": ["...", "..."], // two alternates, each ≤5 words 
#    "top_positive_vars": [{"var":"...", "loading":0.0}, ...], 
#    "top_negative_vars": [{"var":"...", "loading":0.0}, ...], 
#    "unknown_vars": ["..."], // any selected vars not found in codebook 
#    "confidence": "High|Medium|Low", 
#    "rationale": "1–2 sentences explaining the theme using the codebook meanings."
#  } 
# ]
# 
# Style constraints 
# - Don’t moralize or editorialize; 
# keep labels descriptive. 
# - Avoid jargon unless it’s standard in social science. 
# - If the theme is mainly one domain (e.g., “religiosity”), the axis label may 
#   be a single concept; otherwise use a contrast label like “X vs Y” (counts as words).
# - Double-check the ≤5-word rule before finalizing labels. 
#
# INPUT STARTS NOW.




## ------------------------------------------------------------------
## FULL grid: ALL PCs + ALL variables
##   - tile fill = abs(loading) (0..1), white -> orangered
##   - text label = ORIGINAL signed loading (NOT abs)
## ------------------------------------------------------------------

plot_pc_loadings_grid_full_signed_labels <- function(gss_basis,
                                                     value_digits = 2,
                                                     order_vars_by = c("max_abs", "original"),
                                                     order_pcs_by  = c("eigen_order", "original")) {
  stopifnot(!is.null(gss_basis$U), !is.null(gss_basis$varnames))
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Need data.table.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Need ggplot2.")
  if (!requireNamespace("scales", quietly = TRUE)) stop("Need scales.")
  
  order_vars_by <- match.arg(order_vars_by)
  order_pcs_by  <- match.arg(order_pcs_by)
  
  U <- gss_basis$U
  pcs <- seq_len(ncol(U))
  colnames(U) <- paste0("PC", pcs)
  
  dt <- data.table::as.data.table(U)
  dt[, var := gss_basis$varnames]
  
  dt_long <- data.table::melt(
    dt,
    id.vars = "var",
    variable.name = "pc",
    value.name = "loading"
  )
  
  # Color by absolute magnitude, but keep signed values for labels
  dt_long[, abs_loading := abs(loading)]
  
  # Variable ordering
  if (order_vars_by == "max_abs") {
    vr <- dt_long[, .(max_abs = max(abs_loading, na.rm = TRUE)), by = var][order(-max_abs)]
    dt_long[, var := factor(var, levels = rev(vr$var))]  # strongest at top
  } else {
    dt_long[, var := factor(var, levels = rev(gss_basis$varnames))]
  }
  
  # PC ordering
  if (order_pcs_by == "original") {
    dt_long[, pc := factor(pc, levels = paste0("PC", pcs))]
  } else {
    # eigen() output order (usually decreasing eigenvalue); already PC1..PCn in that order
    dt_long[, pc := factor(pc, levels = paste0("PC", pcs))]
  }
  
  # Labels are the ORIGINAL signed loadings
  fmt <- function(z) formatC(z, format = "f", digits = value_digits)
  dt_long[, label := fmt(loading)]
  
  n_vars <- length(gss_basis$varnames)
  n_pcs  <- ncol(gss_basis$U)
  text_size <- max(0.25, min(2.2, 240 / (n_vars + n_pcs)))  # always draws, but shrinks
  
  ggplot2::ggplot(dt_long, ggplot2::aes(x = pc, y = var)) +
    ggplot2::geom_tile(ggplot2::aes(fill = abs_loading), color = "grey90", linewidth = 0.10) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = text_size) +
    ggplot2::scale_fill_gradient(
      low = "white",
      high = "orangered",
      limits = c(0, 1),         # abs(loadings) are bounded by 1 for eigenvectors
      oob = scales::squish,
      name = "|loading|"
    ) +
    ggplot2::labs(
      title = "GSS PC loadings",
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid = ggplot2::element_blank()
    )
}

# --- Call at end of your script ---
p_loadings_full <- plot_pc_loadings_grid_full_signed_labels(gss_basis, value_digits = 2)

OUTFILE_LOADINGS_PDF <- file.path(
  BASE_VIZ_DIR,
  sprintf("missing_variance_rayleigh_%s_loadings_grid.pdf", YEAR)
)

# Size PDF to grid so labels actually exist on the page
n_vars <- length(gss_basis$varnames)
n_pcs  <- ncol(gss_basis$U)
pdf_w <- min(120, max(16, 0.22 * n_pcs + 4))
pdf_h <- min(120, max(10, 0.15 * n_vars + 4))

grDevices::pdf(OUTFILE_LOADINGS_PDF, width = pdf_w, height = pdf_h, onefile = TRUE)
print(p_loadings_full)
grDevices::dev.off()

message("Wrote FULL loadings grid plot to: ", OUTFILE_LOADINGS_PDF)
message(sprintf("PDF size: %.1f x %.1f inches (n_vars=%d, n_pcs=%d).", pdf_w, pdf_h, n_vars, n_pcs))
