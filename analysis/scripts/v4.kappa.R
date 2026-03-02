################################################################################
# Pairwise agreement heatmap: Cohen's kappa (point estimates only)
# (lower triangle only) — GSS forced to be first in the matrix ordering
#
# Overview
# --------
# This script measures and visualizes pairwise agreement between raters
# (GSS + multiple LLM raters) using Cohen’s kappa point estimates.
#
# We:
#   1) Discover per-rater harmonized data files on disk for YEAR.
#   2) For each rater:
#        - keep non-persona belief variables only (persona_var == FALSE),
#        - collapse multiple runs into a single “modal” answer per
#          (persona_id, variable) by:
#             a) taking the first non-missing answer within each run, then
#             b) taking the random mode across runs (mode_random_ties()).
#   3) Merge all raters into a wide table keyed by (persona_id, variable),
#      keeping GSS as the first rater column (robust to naming differences).
#   4) Construct an additional synthetic rater:
#        - answer.llm_agg = modal response across *all LLM raters* (excluding GSS)
#      and place it immediately after GSS in the matrix ordering.
#   5) Compute a kappa matrix of point estimates (no CIs):
#        - weighted Cohen’s kappa via irr::kappa2(weight=WEIGHT),
#        - fixed category levels per pair using the union of observed responses,
#        - store only the lower triangle (i > j) to avoid redundancy.
#   6) Render a grid-based heatmap to PDF:
#        - diagonal shaded,
#        - upper triangle blank,
#        - cell colors from col_kappa(k) with a -1..1 legend.
#
# No-confidence-interval policy
# -----------------------------
# We do not compute CIs for kappa here. With large n, CI widths tend to shrink
# and are not the limiting factor for interpretation in this workflow.
#
#
# Persona handling
# ---------------
# We exclude persona variables by default:
#   - if persona_var exists, we filter persona_var == FALSE.
# This keeps agreement comparisons focused on belief items rather than persona
# descriptors, matching the harmonized-data convention used elsewhere.
#
# Plot styling notes
# ------------------
# The heatmap is drawn with grid graphics (not ggplot2) to allow:
#   - large margins for long rater labels,
#   - rotation and manual positioning of axis labels,
#   - a compact legend inset in the top-right of the panel.
#
# Palette consistency (heatmap ↔ ggplot)
# --------------------------------------
# This file also defines a ggplot2 fill scale that matches the heatmap palette
# exactly by sampling col_kappa() over a dense grid in [-1, 1] and using
# scale_fill_gradientn(..., values=rescale(k_grid)).
#
# Optional downstream plot: constraints vs agreement
# --------------------------------------------------
# If constraint_point_stats_YEAR.csv is present, we:
#   - merge each LLM’s kappa-vs-GSS with constraint diagnostics (overall group),
#   - run OLS: kappa_vs_gss ~ pc1_var_explained + effective_dependence,
#   - plot (PC1 share, De) with squares filled by kappa-vs-GSS using the same
#     col_kappa() palette and annotate regression coefficients + R².
#
# Prerequisites (provided elsewhere in the project)
# ------------------------------------------------
# Required globals:
#   - BASE_OUT_DIR
#   - BASE_VIZ_DIR
#
# Required packages
# -----------------
#   - data.table
#   - grid
#   - grDevices
#   - irr
#   - ggplot2 (only for the constraints-vs-agreement scatter)
#   - scales   (only for the constraints-vs-agreement scatter)
#
# Output
# ------
# Writes a PDF heatmap to:
#   OUTFILE_PDF = file.path(BASE_VIZ_DIR, sprintf("pairwise_kappa_%s.pdf", YEAR))
#
# Optionally writes a 3-way scatter PDF to:
#   OUTFILE_PDF_3WAY = file.path(BASE_VIZ_DIR, sprintf("constraints_pc1_De_kappa_%s.pdf", YEAR))
################################################################################


################################################################################
# Pairwise agreement heatmap: Cohen's kappa (point estimates only)
# (lower triangle only) — GSS forced to be first in the matrix ordering
# NOTE: We do not compute confidence intervals; with large n, CI width is expected
#       to shrink toward 0 and is not limiting for interpretation here.
################################################################################

YEAR        <- get0("YEAR", ifnotfound = 2024L)
RATER_GSS   <- "gss"

OUTFILE_PDF <- file.path(BASE_VIZ_DIR, sprintf("pairwise_kappa_%s.pdf", YEAR))
PDF_WIDTH   <- 12.5
PDF_HEIGHT  <- 12.5

WEIGHT      <- "squared"

# ---- toggle: show synthetic LLM aggregate in heatmap ------------------------
# If FALSE: we can still (optionally) compute answer.llm_agg, but we do NOT
# include it in the kappa matrix/heatmap ordering.
SHOW_LLM_AGG <- get0("SHOW_LLM_AGG", ifnotfound = FALSE)

# If you want to compute it even when not showing it (e.g., for other downstream uses),
# leave this as-is. If you want "compute only when showing", set:
# CREATE_LLM_AGG <- SHOW_LLM_AGG
CREATE_LLM_AGG <- get0("CREATE_LLM_AGG", ifnotfound = FALSE)


CELL_TEXT_FONTSIZE <- 9.6
AXIS_FONTSIZE      <- 9.5
TITLE_FONTSIZE     <- 15
NOTE_FONTSIZE      <- 10

## Increased margins
LEFT_MARGIN_MM   <- 75
BOTTOM_MARGIN_MM <- 80
TOP_MARGIN_MM    <- 18
RIGHT_MARGIN_MM  <- 10

# ---- helpers ---------------------------------------------------------------

rater_from_path <- function(f) basename(dirname(f))

mode_random_ties <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_character_)
  
  tab <- table(x)
  m   <- max(tab)
  top <- names(tab)[tab == m]
  
  if (length(top) == 1L) top else sample(top, 1L)
}


clean_rater_label <- function(x) {
  x <- gsub(as.character(YEAR), "", x, fixed = TRUE)
  x <- gsub("^[. _-]+", "", x, perl = TRUE)
  x <- gsub("[. _-]+$", "", x, perl = TRUE)
  x <- gsub("[. _-]{2,}", "_", x, perl = TRUE)
  x
}

# ---- build ratings_wide_modal ---------------------------------------------

required_objects <- c("BASE_OUT_DIR", "BASE_VIZ_DIR")
missing_objects  <- required_objects[!vapply(required_objects, exists, logical(1), inherits = TRUE)]
if (length(missing_objects)) {
  stop("Missing required objects in the workspace: ", paste(missing_objects, collapse = ", "))
}

all_rds <- list.files(
  path = BASE_OUT_DIR,
  pattern = "\\.rds$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

rds_files <- all_rds[
  grepl("^harmoni(s|z)ed_data\\.rds$", basename(all_rds), ignore.case = TRUE)
]

extract_rater_modal_dt <- function(f, keep_persona_vars = FALSE) {
  dt <- as.data.table(readRDS(f))
  
  if (!keep_persona_vars && "persona_var" %chin% names(dt)) dt <- dt[persona_var == FALSE]
  if (!("run" %chin% names(dt))) dt[, run := 1L]
  
  dt <- dt[, .(persona_id, variable, run, answer = as.character(answer))]
  setorder(dt, persona_id, variable, run)
  
  # first per run, then modal across runs
  dt <- dt[, .(answer = {
    a <- answer[!is.na(answer)]
    if (length(a)) a[1L] else NA_character_
  }), by = .(persona_id, variable, run)][
    , .(answer = mode_random_ties(answer)), by = .(persona_id, variable)
  ]
  
  setnames(dt, "answer", make.names(paste0("answer.", rater_from_path(f))))
  dt
}

dt_list <- lapply(rds_files, extract_rater_modal_dt, keep_persona_vars = FALSE)

# put GSS first in merge list (nice-to-have; matrix order enforced below)
rater_names <- sapply(rds_files, rater_from_path)
is_gss      <- tolower(rater_names) == tolower(RATER_GSS)
dt_list     <- c(dt_list[is_gss], dt_list[!is_gss])

keys <- c("persona_id", "variable")
ratings_wide_modal <- Reduce(function(x, y) merge(x, y, by = keys, all.x = TRUE), dt_list)

# ---- enforce column order: GSS FIRST (robust to naming differences) ---------

ans_cols_all <- grep("^answer\\.", names(ratings_wide_modal), value = TRUE)

raters_raw   <- sub("^answer\\.", "", ans_cols_all)
raters_clean <- tolower(vapply(raters_raw, clean_rater_label, character(1)))

gss_pos <- which(raters_clean == tolower(RATER_GSS))
gss_col <- if (length(gss_pos)) ans_cols_all[gss_pos[1L]] else character(0)

ans_cols <- if (length(gss_col)) c(gss_col, setdiff(ans_cols_all, gss_col)) else ans_cols_all

# ---- LLM columns (exclude GSS) ---------------------------------------------

llm_cols <- if (length(gss_col)) setdiff(ans_cols_all, gss_col) else ans_cols_all

# ---- (optional) aggregate LLM prediction -----------------------------------

if (CREATE_LLM_AGG) {
  if (length(llm_cols) == 0L) {
    warning("CREATE_LLM_AGG=TRUE but found 0 LLM columns (nothing to aggregate).")
  } else {
    ratings_wide_modal[, answer.llm_agg := {
      v <- unlist(.SD, use.names = FALSE)
      mode_random_ties(v)
    }, by = .(persona_id, variable), .SDcols = llm_cols]
  }
}

# ---- enforce matrix order: GSS first, optional llm_agg second --------------

if (length(gss_col)) {
  ans_cols <- c(gss_col, setdiff(ans_cols_all, gss_col))
} else {
  ans_cols <- ans_cols_all
}

if (SHOW_LLM_AGG && ("answer.llm_agg" %in% names(ratings_wide_modal))) {
  ans_cols <- if (length(gss_col)) {
    c(gss_col, "answer.llm_agg", setdiff(ans_cols_all, gss_col))
  } else {
    c("answer.llm_agg", ans_cols_all)
  }
}

message("SHOW_LLM_AGG = ", SHOW_LLM_AGG)
message("Matrix order: ", paste(ans_cols, collapse = " | "))


# ---- kappa point estimate for one pair (with fixed levels) -----------------

kappa_pair <- function(dt, col_x, col_y, lev, weight = "squared") {
  x <- dt[[col_x]]
  y <- dt[[col_y]]
  keep <- !is.na(x) & !is.na(y)
  if (sum(keep) < 2L) return(NA_real_)
  
  xx <- factor(x[keep], levels = lev, ordered = TRUE)
  yy <- factor(y[keep], levels = lev, ordered = TRUE)
  
  ik <- tryCatch(
    irr::kappa2(cbind(xx, yy), weight = weight),
    error = function(e) NULL
  )
  if (is.null(ik)) return(NA_real_)
  as.numeric(ik$value)
}

# ---- kappa matrix: point estimates only (lower triangle) --------------------

kappa_mat_point <- function(wide_dt, cols, weight = "squared") {
  n <- length(cols)
  K0 <- matrix(NA_real_, n, n, dimnames = list(cols, cols))
  
  for (i in 2:n) for (j in 1:(i - 1)) {
    x <- wide_dt[[cols[i]]]
    y <- wide_dt[[cols[j]]]
    keep <- !is.na(x) & !is.na(y)
    if (sum(keep) < 2L) next
    
    lev <- sort(unique(c(x[keep], y[keep])))
    K0[i, j] <- kappa_pair(wide_dt, cols[i], cols[j], lev, weight = weight)
  }
  
  K0
}

# ---- plot helpers ----------------------------------------------------------

col_kappa <- function(k) {
  if (is.na(k)) return("grey95")
  k <- max(-1, min(1, k))
  r1 <- grDevices::colorRamp(c("purple", "orangered"))
  r2 <- grDevices::colorRamp(c("orangered", "springgreen"))
  grDevices::rgb((if (k <= 0) r1(k + 1) else r2(k)) / 255)
}

plot_kappa <- function(K,
                       title,
                       note = "",
                       cell_fontsize = 9.6,
                       axis_fontsize = 9.5,
                       title_fontsize = 15,
                       note_fontsize = 10) {
  n <- nrow(K)
  labs <- vapply(sub("^answer\\.", "", rownames(K)), clean_rater_label, character(1))
  
  grid.newpage()
  
  grid.text(
    title,
    x = unit(0, "npc"), y = unit(1, "npc") - unit(2, "mm"),
    just = c("left", "top"),
    gp = gpar(fontsize = title_fontsize, fontface = "bold")
  )
  
  grid.text(
    note,
    x = unit(0, "npc"), y = unit(1, "npc") - unit(10, "mm"),
    just = c("left", "top"),
    gp = gpar(fontsize = note_fontsize)
  )
  
  ml <- unit(LEFT_MARGIN_MM, "mm")
  mb <- unit(BOTTOM_MARGIN_MM, "mm")
  mt <- unit(TOP_MARGIN_MM, "mm")
  mr <- unit(RIGHT_MARGIN_MM, "mm")
  
  pushViewport(viewport(
    x = ml, y = mb,
    width  = unit(1, "npc") - ml - mr,
    height = unit(1, "npc") - mb - mt,
    just = c("left", "bottom"),
    xscale = c(0, n), yscale = c(0, n)
  ))
  
  for (i in seq_len(n)) for (j in seq_len(n)) {
    x <- j - .5
    y <- n - i + .5
    
    if (i > j) {
      kv <- K[i, j]
      grid.rect(x, y, 1, 1, default.units = "native",
                gp = gpar(fill = col_kappa(kv), col = "grey85", lwd = .55))
      
      if (is.finite(kv)) {
        grid.text(sprintf("%.2f", kv), x, y, default.units = "native",
                  gp = gpar(fontsize = cell_fontsize, fontface = "bold"))
      }
    } else if (i == j) {
      grid.rect(x, y, 1, 1, default.units = "native",
                gp = gpar(fill = "grey92", col = "grey85", lwd = .55))
    } else {
      grid.rect(x, y, 1, 1, default.units = "native",
                gp = gpar(fill = "white", col = NA))
    }
  }
  
  # Y labels
  for (i in seq_len(n)) {
    grid.text(labs[i],
              x = unit(-2, "mm"),
              y = n - i + .5,
              default.units = "native",
              just = "right",
              gp = gpar(fontsize = axis_fontsize))
  }
  
  # X labels
  y_lab <- unit(0, "native") - unit(1.5, "mm")
  for (i in seq_len(n)) {
    grid.text(labs[i],
              x = i - .5, y = y_lab,
              default.units = "native",
              rot = 90, just = "right",
              gp = gpar(fontsize = axis_fontsize))
  }
  
  # Legend inside top-right
  pushViewport(viewport(
    x = unit(0.985, "npc"),
    y = unit(0.985, "npc"),
    width  = unit(22, "mm"),
    height = unit(75, "mm"),
    just = c("right", "top")
  ))
  
  grid.text("kappa", x = unit(0.5, "npc"), y = unit(1, "npc") - unit(1, "mm"),
            just = c("center", "top"),
            gp = gpar(fontsize = 11, fontface = "bold"))
  
  steps <- 140
  bar_w <- unit(5.0, "mm")
  bar_h <- unit(58, "mm")
  xmid  <- unit(0.35, "npc")
  ytop  <- unit(1, "npc") - unit(10, "mm")
  
  for (s in seq_len(steps)) {
    t <- (s - 1) / (steps - 1)
    k <- 1 - 2 * t
    grid.rect(xmid, ytop - bar_h * t, bar_w, bar_h / steps,
              just = "center",
              gp = gpar(fill = col_kappa(k), col = NA))
  }
  
  tx <- xmid + unit(7.5, "mm")
  grid.text("1",  tx, ytop,           just = c("left","center"), gp = gpar(fontsize = 9))
  grid.text("0",  tx, ytop - bar_h/2, just = c("left","center"), gp = gpar(fontsize = 9))
  grid.text("-1", tx, ytop - bar_h,   just = c("left","center"), gp = gpar(fontsize = 9))
  
  upViewport(2)
}

# ---- run + save ------------------------------------------------------------

K <- kappa_mat_point(
  ratings_wide_modal,
  ans_cols,
  weight = WEIGHT
)

pdf(OUTFILE_PDF, width = PDF_WIDTH, height = PDF_HEIGHT)
plot_kappa(
  K,
  title = sprintf("Pairwise agreement: Cohen's kappa (weight=%s)", WEIGHT),
  cell_fontsize = CELL_TEXT_FONTSIZE,
  axis_fontsize = AXIS_FONTSIZE,
  title_fontsize = TITLE_FONTSIZE,
  note_fontsize = NOTE_FONTSIZE
)
dev.off()

message("Wrote kappa plot to: ", OUTFILE_PDF)




# --- make ggplot fill scale match the heatmap's col_kappa() exactly ---------

# NOTE: col_kappa() includes dodgerblue between orangered and springgreen (positive side).
col_kappa <- function(k) {
  if (is.na(k)) return("grey95")
  k <- max(-1, min(1, k))
  
  rneg <- grDevices::colorRamp(c("purple", "orangered"))
  rpos <- grDevices::colorRamp(c("orangered", "dodgerblue", "springgreen"))
  
  rgb <- if (k <= 0) rneg(k + 1) else rpos(k)
  grDevices::rgb(rgb / 255)
}

col_kappa_vec <- Vectorize(col_kappa, USE.NAMES = FALSE)

k_grid   <- seq(-1, 1, length.out = 401)
pal_cols <- col_kappa_vec(k_grid)


# point estimates from script 1
constraint_file <- file.path(BASE_OUT_DIR, sprintf("constraint_point_stats_%s.csv", YEAR))
if (!file.exists(constraint_file)) {
  warning("Missing constraint point-estimate file: ", constraint_file,
          "\nRun the constraint diagnostics script first.")
} else {
  
  constraint_point <- data.table::fread(constraint_file)
  
  # enforce merge key (in case you forgot to add it in script 1)
  if (!("rater_clean" %in% names(constraint_point))) {
    constraint_point[, rater_clean := tolower(clean_rater_label(rater))]
  }
  
  # use overall by default (change if you want facets by edu_group)
  constraint_point <- constraint_point[edu_group == "overall"]
  
  # build kappa-vs-gss table from K
  rn <- rownames(K)
  cn <- colnames(K)
  
  gss_idx <- if (length(gss_col)) match(gss_col, cn) else NA_integer_
  if (!is.finite(gss_idx)) {
    warning("Could not find GSS column in kappa matrix; cannot compute kappa-vs-gss.")
  } else {
    
    kappa_vs_gss <- data.table::data.table(
      ans_col      = rn,
      rater_raw    = sub("^answer\\.", "", rn),
      rater_clean  = tolower(vapply(sub("^answer\\.", "", rn), clean_rater_label, character(1))),
      kappa_vs_gss = as.numeric(K[, gss_idx])
    )
    
    # keep LLMs only (exclude gss; optionally exclude llm_agg)
    kappa_vs_gss <- kappa_vs_gss[
      rater_clean != tolower(RATER_GSS) &
        rater_clean != "llm_agg"
    ]
    
    # merge constraints + kappa
    constraint_point[, rater_clean := gsub("-", ".", rater_clean)]
    plot_dt <- merge(
      constraint_point,
      kappa_vs_gss[, .(rater_clean, kappa_vs_gss)],
      by  = "rater_clean",
      all = FALSE
    )
    
    plot_dt <- plot_dt[
      is.finite(kappa_vs_gss) &
        is.finite(pc1_var_explained) &
        is.finite(effective_dependence)
    ]
    
    # ---- linear regression: agreement ~ constraints ------------------------
    
    fit <- stats::lm(kappa_vs_gss ~ pc1_var_explained + effective_dependence, data = plot_dt)
    cf  <- summary(fit)$coefficients
    
    b1 <- unname(cf["pc1_var_explained",    "Estimate"])
    p1 <- unname(cf["pc1_var_explained",    "Pr(>|t|)"])
    
    b2 <- unname(cf["effective_dependence", "Estimate"])
    p2 <- unname(cf["effective_dependence", "Pr(>|t|)"])
    
    r2 <- summary(fit)$r.squared
    
    fmt_p <- function(p) {
      if (!is.finite(p)) return("NA")
      if (p < 1e-3) "<0.001" else formatC(p, format = "f", digits = 3)
    }
    
    coef_lab <- sprintf(
      "OLS: kappa ~ PC1 + De\nb_PC1 = %.3f (p=%s)\nb_De  = %.3f (p=%s)\nR² = %.2f",
      b1, fmt_p(p1), b2, fmt_p(p2), r2
    )
    
    # OLS label positioned bottom-right *inside* panel
    ann_dt <- data.table::data.table(x = 0.98, y = 0.02, lab = coef_lab)
    
    # output
    OUTFILE_PDF_3WAY <- file.path(BASE_VIZ_DIR, sprintf("constraints_pc1_De_kappa_%s.pdf", YEAR))
    
    p <- ggplot2::ggplot(plot_dt, ggplot2::aes(x = pc1_var_explained, y = effective_dependence)) +
      ggplot2::geom_point(
        ggplot2::aes(fill = kappa_vs_gss),
        shape  = 22,
        size   = 2,
        colour = "black",
        stroke = 0.25,
        alpha  = 0.5
      ) +
      ggplot2::scale_x_continuous(
        limits = c(0, 1),
        breaks = scales::pretty_breaks(5),
        labels = scales::label_number(accuracy = 0.01),
        name   = "PC1 share of variance"
      ) +
      ggplot2::scale_y_continuous(
        limits = c(0, 1),
        breaks = scales::pretty_breaks(5),
        labels = scales::label_number(accuracy = 0.01),
        name   = "Effective dependence De (1 - |R|^(1/p))"
      ) +
      ggplot2::scale_fill_gradientn(
        colours = pal_cols,
        values  = scales::rescale(k_grid, to = c(0, 1), from = c(-1, 1)),
        limits  = c(-1, 1),
        breaks  = c(-1, -0.5, 0, 0.5, 1),
        name    = "Cohen's kappa\n(vs GSS)"
      ) +
      ggplot2::ggtitle("Constraints vs Agreement (overall)") +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        legend.position   = "right",
        legend.title      = ggplot2::element_text(size = 8),
        legend.text       = ggplot2::element_text(size = 7),
        legend.key.height = grid::unit(3.0, "mm"),
        legend.key.width  = grid::unit(3.0, "mm")
      ) +
      ggplot2::guides(fill = ggplot2::guide_colorbar(
        barheight = grid::unit(40, "mm"),
        barwidth  = grid::unit(5, "mm"),
        title.position = "top",
        title.hjust = 0.5
      )) +
      ggplot2::geom_label(
        data = ann_dt,
        ggplot2::aes(x = x, y = y, label = lab),
        inherit.aes = FALSE,
        hjust = 1, vjust = 0,
        size = 3.2,
        linewidth = 0.25,
        fill = "white",
        alpha = 0.85
      )
    
    ggplot2::ggsave(OUTFILE_PDF_3WAY, plot = p, width = 5, height = 4, units = "in")
    message("Wrote 3-way plot to: ", OUTFILE_PDF_3WAY)
  }
}



