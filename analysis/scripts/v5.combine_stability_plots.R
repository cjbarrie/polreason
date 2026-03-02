################################################################################
# Combine corr_mice_stability_check.pdf plots from all models into one figure
################################################################################

suppressPackageStartupMessages({
  library(magick)
  library(grid)
  library(gridExtra)
})

combine_stability_plots <- function(
    base_viz_dir    = "analysis/viz",
    year            = 2024,
    output_dir      = NULL,
    ncol            = 2,
    nrow_per_page   = 3,
    scale           = 0.5,
    density         = 300   # DPI for reading/writing PDFs
) {

  message("\n=== Combining Stability Check Plots ===")

  # Find all model directories with the stability check PDF
  model_dirs <- list.dirs(base_viz_dir, recursive = FALSE, full.names = TRUE)
  model_dirs <- model_dirs[grepl(sprintf("-%d$", year), model_dirs)]

  # Find PDFs
  pdf_files <- file.path(model_dirs, "corr_mice_stability_check.pdf")
  pdf_files <- pdf_files[file.exists(pdf_files)]

  message("Found ", length(pdf_files), " stability check PDFs")

  if (length(pdf_files) == 0) {
    stop("No corr_mice_stability_check.pdf files found")
  }

  # Extract model names from paths
  model_names <- basename(dirname(pdf_files))
  model_names <- sub(sprintf("-%d$", year), "", model_names)

  # Read all PDFs as images
  message("Reading PDFs at ", density, " DPI...")
  images <- lapply(seq_along(pdf_files), function(i) {
    message("  ", i, "/", length(pdf_files), ": ", model_names[i])
    img <- image_read_pdf(pdf_files[i], density = density)
    img
  })

  # Add model name as label to each image
  message("\nAdding labels...")
  labeled_images <- lapply(seq_along(images), function(i) {
    img <- images[[i]]
    label <- model_names[i]

    img_labeled <- image_annotate(
      img,
      label,
      size = 24,
      gravity = "north",
      color = "black",
      boxcolor = "white",
      font = "Helvetica"
    )

    img_labeled
  })

  # Resize images
  resized_images <- lapply(labeled_images, function(img) {
    info <- image_info(img)
    new_width <- round(info$width * scale)
    image_resize(img, geometry = sprintf("%dx", new_width))
  })

  # Calculate pagination
  models_per_page <- ncol * nrow_per_page
  n_images <- length(resized_images)
  n_pages <- ceiling(n_images / models_per_page)

  message("\nCreating ", n_pages, " pages with ", ncol, " cols x ", nrow_per_page, " rows each")

  # Output directory
  if (is.null(output_dir)) {
    output_dir <- base_viz_dir
  }

  # Create blank image for padding
  info <- image_info(resized_images[[1]])
  blank <- image_blank(info$width, info$height, color = "white")

  # Generate each page
  for (page in seq_len(n_pages)) {
    start_idx <- (page - 1) * models_per_page + 1
    end_idx <- min(page * models_per_page, n_images)

    page_images <- resized_images[start_idx:end_idx]

    # Pad to fill page if needed
    n_pad <- models_per_page - length(page_images)
    if (n_pad > 0) {
      for (i in seq_len(n_pad)) {
        page_images <- c(page_images, list(blank))
      }
    }

    # Combine into rows
    rows <- list()
    for (r in seq_len(nrow_per_page)) {
      row_start <- (r - 1) * ncol + 1
      row_end <- r * ncol
      row_images <- page_images[row_start:row_end]

      row_combined <- image_append(image_join(row_images), stack = FALSE)
      rows[[r]] <- row_combined
    }

    # Combine rows vertically
    page_image <- image_append(image_join(rows), stack = TRUE)

    # Save page
    out_file <- file.path(output_dir, sprintf("stability_check_%d_%d.pdf", year, page))

    message("Saving page ", page, "/", n_pages, ": ", out_file,
            " (models ", start_idx, "-", end_idx, ")")

    image_write(page_image, out_file, format = "pdf", density = density)
  }

  # Report dimensions
  final_info <- image_info(page_image)
  message("\nDimensions per page: ", final_info$width, " x ", final_info$height, " pixels")

  message("\n=== Done ===\n")

  invisible(NULL)
}

################################################################################
# MAIN EXECUTION
################################################################################

if (exists("BASE_VIZ_DIR") && exists("YEAR")) {
  combine_stability_plots(
    base_viz_dir  = BASE_VIZ_DIR,
    year          = YEAR,
    ncol          = 2,
    nrow_per_page = 3,
    scale         = 0.5
  )
} else {
  # Run with defaults: 2 cols x 3 rows = 6 models per page
  combine_stability_plots(
    base_viz_dir  = "analysis/viz",
    year          = 2024,
    ncol          = 2,
    nrow_per_page = 3,
    scale         = 0.5
  )
}
