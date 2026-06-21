fdm_write_tsv <- function(x, path) {
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

fdm_read_tsv <- function(path) {
  utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}

fdm_wrap_label <- function(x, width = 18) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1))
}

fdm_state_palette <- function() {
  c(
    resident = "#4E79A7",
    inflammatory = "#E15759",
    myofibroblast = "#59A14F",
    "ECM-remodeling" = "#B07AA1",
    "antigen-presenting" = "#F28E2B",
    "IFN-stress" = "#76B7B2"
  )
}

fdm_condition_palette <- function() {
  c(
    normal = "#4E79A7",
    disease = "#E15759",
    control = "#6B7280",
    perturbed = "#E15759"
  )
}

fdm_method_palette <- function() {
  c(
    marker_scoring = "#6B7280",
    fibrodynmix_initializer = "#1B998B",
    fibrodynmix_nb = "#386CB0",
    fibrodynmix_vi = "#984EA3"
  )
}

fdm_evidence_palette <- function() {
  c(
    model = "#E8EEF6",
    simulation = "#E2F0D9",
    public_data = "#FDE7C7",
    uncertainty = "#EADCF8",
    boundary = "#ECECEC"
  )
}

fdm_theme <- function(base_size = 8, base_family = "Helvetica") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = base_size, colour = "#222222"),
      axis.text = ggplot2::element_text(size = base_size - 1, colour = "#222222"),
      axis.line = ggplot2::element_line(linewidth = 0.25, colour = "#333333"),
      axis.ticks = ggplot2::element_line(linewidth = 0.25, colour = "#333333"),
      plot.title = ggplot2::element_text(size = base_size + 1, face = "bold", hjust = 0, colour = "#111111"),
      plot.subtitle = ggplot2::element_text(size = base_size - 0.5, hjust = 0, colour = "#555555"),
      plot.caption = ggplot2::element_text(size = base_size - 1.5, hjust = 0, colour = "#666666"),
      legend.title = ggplot2::element_text(size = base_size - 1),
      legend.text = ggplot2::element_text(size = base_size - 1),
      legend.key.height = grid::unit(0.13, "in"),
      legend.key.width = grid::unit(0.18, "in"),
      legend.spacing.x = grid::unit(0.04, "in"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = base_size, face = "bold", colour = "#222222"),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

fdm_blank_theme <- function(base_size = 8, base_family = "Helvetica") {
  ggplot2::theme_void(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = base_size + 1, face = "bold", hjust = 0, colour = "#111111"),
      plot.subtitle = ggplot2::element_text(size = base_size - 0.5, hjust = 0, colour = "#555555"),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.title = ggplot2::element_text(size = base_size - 1),
      legend.text = ggplot2::element_text(size = base_size - 1)
    )
}

fdm_panel_tag_theme <- function() {
  ggplot2::theme(
    plot.tag = ggplot2::element_text(size = 13, face = "bold", colour = "#111111"),
    plot.tag.position = c(0, 1)
  )
}

fdm_save_plot <- function(plot, exports_dir, stem, width, height, dpi = 450) {
  ggplot2::ggsave(file.path(exports_dir, paste0(stem, ".pdf")), plot, width = width, height = height, device = grDevices::pdf)
  ggplot2::ggsave(file.path(exports_dir, paste0(stem, ".svg")), plot, width = width, height = height, device = svglite::svglite)
  ggplot2::ggsave(file.path(exports_dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = dpi, device = ragg::agg_png)
  ggplot2::ggsave(file.path(exports_dir, paste0(stem, ".tiff")), plot, width = width, height = height, dpi = dpi, device = "tiff", compression = "lzw")
}

fdm_write_export_qc <- function(exports_dir, qc_dir, stem) {
  if (requireNamespace("magick", quietly = TRUE)) {
    img <- magick::image_read(file.path(exports_dir, paste0(stem, ".png")))
    magick::image_write(img, file.path(exports_dir, "contact_sheet.png"))
  }

  qc_files <- file.path(exports_dir, c(
    paste0(stem, ".pdf"),
    paste0(stem, ".svg"),
    paste0(stem, ".png"),
    paste0(stem, ".tiff"),
    "contact_sheet.png"
  ))
  qc <- data.frame(
    file = basename(qc_files),
    exists = file.exists(qc_files),
    bytes = ifelse(file.exists(qc_files), file.info(qc_files)$size, NA_real_),
    stringsAsFactors = FALSE
  )
  if (requireNamespace("magick", quietly = TRUE)) {
    dims <- lapply(qc_files, function(path) {
      if (!file.exists(path) || !grepl("\\.(png|tiff)$", path)) {
        return(c(width = NA_real_, height = NA_real_))
      }
      info <- magick::image_info(magick::image_read(path))
      c(width = info$width[1], height = info$height[1])
    })
    dims <- do.call(rbind, dims)
    qc$width <- dims[, "width"]
    qc$height <- dims[, "height"]
  }
  fdm_write_tsv(qc, file.path(qc_dir, "export_image_qc.tsv"))
}
