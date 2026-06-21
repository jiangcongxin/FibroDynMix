#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else "scripts/make_figure2.R"
ROOT <- normalizePath(file.path(dirname(script_path), ".."))

Sys.setenv(
  FDM_FIGURE_OUT = file.path(ROOT, "figures", "figure2"),
  FDM_FIGURE_STEM = "figure2",
  FDM_FIGURE_LABEL = "Figure 2"
)

source(file.path(ROOT, "scripts", "make_figure2_redesign.R"))
