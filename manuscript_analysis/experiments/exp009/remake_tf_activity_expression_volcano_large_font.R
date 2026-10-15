#!/usr/bin/env Rscript
# High-font variant of the standard EXP009 volcano plots. Results go in a
# separate directory so the standard-size figures remain untouched.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Cannot locate this script's path")
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]))
script_dir <- dirname(script_path)
project_root <- "/project/zduren/durenlab/naqing_palmetto/Alzheimer_project_explore"

Sys.setenv(
  EXP009_OUTPUT_ROOT = file.path(project_root, "experiments/EXP009_results/large_font"),
  EXP009_FIGURE_WIDTH = "11",
  EXP009_FIGURE_HEIGHT = "8",
  EXP009_PNG_DPI = "300",
  EXP009_AXIS_TEXT_SIZE = "25",
  EXP009_AXIS_TITLE_SIZE = "30",
  EXP009_PLOT_TITLE_SIZE = "30",
  EXP009_TF_LABEL_SIZE = "4.8",
  EXP009_LEGEND_TEXT_SIZE = "18",
  EXP009_CAPTION_SIZE = "15",
  EXP009_EMPTY_LABEL_SIZE = "7"
)

sys.source(file.path(script_dir, "remake_tf_activity_expression_volcano.R"),
           envir = globalenv())
