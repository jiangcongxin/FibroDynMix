read_selected_nb_defaults <- function(path = file.path("analysis", "validation_aware_nb_selection", "validation_aware_nb_selected_summary.tsv")) {
  if (!file.exists(path) || file.info(path)$size == 0L) {
    stop(sprintf("Missing selected NB defaults table: %s", path), call. = FALSE)
  }
  utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}

selected_nb_default <- function(defaults, scenario, variant = "nb", fallback = 2L) {
  required <- c("scenario", "variant", "selected_n_outer")
  missing <- setdiff(required, colnames(defaults))
  if (length(missing) > 0L) {
    stop(sprintf("Selected NB defaults table is missing columns: %s", paste(missing, collapse = ", ")), call. = FALSE)
  }
  hit <- defaults$selected_n_outer[defaults$scenario == scenario & defaults$variant == variant]
  hit <- suppressWarnings(as.integer(hit))
  hit <- hit[!is.na(hit) & hit > 0L]
  if (length(hit) == 0L) {
    return(as.integer(fallback))
  }
  hit[[1L]]
}

selected_nb_defaults_used <- function(defaults) {
  data.frame(
    task = c(
      "vi_continuous",
      "vi_batch_confounding",
      "vi_rare_transition",
      "cross_cohort_transfer_batch",
      "realdata_multistudy_downstream"
    ),
    scenario = c(
      "continuous",
      "batch_confounding",
      "rare_transition",
      "batch_confounding",
      "batch_confounding"
    ),
    variant = c("nb", "nb", "nb", "nb", "nb_study"),
    selected_n_outer = c(
      selected_nb_default(defaults, "continuous", "nb", 2L),
      selected_nb_default(defaults, "batch_confounding", "nb", 10L),
      selected_nb_default(defaults, "rare_transition", "nb", 2L),
      selected_nb_default(defaults, "batch_confounding", "nb", 10L),
      selected_nb_default(defaults, "batch_confounding", "nb_study", 2L)
    ),
    selection_source = "analysis/validation_aware_nb_selection/validation_aware_nb_selected_summary.tsv",
    stringsAsFactors = FALSE
  )
}
