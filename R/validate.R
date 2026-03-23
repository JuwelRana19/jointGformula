validate_gformula_data <- function(data, config) {
  dt <- as.data.table(copy(data))

  required_cols <- unique(c(
    config$id,
    config$time_in,
    config$time_out,
    config$outcome,
    config$exposures,
    config$intervention_exposures,
    config$exposure_lags,
    config$baseline_covariates,
    config$tvc,
    config$tvc_lags
  ))

  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0L) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  dup_dt <- dt[, .N, by = c(config$id, config$time_in)][N > 1L]
  if (nrow(dup_dt) > 0L) {
    stop("ID-time combinations are not unique.")
  }

  bad_followup <- dt[get(config$time_out) <= get(config$time_in)]
  if (nrow(bad_followup) > 0L) {
    stop("Found rows with non-positive follow-up intervals.")
  }

  numeric_like <- c("normal", "bounded normal", "bounded_normal", "lognormal")
  binary_like <- c("binary", "binomial")
  categorical_like <- c("categorical", "multinomial")

  for (v in config$exposures) {
    if (!is.numeric(dt[[v]])) {
      stop("Exposure must be numeric: ", v)
    }
  }

  for (v in config$tvc) {
    vtype <- tolower(config$covtypes[[v]])
    if (vtype %in% numeric_like && !is.numeric(dt[[v]])) {
      stop("Numeric TVC must be numeric: ", v)
    }
    if (vtype %in% binary_like) {
      vals <- unique(stats::na.omit(dt[[v]]))
      if (!all(vals %in% c(0, 1))) {
        stop("Binary TVC must contain only 0/1 values: ", v)
      }
    }
    if (vtype %in% categorical_like && !is.factor(dt[[v]])) {
      stop("Categorical TVC must be prepared as a factor before modeling: ", v)
    }
  }

  outcome_vals <- unique(stats::na.omit(dt[[config$outcome]]))
  if (!all(outcome_vals %in% c(0, 1))) {
    stop("Outcome must be binary 0/1: ", config$outcome)
  }

  invisible(TRUE)
}
