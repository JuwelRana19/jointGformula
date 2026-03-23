make_lags <- function(data, config) {
  dt <- as.data.table(copy(data))
  data.table::setorderv(dt, c(config$id, config$time_in))

  vars_to_lag <- unique(c(config$exposures, config$tvc))
  lag_names <- c(config$exposure_lags, config$tvc_lags)
  names(lag_names) <- vars_to_lag

  for (v in vars_to_lag) {
    lag_v <- lag_names[[v]]
    dt[, (lag_v) := data.table::shift(get(v), type = "lag"), by = c(config$id)]
    dt[, (lag_v) := data.table::fcase(
      is.na(get(lag_v)), get(v)[1L],
      default = get(lag_v)
    ), by = c(config$id)]
  }

  dt[]
}


