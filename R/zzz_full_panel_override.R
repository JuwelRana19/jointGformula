complete_mc_data <- function(data, config, mc_size = 10000L, seed = 1234L, replace = NULL) {
  dt <- as.data.table(copy(data))
  validate_gformula_data(dt, config)
  data.table::setorderv(dt, c(config$id, config$time_in))

  set.seed(seed)

  mintime <- min(dt[[config$time_in]], na.rm = TRUE)
  times <- sort(unique(dt[[config$time_in]]))
  n_times <- length(times)

  tvc_levels <- lapply(config$tvc, function(v) if (is.factor(dt[[v]])) levels(dt[[v]]) else NULL)
  names(tvc_levels) <- config$tvc

  base_cols <- unique(c(config$id, config$baseline_covariates))
  base_dt <- dt[, ..base_cols][, .SD[1L], by = c(config$id)]
  base_ids <- base_dt[[config$id]]

  if (is.null(replace)) {
    replace <- mc_size > length(base_ids)
  }

  t0_tvc <- dt[get(config$time_in) == mintime, c(config$id, config$tvc), with = FALSE]
  setnames(t0_tvc, config$tvc, paste0(config$tvc, "_t0"))

  t0_exp <- dt[get(config$time_in) == mintime, c(config$id, config$exposures), with = FALSE]
  setnames(t0_exp, config$exposures, paste0(config$exposures, "_t0"))

  setkeyv(base_dt, config$id)
  setkeyv(t0_tvc, config$id)
  setkeyv(t0_exp, config$id)
  base_dt <- t0_tvc[base_dt]
  base_dt <- t0_exp[base_dt]

  draw_ids <- sample(base_ids, size = mc_size, replace = replace)
  sim_base <- base_dt[J(draw_ids), nomatch = 0L, allow.cartesian = TRUE]
  sim_base[, sim_person := .I]
  sim_base[, uid := paste0("sim_", sim_person)]
  sim_base[, original_id := get(config$id)]

  sim_dt <- sim_base[rep(seq_len(.N), each = n_times)]
  sim_dt[, (config$time_in) := rep(times, times = mc_size)]
  sim_dt[, (config$time_out) := get(config$time_in) + 1L]
  sim_dt[, (config$id) := uid]

  intervention_cols <- intersect(config$intervention_exposures, names(dt))
  if (length(intervention_cols) > 0L) {
    int_dt <- dt[, c(config$id, config$time_in, intervention_cols), with = FALSE]
    setnames(int_dt, config$id, "original_id")
    sim_dt <- merge(sim_dt, int_dt, by = c("original_id", config$time_in), all.x = TRUE, sort = FALSE)
  } else {
    for (nm in config$intervention_exposures) {
      sim_dt[, (nm) := NA_real_]
    }
  }

  for (v in config$tvc) {
    vtype <- tolower(config$covtypes[[v]])
    if (vtype %in% c("categorical", "multinomial")) {
      sim_dt[, (v) := factor(NA, levels = tvc_levels[[v]])]
    } else if (vtype %in% c("binary", "binomial")) {
      sim_dt[, (v) := NA_integer_]
    } else {
      sim_dt[, (v) := NA_real_]
    }
    sim_dt[get(config$time_in) == mintime, (v) := get(paste0(v, "_t0"))]
  }

  for (v in config$exposures) {
    sim_dt[, (v) := NA_real_]
    sim_dt[get(config$time_in) == mintime, (v) := get(paste0(v, "_t0"))]
  }

  for (i in seq_along(config$tvc)) {
    v <- config$tvc[i]
    lag_name <- config$tvc_lags[i]
    vtype <- tolower(config$covtypes[[v]])
    if (vtype %in% c("binary", "binomial")) {
      sim_dt[, (lag_name) := NA_integer_]
    } else {
      sim_dt[, (lag_name) := NA_real_]
    }
  }

  for (lag_name in config$exposure_lags) {
    sim_dt[, (lag_name) := NA_real_]
  }

  for (nm in c(paste0(config$tvc, "_t0"), paste0(config$exposures, "_t0"))) {
    if (nm %in% names(sim_dt)) sim_dt[, (nm) := NULL]
  }

  sim_dt[, (config$outcome) := NA_real_]
  sim_dt[, hazard := NA_real_]
  sim_dt[, cumrisk := 0]
  sim_dt[, survival_prob := 1]

  data.table::setorderv(sim_dt, c("uid", config$time_in))
  standardize_output_columns(sim_dt, config, prediction_cols = c("hazard", "cumrisk", "survival_prob"))
}

simulate_joint_gformula <- function(mc_data, config, model_fit, temp = 1.25) {

infer_exposure_order <- function(config) {
  exposures <- config$exposures
  formulas <- config$exposure_formulas[exposures]
  deps <- lapply(exposures, function(v) {
    vars <- all.vars(formulas[[v]])
    dep <- intersect(exposures, vars)
    setdiff(dep, v)
  })
  names(deps) <- exposures

  resolved <- character(0)
  remaining <- exposures
  while (length(remaining) > 0L) {
    ready <- remaining[vapply(remaining, function(v) all(deps[[v]] %in% resolved), logical(1))]
    if (length(ready) == 0L) {
      return(config$exposure_order)
    }
    resolved <- c(resolved, ready)
    remaining <- setdiff(remaining, ready)
  }
  resolved
}

  dt <- as.data.table(copy(mc_data))
  data.table::setorderv(dt, c("uid", config$time_in))

  mintime <- min(dt[[config$time_in]], na.rm = TRUE)
  times <- sort(unique(dt[[config$time_in]]))

  set_factor_from_index <- function(x, idx) {
    factor(levels(x)[idx], levels = levels(x), ordered = is.ordered(x))
  }

  prev_tvc_cols <- c("uid", config$tvc)
  prev_tvc <- dt[get(config$time_in) == mintime, ..prev_tvc_cols]
  setnames(prev_tvc, config$tvc, paste0(config$tvc, "_prev"))
  setkeyv(prev_tvc, "uid")

  prev_exp_cols <- c("uid", config$exposures)
  prev_exp <- dt[get(config$time_in) == mintime, ..prev_exp_cols]
  setnames(prev_exp, config$exposures, paste0(config$exposures, "_prev"))
  setkeyv(prev_exp, "uid")

  for (tt in times) {
    if (tt > mintime) {
      tvc_lag_map <- stats::setNames(paste0(config$tvc, "_prev"), config$tvc_lags)
      for (i in seq_along(config$tvc)) {
        v <- config$tvc[i]
        lag_name <- config$tvc_lags[i]
        prev_name <- tvc_lag_map[[lag_name]]
        vtype <- tolower(config$covtypes[[v]])
        vals <- prev_tvc[dt[get(config$time_in) == tt, .(uid)], on = "uid", get(prev_name)]
        if (vtype %in% c("categorical", "multinomial")) {
          dt[get(config$time_in) == tt, (lag_name) := as.numeric(as.character(vals))]
        } else if (vtype %in% c("binary", "binomial")) {
          dt[get(config$time_in) == tt, (lag_name) := as.integer(vals)]
        } else {
          dt[get(config$time_in) == tt, (lag_name) := as.numeric(vals)]
        }
      }

      exp_lag_map <- stats::setNames(paste0(config$exposures, "_prev"), config$exposure_lags)
      for (lag_name in names(exp_lag_map)) {
        prev_name <- exp_lag_map[[lag_name]]
        dt[get(config$time_in) == tt, (lag_name) := prev_exp[.SD, on = "uid", get(prev_name)]]
      }

      idx <- which(dt[[config$time_in]] == tt)
      id_example <- dt$uid[idx][1]

      for (v in config$tvc) {
        vtype <- tolower(config$covtypes[[v]])
        if (vtype %in% c("categorical", "multinomial")) {
          p <- safe_predict_probs(model_fit$tvc_models[[v]], dt[idx], v, tt, id_example)
          p <- as.matrix(p)
          if (!is.null(colnames(p))) {
            miss <- setdiff(levels(dt[[v]]), colnames(p))
            for (m in miss) {
              p <- cbind(p, stats::setNames(rep(0, nrow(p)), m))
            }
            p <- p[, levels(dt[[v]]), drop = FALSE]
          }
          k <- draw_from_probs_fast(p, temp = temp)
          dt[idx, (v) := set_factor_from_index(get(v), k)]
        } else if (vtype %in% c("binary", "binomial")) {
          dt[idx, (v) := simulate_binary(model_fit$tvc_models[[v]], dt[idx])]
        } else if (vtype %in% c("normal", "bounded normal", "bounded_normal")) {
          dt[idx, (v) := simulate_continuous(model_fit$tvc_models[[v]], dt[idx], model_fit$tvc_support[[v]], vtype)]
        } else {
          stop("Unsupported covtype for simulation: ", config$covtypes[[v]])
        }
      }

      exposure_order <- infer_exposure_order(config)
      for (exposure in exposure_order) {
        vtype <- tolower(config$exposure_types[[exposure]])
        if (vtype %in% c("binary", "binomial")) {
          dt[idx, (exposure) := simulate_binary(model_fit$exposure_models[[exposure]], dt[idx])]
        } else if (vtype %in% c("normal", "bounded normal", "bounded_normal", "lognormal")) {
          dt[idx, (exposure) := simulate_continuous(model_fit$exposure_models[[exposure]], dt[idx], model_fit$exposure_support[[exposure]], vtype)]
        } else {
          stop("Unsupported exposure type for simulation: ", config$exposure_types[[exposure]])
        }
      }
    }

    prev_tvc <- dt[get(config$time_in) == tt, ..prev_tvc_cols]
    setnames(prev_tvc, config$tvc, paste0(config$tvc, "_prev"))
    setkeyv(prev_tvc, "uid")

    prev_exp <- dt[get(config$time_in) == tt, ..prev_exp_cols]
    setnames(prev_exp, config$exposures, paste0(config$exposures, "_prev"))
    setkeyv(prev_exp, "uid")
  }

  standardize_output_columns(dt, config, prediction_cols = c("hazard", "cumrisk", "survival_prob"))
}

predict_joint_outcome <- function(sim_data, config, outcome_model, use_intervention = FALSE) {
  dt <- as.data.table(copy(sim_data))
  group_id <- if ("uid" %in% names(dt)) "uid" else config$id

  if (isTRUE(use_intervention)) {
    for (i in seq_along(config$exposures)) {
      exposure_name <- config$exposures[i]
      intervention_name <- config$intervention_exposures[i]
      if (intervention_name %in% names(dt)) {
        idx <- !is.na(dt[[intervention_name]])
        if (any(idx)) {
          dt[idx, (exposure_name) := get(intervention_name)]
        }
      }
    }
  }

  dt[, hazard := stats::predict(outcome_model, newdata = dt, type = "response")]
  dt[, hazard := pmin(pmax(hazard, 0), 1)]
  data.table::setorderv(dt, c(group_id, config$time_in))
  dt[, survival_prob := cumprod(1 - hazard), by = c(group_id)]
  dt[, cumrisk := 1 - survival_prob]

  standardize_output_columns(dt, config, prediction_cols = c("hazard", "survival_prob", "cumrisk"))
}

run_joint_gformula <- function(data,
                               config,
                               mc_size = 10000L,
                               seed = 1234L,
                               replace_mc = NULL,
                               temp = 1.25,
                               natural_course = config$natural_course %||% "modeled") {
  dt <- as.data.table(copy(data))
  validate_gformula_data(dt, config)

  natural_course <- unique(match.arg(natural_course, c("observed", "modeled"), several.ok = TRUE))

  model_fit <- fit_joint_models(dt, config)
  mc_data <- complete_mc_data(dt, config, mc_size = mc_size, seed = seed, replace = replace_mc)
  sim_data <- simulate_joint_gformula(mc_data, config, model_fit, temp = temp)

  natural_observed <- NULL
  natural_modeled <- NULL
  if ("observed" %in% natural_course) {
    natural_observed <- predict_joint_outcome(dt, config, model_fit$outcome_model, use_intervention = FALSE)
  }
  if ("modeled" %in% natural_course) {
    natural_modeled <- predict_joint_outcome(sim_data, config, model_fit$outcome_model, use_intervention = FALSE)
  }
  natural_default <- if ("observed" %in% natural_course) natural_observed else natural_modeled
  intervention <- predict_joint_outcome(sim_data, config, model_fit$outcome_model, use_intervention = TRUE)

  list(
    estimates = estimate_joint_effects(natural_default, intervention, config),
    observed_data = dt,
    mc_data = mc_data,
    sim_data = sim_data,
    natural = natural_default,
    natural_observed = natural_observed,
    natural_modeled = natural_modeled,
    intervention = intervention,
    model_fit = model_fit,
    config = config
  )
}

create_mc_skeleton <- function(data, config, mc_size = 10000L, seed = 1234L, replace_mc = NULL) {
  complete_mc_data(
    data = data,
    config = config,
    mc_size = mc_size,
    seed = seed,
    replace = replace_mc
  )
}
