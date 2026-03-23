draw_from_probs_fast <- function(pmat, temp = 1, fallback = "uniform") {
  pmat <- as.matrix(pmat)
  if (is.null(dim(pmat))) {
    pmat <- cbind(pmat)
  }

  pmat[!is.finite(pmat)] <- 0
  pmat <- pmax(pmat, 0)

  if (!isTRUE(all.equal(temp, 1))) {
    pmat <- pmat^(1 / temp)
  }

  rs <- rowSums(pmat)
  bad <- rs <= 0

  if (any(bad)) {
    if (identical(fallback, "uniform")) {
      pmat[bad, ] <- 1
      rs[bad] <- ncol(pmat)
    } else {
      stop("Bad probability rows in draw_from_probs_fast().")
    }
  }

  pmat <- pmat / rs
  cs <- pmat
  if (ncol(pmat) > 1L) {
    for (j in 2:ncol(pmat)) {
      cs[, j] <- cs[, j] + cs[, j - 1L]
    }
  }
  cs[, ncol(pmat)] <- 1

  u <- runif(nrow(cs))
  max.col(u <= cs, ties.method = "first")
}


safe_predict_probs <- function(model, newdata, varname, time_value, id_value = NULL) {
  p <- tryCatch(
    predict(model, newdata = newdata, type = "probs"),
    error = function(e) {
      stop(
        "Predict error for ", varname,
        " at time ", time_value,
        if (!is.null(id_value)) paste0(" (id=", id_value, ")") else "",
        ": ", conditionMessage(e)
      )
    }
  )

  p <- as.matrix(p)
  bad_rows <- !is.finite(rowSums(p, na.rm = TRUE)) | rowSums(p, na.rm = TRUE) <= 0
  if (any(bad_rows)) {
    stop("Invalid probability rows for ", varname, " at time ", time_value)
  }

  p
}


fit_gaussian_support <- function(model, observed, value_type, range = NULL) {
  out <- list(
    type = value_type,
    sigma = stats::sd(stats::residuals(model, type = "response"), na.rm = TRUE),
    min = NA_real_,
    max = NA_real_
  )

  if (tolower(value_type) %in% c("bounded normal", "bounded_normal", "normal", "lognormal")) {
    if (!is.null(range)) {
      out$min <- range[[1]]
      out$max <- range[[2]]
    } else {
      out$min <- min(observed, na.rm = TRUE)
      out$max <- max(observed, na.rm = TRUE)
    }
  }

  out
}


simulate_binary <- function(model, newdata) {
  p <- as.numeric(stats::predict(model, newdata = newdata, type = "response"))
  stats::rbinom(length(p), size = 1L, prob = pmin(pmax(p, 0), 1))
}


simulate_continuous <- function(model, newdata, support, value_type) {
  pred <- as.numeric(stats::predict(model, newdata = newdata, type = "response"))
  value_type <- tolower(value_type)

  if (value_type == "lognormal") {
    log_sim <- pred - 0.5 * support$sigma^2 + stats::rnorm(length(pred), sd = support$sigma)
    vals <- exp(log_sim)
  } else {
    vals <- pred + stats::rnorm(length(pred), sd = support$sigma)
  }

  if (value_type %in% c("bounded normal", "bounded_normal", "lognormal")) {
    vals <- pmax(pmin(vals, support$max), support$min)
  }

  vals
}




force_levels_dt <- function(dt, reference_dt, vars) {
  for (v in vars) {
    if (v %in% names(dt) && v %in% names(reference_dt) && is.factor(reference_dt[[v]])) {
      dt[, (v) := factor(
        as.character(get(v)),
        levels = levels(reference_dt[[v]]),
        ordered = is.ordered(reference_dt[[v]])
      )]
    }
  }
  dt[]
}

standardize_output_columns <- function(dt, config, prediction_cols = NULL) {
  id_cols <- intersect(c("uid", "sim_person", "original_id", config$id, config$time_in, config$time_out), names(dt))
  base_cols <- intersect(config$baseline_covariates, names(dt))
  tvc_cols <- intersect(config$tvc, names(dt))
  tvc_lag_cols <- intersect(config$tvc_lags, names(dt))
  exposure_lag_cols <- intersect(config$exposure_lags, names(dt))
  exposure_cols <- intersect(config$exposures, names(dt))
  outcome_cols <- intersect(config$outcome, names(dt))
  prediction_cols <- intersect(prediction_cols, names(dt))

  ordered_cols <- unique(c(
    id_cols,
    base_cols,
    tvc_cols,
    tvc_lag_cols,
    exposure_lag_cols,
    exposure_cols,
    outcome_cols,
    prediction_cols
  ))

  data.table::setcolorder(dt, c(ordered_cols, setdiff(names(dt), ordered_cols)))
  dt[]
}
complete_mc_data <- function(data, config, mc_size = 10000L, seed = 1234L, replace = NULL) {
  dt <- as.data.table(copy(data))
  validate_gformula_data(dt, config)
  data.table::setorderv(dt, c(config, config))

  set.seed(seed)

  mintime <- min(dt[[config]], na.rm = TRUE)
  times <- sort(unique(dt[[config]]))
  n_times <- length(times)

  tvc_levels <- lapply(config, function(v) if (is.factor(dt[[v]])) levels(dt[[v]]) else NULL)
  names(tvc_levels) <- config

  base_cols <- unique(c(config, config))
  base_dt <- dt[, ..base_cols][, .SD[1L], by = c(config)]
  base_ids <- base_dt[[config]]

  if (is.null(replace)) {
    replace <- mc_size > length(base_ids)
  }

  t0_tvc <- dt[get(config) == mintime, c(config, config), with = FALSE]
  data.table::setnames(t0_tvc, config, paste0(config, "_t0"))

  t0_exp <- dt[get(config) == mintime, c(config, config), with = FALSE]
  data.table::setnames(t0_exp, config, paste0(config, "_t0"))

  data.table::setkeyv(base_dt, config)
  data.table::setkeyv(t0_tvc, config)
  data.table::setkeyv(t0_exp, config)
  base_dt <- t0_tvc[base_dt]
  base_dt <- t0_exp[base_dt]

  draw_ids <- sample(base_ids, size = mc_size, replace = replace)
  sim_base <- base_dt[J(draw_ids), nomatch = 0L, allow.cartesian = TRUE]
  sim_base[, sim_person := .I]
  sim_base[, uid := paste0("sim_", sim_person)]
  sim_base[, original_id := get(config)]

  sim_dt <- sim_base[rep(seq_len(.N), each = n_times)]
  sim_dt[, (config) := rep(times, times = mc_size)]
  sim_dt[, (config) := get(config) + 1L]
  sim_dt[, (config) := uid]

  intervention_cols <- intersect(config, names(dt))
  if (length(intervention_cols) > 0L) {
    int_dt <- dt[, c(config, config, intervention_cols), with = FALSE]
    data.table::setnames(int_dt, config, "original_id")
    sim_dt <- merge(sim_dt, int_dt, by = c("original_id", config), all.x = TRUE, sort = FALSE)
  } else {
    for (nm in config) {
      sim_dt[, (nm) := NA_real_]
    }
  }

  for (v in config) {
    vtype <- tolower(config[[v]])
    if (vtype %in% c("categorical", "multinomial")) {
      sim_dt[, (v) := factor(NA, levels = tvc_levels[[v]])]
    } else if (vtype %in% c("binary", "binomial")) {
      sim_dt[, (v) := NA_integer_]
    } else {
      sim_dt[, (v) := NA_real_]
    }
    sim_dt[get(config) == mintime, (v) := get(paste0(v, "_t0"))]
  }

  for (v in config) {
    sim_dt[, (v) := NA_real_]
    sim_dt[get(config) == mintime, (v) := get(paste0(v, "_t0"))]
  }

  for (i in seq_along(config)) {
    v <- config[i]
    lag_name <- config[i]
    vtype <- tolower(config[[v]])
    if (vtype %in% c("categorical", "multinomial")) {
      sim_dt[, (lag_name) := factor(NA, levels = tvc_levels[[v]])]
    } else if (vtype %in% c("binary", "binomial")) {
      sim_dt[, (lag_name) := NA_integer_]
    } else {
      sim_dt[, (lag_name) := NA_real_]
    }
  }

  for (lag_name in config) {
    sim_dt[, (lag_name) := NA_real_]
  }

  t0_names <- c(paste0(config, "_t0"), paste0(config, "_t0"))
  for (nm in t0_names) {
    if (nm %in% names(sim_dt)) sim_dt[, (nm) := NULL]
  }

  sim_dt[, (config) := NA_real_]
  sim_dt[, hazard := NA_real_]
  sim_dt[, cumrisk := 0]
  sim_dt[, survival_prob := 1]

  data.table::setorderv(sim_dt, c("uid", config))
  standardize_output_columns(sim_dt, config, prediction_cols = c("hazard", "cumrisk", "survival_prob"))
}
fit_joint_models <- function(data, config) {
  dt <- as.data.table(copy(data))
  validate_gformula_data(dt, config)

  outcome_model <- speedglm(
    formula = config$outcome_formula,
    data = dt,
    sparse = FALSE,
    family = stats::binomial("logit")
  )

  tvc_models <- lapply(names(config$tvc_formulas), function(v) {
    fm <- config$tvc_formulas[[v]]
    vtype <- tolower(config$covtypes[[v]])

    if (vtype %in% c("categorical", "multinomial")) {
      multinom(formula = fm, data = dt, trace = FALSE)
    } else if (vtype %in% c("binary", "binomial")) {
      stats::glm(formula = fm, data = dt, family = stats::binomial())
    } else if (vtype %in% c("normal", "bounded normal", "bounded_normal")) {
      stats::glm(formula = fm, data = dt, family = stats::gaussian())
    } else {
      stop("Unsupported covtype for ", v, ": ", config$covtypes[[v]])
    }
  })
  names(tvc_models) <- names(config$tvc_formulas)

  exposure_models <- lapply(names(config$exposure_formulas), function(v) {
    fm <- config$exposure_formulas[[v]]
    vtype <- tolower(config$exposure_types[[v]])

    if (vtype %in% c("normal", "bounded normal", "bounded_normal", "lognormal")) {
      stats::glm(formula = fm, data = dt, family = stats::gaussian())
    } else if (vtype %in% c("binary", "binomial")) {
      stats::glm(formula = fm, data = dt, family = stats::binomial())
    } else {
      stop("Unsupported exposure type for ", v, ": ", config$exposure_types[[v]])
    }
  })
  names(exposure_models) <- names(config$exposure_formulas)

  tvc_support <- lapply(config$tvc, function(v) {
    vtype <- config$covtypes[[v]]
    model <- tvc_models[[v]]
    if (tolower(vtype) %in% c("normal", "bounded normal", "bounded_normal")) {
      range <- NULL
      if (!is.null(config$covranges) && !is.null(config$covranges[[v]])) {
        range <- config$covranges[[v]]
      }
      fit_gaussian_support(model, dt[[v]], vtype, range = range)
    } else {
      NULL
    }
  })
  names(tvc_support) <- config$tvc

  exposure_support <- lapply(config$exposures, function(v) {
    vtype <- config$exposure_types[[v]]
    model <- exposure_models[[v]]
    if (tolower(vtype) %in% c("normal", "bounded normal", "bounded_normal", "lognormal")) {
      range <- NULL
      if (!is.null(config$exposure_ranges) && !is.null(config$exposure_ranges[[v]])) {
        range <- config$exposure_ranges[[v]]
      }
      fit_gaussian_support(model, dt[[v]], vtype, range = range)
    } else {
      NULL
    }
  })
  names(exposure_support) <- config$exposures

  list(
    outcome_model = outcome_model,
    tvc_models = tvc_models,
    exposure_models = exposure_models,
    tvc_support = tvc_support,
    exposure_support = exposure_support
  )
}


simulate_joint_gformula <- function(mc_data, config, model_fit, temp = 1.25) {
  dt <- as.data.table(copy(mc_data))
  data.table::setorderv(dt, c("uid", config))

  mintime <- min(dt[[config]], na.rm = TRUE)
  times <- sort(unique(dt[[config]]))

  set_factor_from_index <- function(x, idx) {
    factor(levels(x)[idx], levels = levels(x), ordered = is.ordered(x))
  }

  prev_tvc_cols <- c("uid", config)
  prev_tvc <- dt[get(config) == mintime, ..prev_tvc_cols]
  data.table::setnames(prev_tvc, config, paste0(config, "_prev"))
  data.table::setkeyv(prev_tvc, "uid")

  prev_exp_cols <- c("uid", config)
  prev_exp <- dt[get(config) == mintime, ..prev_exp_cols]
  data.table::setnames(prev_exp, config, paste0(config, "_prev"))
  data.table::setkeyv(prev_exp, "uid")

  for (tt in times) {
    if (tt > mintime) {
      tvc_lag_map <- stats::setNames(paste0(config, "_prev"), config)
      for (lag_name in names(tvc_lag_map)) {
        prev_name <- tvc_lag_map[[lag_name]]
        dt[get(config) == tt, (lag_name) := prev_tvc[.SD, on = "uid", get(prev_name)]]
      }

      exp_lag_map <- stats::setNames(paste0(config, "_prev"), config)
      for (lag_name in names(exp_lag_map)) {
        prev_name <- exp_lag_map[[lag_name]]
        dt[get(config) == tt, (lag_name) := prev_exp[.SD, on = "uid", get(prev_name)]]
      }

      idx <- which(dt[[config]] == tt)
      id_example <- dt[idx][1]

      for (v in config) {
        vtype <- tolower(config[[v]])
        if (vtype %in% c("categorical", "multinomial")) {
          p <- safe_predict_probs(model_fit[[v]], dt[idx], v, tt, id_example)
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
          dt[idx, (v) := simulate_binary(model_fit[[v]], dt[idx])]
        } else if (vtype %in% c("normal", "bounded normal", "bounded_normal")) {
          dt[idx, (v) := simulate_continuous(model_fit[[v]], dt[idx], model_fit[[v]], vtype)]
        } else {
          stop("Unsupported covtype for simulation: ", config[[v]])
        }
      }

      for (exposure in config) {
        vtype <- tolower(config[[exposure]])
        if (vtype %in% c("binary", "binomial")) {
          dt[idx, (exposure) := simulate_binary(model_fit[[exposure]], dt[idx])]
        } else if (vtype %in% c("normal", "bounded normal", "bounded_normal", "lognormal")) {
          dt[idx, (exposure) := simulate_continuous(model_fit[[exposure]], dt[idx], model_fit[[exposure]], vtype)]
        } else {
          stop("Unsupported exposure type for simulation: ", config[[exposure]])
        }
      }
    }

    prev_tvc <- dt[get(config) == tt, ..prev_tvc_cols]
    data.table::setnames(prev_tvc, config, paste0(config, "_prev"))
    data.table::setkeyv(prev_tvc, "uid")

    prev_exp <- dt[get(config) == tt, ..prev_exp_cols]
    data.table::setnames(prev_exp, config, paste0(config, "_prev"))
    data.table::setkeyv(prev_exp, "uid")
  }

  standardize_output_columns(dt, config, prediction_cols = c("hazard", "cumrisk", "survival_prob"))
}
predict_joint_outcome <- function(sim_data, config, outcome_model, use_intervention = FALSE) {
  dt <- as.data.table(copy(sim_data))
  group_id <- if ("uid" %in% names(dt)) "uid" else config

  if (isTRUE(use_intervention)) {
    for (i in seq_along(config)) {
      exposure_name <- config[i]
      intervention_name <- config[i]
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
  data.table::setorderv(dt, c(group_id, config))
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




format_progress_time <- function(seconds) {
  if (!is.finite(seconds) || is.na(seconds)) return('NA')
  seconds <- max(0, round(seconds))
  hrs <- seconds %/% 3600
  mins <- (seconds %% 3600) %/% 60
  secs <- seconds %% 60
  if (hrs > 0) {
    sprintf('%02d:%02d:%02d', hrs, mins, secs)
  } else {
    sprintf('%02d:%02d', mins, secs)
  }
}
print_progress_bar <- function(done, total, start_time, label = 'Bootstraps') {
  total <- as.integer(total)
  done <- as.integer(done)
  width <- 24L
  frac <- if (total > 0L) max(0, min(1, done / total)) else 1
  filled <- as.integer(floor(width * frac))
  empty <- width - filled
  elapsed_sec <- as.numeric(difftime(Sys.time(), start_time, units = 'secs'))
  eta_sec <- if (done > 0L) elapsed_sec / done * (total - done) else NA_real_
  green_on <- '\033[32m'
  color_off <- '\033[39m'
  bar <- paste0(green_on, strrep('=', filled), color_off, strrep('-', empty))
  pct <- sprintf('%3d%%', round(frac * 100))
  cat(sprintf('%s [%s] %s %d/%d | elapsed %s | ETA %s\n',
              label,
              bar,
              pct,
              done,
              total,
              format_progress_time(elapsed_sec),
              format_progress_time(eta_sec)))
}

bootstrap_joint_gformula <- function(data,
                                     config,
                                     n_boot = 200L,
                                     mc_size = 10000L,
                                     seed = 1234L,
                                     replace_mc = NULL,
                                     temp = 1.25,
                                     natural_course = config$natural_course %||% "modeled",
                                     parallel = FALSE,
                                     n_workers = 2L,
                                     batch_size = NULL,
                                     checkpoint_file = NULL,
                                     verbose = TRUE,
                                     audit_mode = FALSE,
                                     keep_boot_objects = FALSE) {
  dt <- as.data.table(copy(data))
  validate_gformula_data(dt, config)
  data.table::setorderv(dt, c(config$id, config$time_in))

  ids <- unique(dt[[config$id]])
  level_vars <- intersect(config$factor_vars, names(dt))
  start_time <- Sys.time()
  keep_boot_objects <- isTRUE(keep_boot_objects) || isTRUE(audit_mode)

  runner <- function(b) {
    set.seed(seed + b)
    boot_ids <- sample(ids, size = length(ids), replace = TRUE)
    boot_map <- data.table(tmp_id = boot_ids, boot_id = seq_along(boot_ids))
    setnames(boot_map, "tmp_id", config$id)
    boot_dt <- merge(boot_map, dt, by = config$id, allow.cartesian = TRUE)
    boot_dt[, original_id := get(config$id)]
    boot_dt[, (config$id) := paste0(boot_id, "_", original_id)]
    boot_dt <- force_levels_dt(boot_dt, dt, level_vars)

    out <- tryCatch(
      run_joint_gformula(
        data = boot_dt,
        config = config,
        mc_size = mc_size,
        seed = seed + 100000L + b,
        replace_mc = replace_mc,
        temp = temp,
        natural_course = natural_course
      ),
      error = function(e) e
    )

    if (inherits(out, "error")) {
      return(list(
        result = data.table(
          boot = b,
          risk_natural = NA_real_,
          risk_intervention = NA_real_,
          RD = NA_real_,
          RR = NA_real_,
          status = paste("error:", conditionMessage(out))
        ),
        object = if (keep_boot_objects) list(boot_data = boot_dt, error = conditionMessage(out)) else NULL
      ))
    }

    est <- out$estimates
    result_dt <- data.table(
      boot = b,
      risk_natural = est$risk_natural,
      risk_intervention = est$risk_intervention,
      RD = est$RD,
      RR = est$RR,
      status = "success"
    )

    audit_object <- NULL
    if (keep_boot_objects) {
      audit_object <- list(
        boot_data = boot_dt,
        mc_data = out$mc_data,
        sim_data = out$sim_data,
        natural = out$natural,
        intervention = out$intervention,
        estimates = out$estimates
      )
    }

    list(result = result_dt[], object = audit_object)
  }

  results <- vector("list", n_boot)
  audit_objects <- if (keep_boot_objects) vector("list", n_boot) else NULL

  save_checkpoint <- function(last_completed) {
    if (!is.null(checkpoint_file)) {
      payload <- list(results = results, last_completed = last_completed)
      if (keep_boot_objects) {
        payload$audit_objects <- audit_objects
      }
      saveRDS(payload, checkpoint_file)
    }
  }

  if (isTRUE(parallel)) {
    if (!requireNamespace("future", quietly = TRUE) || !requireNamespace("future.apply", quietly = TRUE)) {
      stop("Parallel bootstrap requires the 'future' and 'future.apply' packages.")
    }
    if (is.null(batch_size)) {
      batch_size <- max(1L, as.integer(n_workers))
    }
    future::plan(future::multisession, workers = n_workers)
    on.exit(future::plan(future::sequential), add = TRUE)

    batches <- split(seq_len(n_boot), ceiling(seq_len(n_boot) / batch_size))
    for (i in seq_along(batches)) {
      batch_ids <- batches[[i]]
      if (verbose) {
        cat("Parallel batch", i, "of", length(batches), "\n")
      }
      batch_res <- future.apply::future_lapply(batch_ids, runner, future.seed = TRUE)
      for (j in seq_along(batch_ids)) {
        results[[batch_ids[j]]] <- batch_res[[j]]$result
        if (keep_boot_objects) {
          audit_objects[[batch_ids[j]]] <- batch_res[[j]]$object
        }
      }
      save_checkpoint(max(batch_ids))
      if (verbose) {
        done <- max(batch_ids)
        print_progress_bar(done, n_boot, start_time, label = "Bootstraps")
      }
    }
  } else {
    for (b in seq_len(n_boot)) {
      if (verbose) {
        cat("Bootstrap", b, "of", n_boot, "\n")
      }
      run_res <- runner(b)
      results[[b]] <- run_res$result
      if (keep_boot_objects) {
        audit_objects[[b]] <- run_res$object
      }
      save_checkpoint(b)
      if (verbose) {
        print_progress_bar(b, n_boot, start_time, label = "Bootstraps")
      }
    }
  }

  out <- rbindlist(results, fill = TRUE)
  if (!is.null(checkpoint_file)) {
    if (keep_boot_objects) {
      saveRDS(list(results = out, audit_objects = audit_objects), checkpoint_file)
    } else {
      saveRDS(out, checkpoint_file)
    }
  }

  if (keep_boot_objects) {
    return(list(results = out[], audit_objects = audit_objects, audit_mode = isTRUE(audit_mode)))
  }
  out[]
}


summarize_joint_bootstrap <- function(boot_res) {
  br <- if (is.list(boot_res) && !data.table::is.data.table(boot_res) && !is.null(boot_res$results)) boot_res$results else boot_res
  br <- as.data.table(copy(br))
  br <- br[status == "success"]

  qfun <- function(x) as.numeric(stats::quantile(x, c(0.025, 0.975), na.rm = TRUE, names = FALSE))

  data.table(
    metric = c("risk_natural", "risk_intervention", "RD", "RR"),
    estimate = c(
      mean(br$risk_natural, na.rm = TRUE),
      mean(br$risk_intervention, na.rm = TRUE),
      mean(br$RD, na.rm = TRUE),
      mean(br$RR, na.rm = TRUE)
    ),
    lower = c(
      qfun(br$risk_natural)[1],
      qfun(br$risk_intervention)[1],
      qfun(br$RD)[1],
      qfun(br$RR)[1]
    ),
    upper = c(
      qfun(br$risk_natural)[2],
      qfun(br$risk_intervention)[2],
      qfun(br$RD)[2],
      qfun(br$RR)[2]
    )
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


assign_natural_course <- function(sim_data) {
  as.data.table(copy(sim_data))
}


assign_joint_intervention <- function(sim_data, config, intervention_map = NULL) {
  dt <- as.data.table(copy(sim_data))

  if (is.null(intervention_map)) {
    intervention_map <- stats::setNames(config$intervention_exposures, config$exposures)
  }

  for (nm in names(intervention_map)) {
    dt[, (nm) := get(intervention_map[[nm]])]
  }

  dt[]
}


estimate_joint_effects <- function(natural, intervention, config) {
  final_time <- max(natural[[config$time_in]], na.rm = TRUE)
  risk_natural <- natural[get(config$time_in) == final_time, mean(cumrisk, na.rm = TRUE)]
  risk_intervention <- intervention[get(config$time_in) == final_time, mean(cumrisk, na.rm = TRUE)]

  data.table(
    risk_natural = risk_natural,
    risk_intervention = risk_intervention,
    RD = risk_natural - risk_intervention,
    RR = risk_natural / risk_intervention,
    risk_natural_per_1000 = risk_natural * 1000,
    risk_intervention_per_1000 = risk_intervention * 1000,
    RD_per_1000 = (risk_natural - risk_intervention) * 1000
  )
}


run_joint_pipeline <- function(data,
                               config,
                               mc_size = 10000L,
                               seed = 1234L,
                               replace_mc = NULL,
                               temp = 1.25,
                               bootstrap = FALSE,
                               natural_course = config$natural_course %||% "modeled",
                               n_boot = 200L,
                               parallel = FALSE,
                               n_workers = 2L,
                               batch_size = NULL,
                               checkpoint_file = NULL,
                               verbose = TRUE) {
  if (!isTRUE(bootstrap)) {
    return(run_joint_gformula(
      data = data,
      config = config,
      mc_size = mc_size,
      seed = seed,
      replace_mc = replace_mc,
      temp = temp,
      natural_course = natural_course
    ))
  }

  boot_res <- bootstrap_joint_gformula(
    data = data,
    config = config,
    n_boot = n_boot,
    mc_size = mc_size,
    seed = seed,
    replace_mc = replace_mc,
    temp = temp,
    natural_course = natural_course,
    parallel = parallel,
    n_workers = n_workers,
    batch_size = batch_size,
    checkpoint_file = checkpoint_file,
    verbose = verbose
  )

  list(
    bootstrap_results = boot_res,
    bootstrap_summary = summarize_joint_bootstrap(boot_res),
    config = config
  )
}




describe_model_plan <- function(config) {
  outcome_dt <- data.table(
    component = "outcome",
    variable = config$outcome,
    type = config$outcome_type,
    fit = "speedglm binomial(logit)",
    formula = paste(deparse(config$outcome_formula), collapse = " ")
  )

  tvc_dt <- rbindlist(lapply(names(config$tvc_formulas), function(v) {
    vtype <- tolower(config$covtypes[[v]])
    fit <- if (vtype %in% c("categorical", "multinomial")) {
      "nnet::multinom"
    } else if (vtype %in% c("binary", "binomial")) {
      "stats::glm binomial"
    } else if (vtype %in% c("normal", "bounded normal", "bounded_normal")) {
      "stats::glm gaussian"
    } else {
      "unsupported"
    }

    data.table(
      component = "time_varying_covariate",
      variable = v,
      type = config$covtypes[[v]],
      fit = fit,
      formula = paste(deparse(config$tvc_formulas[[v]]), collapse = " ")
    )
  }), fill = TRUE)

  exposure_dt <- rbindlist(lapply(names(config$exposure_formulas), function(v) {
    vtype <- tolower(config$exposure_types[[v]])
    fit <- if (vtype %in% c("binary", "binomial")) {
      "stats::glm binomial"
    } else if (vtype %in% c("normal", "bounded normal", "bounded_normal", "lognormal")) {
      "stats::glm gaussian"
    } else {
      "unsupported"
    }

    data.table(
      component = "exposure",
      variable = v,
      type = config$exposure_types[[v]],
      fit = fit,
      formula = paste(deparse(config$exposure_formulas[[v]]), collapse = " ")
    )
  }), fill = TRUE)

  rbindlist(list(outcome_dt, tvc_dt, exposure_dt), fill = TRUE)
}

summarize_risk_trajectory <- function(sim_data, config, risk_col = "cumrisk") {
  dt <- as.data.table(copy(sim_data))
  if (!risk_col %in% names(dt)) {
    stop("Risk column not found: ", risk_col)
  }

  dt[, .(
    mean_risk = mean(get(risk_col), na.rm = TRUE),
    mean_risk_per_1000 = 1000 * mean(get(risk_col), na.rm = TRUE)
  ), by = c(config$time_in)]
}



normalize_static_scenario <- function(scenario_values, config, scenario_name = NULL) {
  exposures <- config$exposures

  if (length(scenario_values) == 1L) {
    vals <- rep(as.numeric(scenario_values), length(exposures))
    names(vals) <- exposures
    return(vals)
  }

  if (is.null(names(scenario_values))) {
    if (length(scenario_values) != length(exposures)) {
      stop("Scenario ", scenario_name %||% "", " must have length 1 or match the number of exposures.")
    }
    vals <- as.numeric(scenario_values)
    names(vals) <- exposures
    return(vals)
  }

  if (!all(exposures %in% names(scenario_values))) {
    stop("Scenario ", scenario_name %||% "", " is missing exposure values for: ", paste(setdiff(exposures, names(scenario_values)), collapse = ", "))
  }

  vals <- as.numeric(scenario_values[exposures])
  names(vals) <- exposures
  vals
}

`%||%` <- function(x, y) if (is.null(x)) y else x

prepare_static_intervention_data <- function(mc_data, config, scenario_values, scenario_name = NULL) {
  dt <- as.data.table(copy(mc_data))
  vals <- normalize_static_scenario(scenario_values, config, scenario_name = scenario_name)

  for (i in seq_along(config$exposures)) {
    exp_name <- config$exposures[i]
    lag_name <- config$exposure_lags[i]
    dt[, (exp_name) := vals[[exp_name]]]
    if (lag_name %in% names(dt)) {
      dt[, (lag_name) := vals[[exp_name]]]
    }
  }

  standardize_output_columns(dt, config)
}

simulate_static_intervention <- function(mc_data, config, model_fit, scenario_values, temp = 1.25, scenario_name = NULL) {
  dt <- prepare_static_intervention_data(mc_data, config, scenario_values, scenario_name = scenario_name)
  data.table::setkeyv(dt, c("uid", config$time_in))

  mintime <- min(dt[[config$time_in]], na.rm = TRUE)
  times <- sort(unique(dt[[config$time_in]]))

  set_factor_from_index <- function(x, idx) {
    factor(levels(x)[idx], levels = levels(x), ordered = is.ordered(x))
  }

  prev_cols <- c("uid", config$tvc)
  prev <- dt[get(config$time_in) == mintime, ..prev_cols]
  setnames(prev, config$tvc, paste0(config$tvc, "_prev"))
  data.table::setkeyv(prev, "uid")

  for (tt in times) {
    if (tt > mintime) {
      lag_map <- stats::setNames(paste0(config$tvc, "_prev"), config$tvc_lags)
      for (lag_name in names(lag_map)) {
        prev_name <- lag_map[[lag_name]]
        dt[get(config$time_in) == tt, (lag_name) := prev[.SD, on = "uid", get(prev_name)]]
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
    }

    prev <- dt[get(config$time_in) == tt, ..prev_cols]
    setnames(prev, config$tvc, paste0(config$tvc, "_prev"))
    data.table::setkeyv(prev, "uid")
  }

  standardize_output_columns(dt, config)
}

fit_pooled_hazard_meta_model <- function(simulated_data,
                                         config,
                                         meta_formula = NULL,
                                         meta_family = stats::gaussian()) {
  dt <- as.data.table(copy(simulated_data))

  if (is.null(meta_formula)) {
    meta_formula <- stats::as.formula(
      paste("hazard ~ joint_effect +", config$time_out, "+ I(", config$time_out, "^2)")
    )
  }

  invalid_hazard <- !is.finite(dt$hazard)
  if (any(invalid_hazard)) {
    bad <- dt[invalid_hazard]
    stop(
      paste0(
        "Pooled meta model received ", nrow(bad), " rows with invalid hazard values. ",
        "Check simulated_data for missing predictors before fitting the meta model. ",
        "Affected scenarios: ", paste(unique(bad$scenario), collapse = ", ")
      )
    )
  }

  meta_vars <- all.vars(meta_formula)
  missing_meta <- names(which(vapply(dt[, ..meta_vars], function(x) any(!is.finite(x) | is.na(x)), logical(1))))
  if (length(missing_meta) > 0L) {
    stop(
      paste0(
        "Meta model predictors contain missing or non-finite values: ",
        paste(missing_meta, collapse = ", ")
      )
    )
  }

  stats::glm(formula = meta_formula, data = dt, family = meta_family)
}

summarize_static_scenario_risks <- function(simulated_data, config) {
  dt <- as.data.table(copy(simulated_data))
  dt[, .(
    mean_hazard = mean(hazard, na.rm = TRUE),
    mean_cumrisk = mean(cumrisk, na.rm = TRUE),
    mean_cumrisk_per_1000 = 1000 * mean(cumrisk, na.rm = TRUE)
  ), by = c("scenario", config$time_in)]
}

run_static_intervention_meta <- function(data,
                                         config,
                                         intervention_scenarios,
                                         mc_size = 10000L,
                                         seed = 1234L,
                                         replace_mc = NULL,
                                         temp = 1.25,
                                         joint_effect_fn = NULL,
                                         meta_formula = NULL,
                                         meta_family = stats::gaussian()) {
  dt <- as.data.table(copy(data))
  validate_gformula_data(dt, config)

  model_fit <- fit_joint_models(dt, config)
  mc_data <- complete_mc_data(dt, config, mc_size = mc_size, seed = seed, replace = replace_mc)

  scenario_data <- rbindlist(lapply(names(intervention_scenarios), function(scn) {
    vals <- normalize_static_scenario(intervention_scenarios[[scn]], config, scenario_name = scn)
    sim_dt <- simulate_static_intervention(mc_data, config, model_fit, vals, temp = temp, scenario_name = scn)
    pred_dt <- predict_joint_outcome(sim_dt, config, model_fit$outcome_model, use_intervention = FALSE)

    if (any(!is.finite(pred_dt$hazard))) {
      outcome_vars <- intersect(all.vars(config$outcome_formula), names(pred_dt))
      bad <- pred_dt[!is.finite(hazard)]
      missing_summary <- bad[, vapply(.SD, function(x) sum(is.na(x) | !is.finite(x)), numeric(1)), .SDcols = outcome_vars]
      missing_summary <- missing_summary[missing_summary > 0]
      stop(
        paste0(
          "Static intervention scenario '", scn, "' produced invalid hazard predictions. ",
          "Rows with invalid hazard: ", nrow(bad), ". ",
          if (length(missing_summary) > 0L) paste0("Variables with missing/non-finite values among those rows: ", paste(names(missing_summary), collapse = ", "), ". ") else "No missing predictors were detected in the outcome formula variables. ",
          "This suggests the issue is upstream in scenario construction or model prediction."
        )
      )
    }

    for (exp_name in names(vals)) {
      pred_dt[, (paste0(exp_name, "_assigned")) := vals[[exp_name]]]
    }

    joint_value <- if (is.null(joint_effect_fn)) vals[[config$exposures[1L]]] else joint_effect_fn(vals)
    pred_dt[, joint_effect := joint_value]
    pred_dt[, scenario := scn]
    pred_dt
  }), fill = TRUE)

  meta_model <- fit_pooled_hazard_meta_model(
    simulated_data = scenario_data,
    config = config,
    meta_formula = meta_formula,
    meta_family = meta_family
  )

  list(
    outcome_model = model_fit$outcome_model,
    tvc_models = model_fit$tvc_models,
    exposure_models = model_fit$exposure_models,
    mc_data = mc_data,
    simulated_data = scenario_data,
    scenario_summary = summarize_static_scenario_risks(scenario_data, config),
    meta_model = meta_model,
    config = config
  )
}

bootstrap_static_intervention_meta <- function(data,
                                               config,
                                               intervention_scenarios,
                                               n_boot = 200L,
                                               mc_size = 10000L,
                                               seed = 1234L,
                                               replace_mc = NULL,
                                               temp = 1.25,
                                               joint_effect_fn = NULL,
                                               meta_formula = NULL,
                                               meta_family = stats::gaussian(),
                                               parallel = FALSE,
                                               n_workers = 2L,
                                               batch_size = NULL,
                                               checkpoint_file = NULL,
                                               verbose = TRUE) {
  dt <- as.data.table(copy(data))
  validate_gformula_data(dt, config)
  data.table::setorderv(dt, c(config$id, config$time_in))

  ids <- unique(dt[[config$id]])
  level_vars <- intersect(config$factor_vars, names(dt))
  start_time <- Sys.time()

  runner <- function(b) {
    set.seed(seed + b)
    boot_ids <- sample(ids, size = length(ids), replace = TRUE)
    boot_map <- data.table(tmp_id = boot_ids, boot_id = seq_along(boot_ids))
    setnames(boot_map, "tmp_id", config$id)
    boot_dt <- merge(boot_map, dt, by = config$id, allow.cartesian = TRUE)
    boot_dt[, original_id := get(config$id)]
    boot_dt[, (config$id) := paste0(boot_id, "_", original_id)]
    boot_dt <- force_levels_dt(boot_dt, dt, level_vars)

    out <- tryCatch(
      run_static_intervention_meta(
        data = boot_dt,
        config = config,
        intervention_scenarios = intervention_scenarios,
        mc_size = mc_size,
        seed = seed + 100000L + b,
        replace_mc = replace_mc,
        temp = temp,
        joint_effect_fn = joint_effect_fn,
        meta_formula = meta_formula,
        meta_family = meta_family
      ),
      error = function(e) e
    )

    if (inherits(out, "error")) {
      return(data.table(boot = b, status = paste("error:", conditionMessage(out))))
    }

    coefs <- stats::coef(out$meta_model)
    coef_dt <- as.data.table(as.list(coefs))
    setnames(coef_dt, names(coef_dt), paste0("coef_", names(coef_dt)))
    coef_dt[, boot := b]
    coef_dt[, status := "success"]
    data.table::setcolorder(coef_dt, c("boot", "status", setdiff(names(coef_dt), c("boot", "status"))))
    coef_dt[]
  }

  results <- vector("list", n_boot)

  save_checkpoint <- function(last_completed) {
    if (!is.null(checkpoint_file)) {
      saveRDS(list(results = results, last_completed = last_completed), checkpoint_file)
    }
  }

  if (isTRUE(parallel)) {
    if (!requireNamespace("future", quietly = TRUE) || !requireNamespace("future.apply", quietly = TRUE)) {
      stop("Parallel bootstrap requires the 'future' and 'future.apply' packages.")
    }
    if (is.null(batch_size)) {
      batch_size <- max(1L, as.integer(n_workers))
    }
    future::plan(future::multisession, workers = n_workers)
    on.exit(future::plan(future::sequential), add = TRUE)

    batches <- split(seq_len(n_boot), ceiling(seq_len(n_boot) / batch_size))
    for (i in seq_along(batches)) {
      batch_ids <- batches[[i]]
      if (verbose) {
        cat("Parallel batch", i, "of", length(batches), "\n")
      }
      batch_res <- future.apply::future_lapply(batch_ids, runner, future.seed = TRUE)
      for (j in seq_along(batch_ids)) {
        results[[batch_ids[j]]] <- batch_res[[j]]
      }
      save_checkpoint(max(batch_ids))
      if (verbose) {
        done <- max(batch_ids)
        elapsed_min <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
        cat("Completed", done, "of", n_boot, "meta bootstraps in", round(elapsed_min, 2), "minutes\n")
      }
    }
  } else {
    for (b in seq_len(n_boot)) {
      if (verbose) {
        cat("Meta bootstrap", b, "of", n_boot, "\n")
      }
      results[[b]] <- runner(b)
      save_checkpoint(b)
      if (verbose) {
        elapsed_min <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
        cat("Completed", b, "of", n_boot, "meta bootstraps in", round(elapsed_min, 2), "minutes\n")
      }
    }
  }

  out <- rbindlist(results, fill = TRUE)
  if (!is.null(checkpoint_file)) {
    saveRDS(out, checkpoint_file)
  }
  out[]
}

summarize_static_intervention_bootstrap <- function(boot_res) {
  br <- as.data.table(copy(boot_res))
  br <- br[status == "success"]
  coef_cols <- setdiff(names(br), c("boot", "status"))

  rbindlist(lapply(coef_cols, function(col) {
    vals <- br[[col]]
    data.table(
      term = sub("^coef_", "", col),
      estimate = mean(vals, na.rm = TRUE),
      lower = stats::quantile(vals, probs = 0.025, na.rm = TRUE, names = FALSE),
      upper = stats::quantile(vals, probs = 0.975, na.rm = TRUE, names = FALSE)
    )
  }))
}

diagnose_joint_bootstrap <- function(boot_obj, config = NULL, tolerance = 1e-12) {
  boot_res <- if (is.list(boot_obj) && !data.table::is.data.table(boot_obj)) boot_obj$results else boot_obj
  boot_res <- as.data.table(copy(boot_res))
  success_rate <- mean(boot_res$status == "success")

  out <- list(
    bootstrap_summary = data.table(
      n_boot = nrow(boot_res),
      n_success = sum(boot_res$status == "success"),
      success_rate = success_rate,
      any_missing_estimates = any(!is.finite(boot_res$risk_natural) | !is.finite(boot_res$risk_intervention) | !is.finite(boot_res$RD) | !is.finite(boot_res$RR), na.rm = TRUE)
    )
  )

  if (is.list(boot_obj) && !is.null(boot_obj$audit_objects)) {
    audit <- boot_obj$audit_objects
    audit_checks <- rbindlist(lapply(seq_along(audit), function(i) {
      obj <- audit[[i]]
      if (is.null(obj) || !is.null(obj$error)) {
        return(data.table(boot = i, ok = FALSE, reason = if (is.null(obj)) 'missing audit object' else obj$error))
      }
      mc <- as.data.table(copy(obj$mc_data))
      nat <- as.data.table(copy(obj$natural))
      int <- as.data.table(copy(obj$intervention))
      uid_col <- if ('uid' %in% names(mc)) 'uid' else config$id
      dupes <- mc[, .N, by = c(uid_col, config$time_in)][N > 1L, .N]
      balanced <- {
        counts <- mc[, .N, by = uid_col]$N
        length(unique(counts)) == 1L
      }
      nat_monotone <- nat[order(uid, get(config$time_in)), all(diff(cumrisk) >= -tolerance), by = uid][, all(V1)]
      int_monotone <- int[order(uid, get(config$time_in)), all(diff(cumrisk) >= -tolerance), by = uid][, all(V1)]
      data.table(
        boot = i,
        ok = dupes == 0 && balanced && nat_monotone && int_monotone,
        duplicate_id_time_rows = dupes,
        balanced_mc_rows = balanced,
        natural_monotone = nat_monotone,
        intervention_monotone = int_monotone,
        reason = NA_character_
      )
    }), fill = TRUE)
    out$audit_checks <- audit_checks
  }

  out
}

diagnose_positivity_joint_gformula <- function(observed_data = NULL,
                                               config,
                                               natural = NULL,
                                               intervention = NULL,
                                               gformula_obj = NULL,
                                               model_fit = NULL) {
  if (!is.null(gformula_obj)) {
    if (is.null(observed_data) && !is.null(gformula_obj$observed_data)) observed_data <- gformula_obj$observed_data
    if (is.null(natural) && !is.null(gformula_obj$natural)) natural <- gformula_obj$natural
    if (is.null(intervention) && !is.null(gformula_obj$intervention)) intervention <- gformula_obj$intervention
    if (is.null(model_fit) && !is.null(gformula_obj$model_fit)) model_fit <- gformula_obj$model_fit
  }

  if (is.null(observed_data) || is.null(natural) || is.null(intervention)) {
    stop('Please supply observed_data, natural, and intervention, or pass a run_joint_gformula result via gformula_obj.')
  }

  obs <- as.data.table(as.data.frame(observed_data))
  nat <- as.data.table(as.data.frame(natural))
  int <- as.data.table(as.data.frame(intervention))

  overall_support <- rbindlist(lapply(config$exposures, function(v) {
    obs_min <- min(obs[[v]], na.rm = TRUE)
    obs_max <- max(obs[[v]], na.rm = TRUE)
    nat_vals <- nat[[v]]
    int_vals <- int[[v]]
    data.table(
      exposure = v,
      obs_min = obs_min,
      obs_max = obs_max,
      natural_min = min(nat_vals, na.rm = TRUE),
      natural_max = max(nat_vals, na.rm = TRUE),
      intervention_min = min(int_vals, na.rm = TRUE),
      intervention_max = max(int_vals, na.rm = TRUE),
      prop_natural_outside = mean(nat_vals < obs_min | nat_vals > obs_max, na.rm = TRUE),
      prop_intervention_outside = mean(int_vals < obs_min | int_vals > obs_max, na.rm = TRUE)
    )
  }), fill = TRUE)

  time_support <- rbindlist(lapply(config$exposures, function(v) {
    obs_time <- obs[, .(obs_min = min(get(v), na.rm = TRUE), obs_max = max(get(v), na.rm = TRUE)), by = c(config$time_in)]
    nat_time <- nat[, .(natural_min = min(get(v), na.rm = TRUE), natural_max = max(get(v), na.rm = TRUE)), by = c(config$time_in)]
    int_time <- int[, .(intervention_min = min(get(v), na.rm = TRUE), intervention_max = max(get(v), na.rm = TRUE)), by = c(config$time_in)]

    nat_flag <- merge(nat[, .(time_value = get(config$time_in), value = get(v))],
                      obs_time[, .(time_value = get(config$time_in), obs_min, obs_max)],
                      by = 'time_value', all.x = TRUE)
    int_flag <- merge(int[, .(time_value = get(config$time_in), value = get(v))],
                      obs_time[, .(time_value = get(config$time_in), obs_min, obs_max)],
                      by = 'time_value', all.x = TRUE)

    nat_out <- nat_flag[, .(prop_natural_outside_time_support = mean(value < obs_min | value > obs_max, na.rm = TRUE)), by = time_value]
    int_out <- int_flag[, .(prop_intervention_outside_time_support = mean(value < obs_min | value > obs_max, na.rm = TRUE)), by = time_value]

    out <- merge(obs_time, nat_time, by = config$time_in, all = TRUE)
    out <- merge(out, int_time, by = config$time_in, all = TRUE)
    data.table::setnames(nat_out, 'time_value', config$time_in)
    data.table::setnames(int_out, 'time_value', config$time_in)
    out <- merge(out, nat_out, by = config$time_in, all = TRUE)
    out <- merge(out, int_out, by = config$time_in, all = TRUE)
    out[, exposure := v]
    out
  }), fill = TRUE)

  mean_by_time <- rbindlist(lapply(config$exposures, function(v) {
    nat_mean <- nat[, .(mean_natural = mean(get(v), na.rm = TRUE)), by = c(config$time_in)]
    int_mean <- int[, .(mean_intervention = mean(get(v), na.rm = TRUE)), by = c(config$time_in)]
    out <- merge(nat_mean, int_mean, by = config$time_in, all = TRUE)
    out[, `:=`(exposure = v, mean_diff = mean_intervention - mean_natural)]
    out
  }), fill = TRUE)

  overlap_vars <- unique(intersect(c(config$time_out, config$exposures, config$baseline_covariates, config$tvc), names(obs)))
  overlap_summary <- data.table(status = 'not_run', message = 'Insufficient overlap variables', min_fitted = NA_real_, q25_fitted = NA_real_, median_fitted = NA_real_, mean_fitted = NA_real_, q75_fitted = NA_real_, max_fitted = NA_real_, prop_gt_0_9 = NA_real_, prop_lt_0_1 = NA_real_)
  overlap_by_source <- NULL

  if (length(overlap_vars) > 0L) {
    obs_flag <- obs[, c('source' = 'observed', ..overlap_vars)]
    int_flag <- int[, c('source' = 'intervention', ..overlap_vars)]
    check_dt <- rbindlist(list(obs_flag, int_flag), fill = TRUE)
    overlap_fit <- tryCatch({
      stats::glm(
        stats::as.formula(paste0("I(source == 'intervention') ~ ", paste(overlap_vars, collapse = ' + '))),
        data = as.data.frame(check_dt),
        family = stats::binomial()
      )
    }, error = function(e) e)

    if (inherits(overlap_fit, 'error')) {
      overlap_summary[, `:=`(status = 'failed', message = conditionMessage(overlap_fit))]
    } else {
      fitted_vals <- stats::fitted(overlap_fit)
      q <- stats::quantile(fitted_vals, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE, names = FALSE)
      overlap_summary <- data.table(
        status = 'ok',
        message = NA_character_,
        min_fitted = q[1],
        q25_fitted = q[2],
        median_fitted = q[3],
        mean_fitted = mean(fitted_vals, na.rm = TRUE),
        q75_fitted = q[4],
        max_fitted = q[5],
        prop_gt_0_9 = mean(fitted_vals > 0.9, na.rm = TRUE),
        prop_lt_0_1 = mean(fitted_vals < 0.1, na.rm = TRUE)
      )
      overlap_by_source <- data.table(source = check_dt$source, fitted = fitted_vals)[,
        .(
          mean_fitted = mean(fitted, na.rm = TRUE),
          q25_fitted = stats::quantile(fitted, 0.25, na.rm = TRUE),
          median_fitted = stats::quantile(fitted, 0.5, na.rm = TRUE),
          q75_fitted = stats::quantile(fitted, 0.75, na.rm = TRUE)
        ),
        by = source
      ]
    }
  }

  link_summary <- NULL
  outcome_model <- NULL
  if (!is.null(model_fit) && !is.null(model_fit$outcome_model)) {
    outcome_model <- model_fit$outcome_model
  }
  if (!is.null(outcome_model)) {
    nat_link <- tryCatch(stats::predict(outcome_model, newdata = nat, type = 'link'), error = function(e) NULL)
    int_link <- tryCatch(stats::predict(outcome_model, newdata = int, type = 'link'), error = function(e) NULL)
    if (!is.null(nat_link) && !is.null(int_link)) {
      link_summary <- rbindlist(list(
        data.table(group = 'natural',
                   min_link = min(nat_link, na.rm = TRUE),
                   q25_link = stats::quantile(nat_link, 0.25, na.rm = TRUE),
                   median_link = stats::quantile(nat_link, 0.5, na.rm = TRUE),
                   mean_link = mean(nat_link, na.rm = TRUE),
                   q75_link = stats::quantile(nat_link, 0.75, na.rm = TRUE),
                   max_link = max(nat_link, na.rm = TRUE)),
        data.table(group = 'intervention',
                   min_link = min(int_link, na.rm = TRUE),
                   q25_link = stats::quantile(int_link, 0.25, na.rm = TRUE),
                   median_link = stats::quantile(int_link, 0.5, na.rm = TRUE),
                   mean_link = mean(int_link, na.rm = TRUE),
                   q75_link = stats::quantile(int_link, 0.75, na.rm = TRUE),
                   max_link = max(int_link, na.rm = TRUE))
      ), fill = TRUE)
    }
  }

  list(
    overall_support = overall_support,
    time_support = time_support,
    mean_by_time = mean_by_time,
    overlap_summary = overlap_summary,
    overlap_by_source = overlap_by_source,
    link_summary = link_summary
  )
}




diagnose_internal_validity_joint_gformula <- function(observed_data = NULL,
                                                      config,
                                                      natural = NULL,
                                                      gformula_obj = NULL,
                                                      use_natural = c("auto", "observed", "modeled", "default"),
                                                      risk_tolerance = 0.01,
                                                      continuous_tolerance = 0.1,
                                                      categorical_tolerance = 0.05) {
  use_natural <- match.arg(use_natural)
  selected_natural <- use_natural

  if (!is.null(gformula_obj)) {
    if (is.null(observed_data) && !is.null(gformula_obj$observed_data)) {
      observed_data <- gformula_obj$observed_data
    }

    if (is.null(natural)) {
      if (use_natural == "observed") {
        natural <- gformula_obj$natural_observed
        selected_natural <- "observed"
      } else if (use_natural == "modeled") {
        natural <- gformula_obj$natural_modeled
        if (is.null(natural)) natural <- gformula_obj$natural
        selected_natural <- "modeled"
      } else if (use_natural == "default") {
        natural <- gformula_obj$natural
        selected_natural <- "default"
      } else {
        natural <- gformula_obj$natural_observed
        selected_natural <- "observed"
        if (is.null(natural)) {
          natural <- gformula_obj$natural_modeled
          selected_natural <- "modeled"
        }
        if (is.null(natural)) {
          natural <- gformula_obj$natural
          selected_natural <- "default"
        }
      }
    }
  }

  if (is.null(observed_data) || is.null(natural)) {
    stop("Please supply observed_data and natural, or pass a run_joint_gformula result via gformula_obj.")
  }

  obs <- as.data.table(as.data.frame(observed_data))
  nat <- as.data.table(as.data.frame(natural))

  summarize_observed_risk <- function(dt) {
    id_name <- config$id
    time_name <- config$time_in
    outcome_name <- config$outcome

    ids <- sort(unique(dt[[id_name]]))
    times <- sort(unique(c(dt[[time_name]], nat[[time_name]])))
    panel <- data.table::CJ(id_value = ids, time_value = times)
    data.table::setnames(panel, c("id_value", "time_value"), c(id_name, time_name))

    obs_y <- dt[, c(id_name, time_name, outcome_name), with = FALSE]
    full <- merge(panel, obs_y, by = c(id_name, time_name), all.x = TRUE)
    full[is.na(get(outcome_name)), (outcome_name) := 0]
    data.table::setorderv(full, c(id_name, time_name))
    full[, cumrisk := 1 - cumprod(1 - get(outcome_name)), by = id_name]
    full[, .(
      mean_risk = mean(cumrisk, na.rm = TRUE),
      mean_risk_per_1000 = 1000 * mean(cumrisk, na.rm = TRUE)
    ), by = c(time_name)]
  }

  observed_risk <- summarize_observed_risk(obs)
  natural_risk <- summarize_risk_trajectory(nat, config)

  risk_trajectory <- merge(
    observed_risk[, .(time_value = get(config$time_in), observed_risk = mean_risk)],
    natural_risk[, .(time_value = get(config$time_in), natural_risk = mean_risk)],
    by = "time_value",
    all = TRUE
  )
  data.table::setnames(risk_trajectory, "time_value", config$time_in)
  risk_trajectory[, `:=`(
    observed_survival = 1 - observed_risk,
    natural_survival = 1 - natural_risk,
    risk_diff = observed_risk - natural_risk,
    abs_risk_diff = abs(observed_risk - natural_risk),
    survival_diff = (1 - observed_risk) - (1 - natural_risk),
    abs_survival_diff = abs((1 - observed_risk) - (1 - natural_risk))
  )]

  survival_trajectory <- risk_trajectory[, .(
    observed_survival,
    natural_survival,
    survival_diff,
    abs_survival_diff
  ), by = c(config$time_in)]

  tvc_vars <- intersect(config$tvc, intersect(names(obs), names(nat)))

  continuous_by_time <- rbindlist(lapply(tvc_vars, function(v) {
    vtype <- tolower(config$covtypes[[v]] %||% "continuous")
    if (vtype %in% c("categorical", "multinomial")) return(NULL)

    obs_sum <- obs[, .(
      observed_mean = mean(get(v), na.rm = TRUE),
      observed_sd = stats::sd(get(v), na.rm = TRUE)
    ), by = c(config$time_in)]
    nat_sum <- nat[, .(
      natural_mean = mean(get(v), na.rm = TRUE),
      natural_sd = stats::sd(get(v), na.rm = TRUE)
    ), by = c(config$time_in)]
    out <- merge(obs_sum, nat_sum, by = config$time_in, all = TRUE)
    out[, `:=`(
      variable = v,
      variable_type = vtype,
      mean_diff = observed_mean - natural_mean,
      abs_mean_diff = abs(observed_mean - natural_mean)
    )]
    out
  }), fill = TRUE)

  continuous_summary <- if (nrow(continuous_by_time) > 0L) {
    continuous_by_time[, .(
      variable_type = variable_type[1],
      rmse_mean = sqrt(mean((observed_mean - natural_mean)^2, na.rm = TRUE)),
      max_abs_mean_diff = max(abs_mean_diff, na.rm = TRUE)
    ), by = variable]
  } else {
    data.table(variable = character(), variable_type = character(), rmse_mean = numeric(), max_abs_mean_diff = numeric())
  }

  categorical_by_time <- rbindlist(lapply(tvc_vars, function(v) {
    vtype <- tolower(config$covtypes[[v]] %||% "continuous")
    if (!vtype %in% c("categorical", "multinomial")) return(NULL)

    obs_cat <- obs[, .N, by = c(config$time_in, v)]
    nat_cat <- nat[, .N, by = c(config$time_in, v)]
    data.table::setnames(obs_cat, v, "level")
    data.table::setnames(nat_cat, v, "level")
    obs_cat[, observed_prop := N / sum(N), by = c(config$time_in)]
    nat_cat[, natural_prop := N / sum(N), by = c(config$time_in)]
    obs_cat <- obs_cat[, .(time_value = get(config$time_in), level, observed_prop)]
    nat_cat <- nat_cat[, .(time_value = get(config$time_in), level, natural_prop)]
    out <- merge(obs_cat, nat_cat, by = c("time_value", "level"), all = TRUE)
    data.table::setnames(out, "time_value", config$time_in)
    out[is.na(observed_prop), observed_prop := 0]
    out[is.na(natural_prop), natural_prop := 0]
    out[, `:=`(
      variable = v,
      variable_type = vtype,
      prop_diff = observed_prop - natural_prop,
      abs_prop_diff = abs(observed_prop - natural_prop)
    )]
    out
  }), fill = TRUE)

  categorical_summary <- if (nrow(categorical_by_time) > 0L) {
    categorical_by_time[, .(
      variable_type = variable_type[1],
      mean_abs_prop_diff = mean(abs_prop_diff, na.rm = TRUE),
      max_abs_prop_diff = max(abs_prop_diff, na.rm = TRUE)
    ), by = variable]
  } else {
    data.table(variable = character(), variable_type = character(), mean_abs_prop_diff = numeric(), max_abs_prop_diff = numeric())
  }

  final_abs_risk_diff <- if (nrow(risk_trajectory) > 0L) tail(risk_trajectory[order(get(config$time_in))]$abs_risk_diff, 1L) else NA_real_
  max_risk_gap <- if (nrow(risk_trajectory) > 0L) max(risk_trajectory$abs_risk_diff, na.rm = TRUE) else NA_real_

  flags <- rbindlist(list(
    data.table(
      check = "final_risk_gap",
      ok = is.finite(final_abs_risk_diff) && final_abs_risk_diff <= risk_tolerance,
      value = final_abs_risk_diff,
      threshold = risk_tolerance,
      detail = "Absolute observed-vs-natural cumulative risk gap at the final time point"
    ),
    data.table(
      check = "max_risk_gap",
      ok = is.finite(max_risk_gap) && max_risk_gap <= risk_tolerance,
      value = max_risk_gap,
      threshold = risk_tolerance,
      detail = "Maximum absolute observed-vs-natural cumulative risk gap across time"
    )
  ), fill = TRUE)

  if (nrow(continuous_summary) > 0L) {
    flags <- rbind(
      flags,
      continuous_summary[, .(
        check = paste0("continuous_tvc_rmse:", variable),
        ok = is.finite(rmse_mean) & rmse_mean <= continuous_tolerance,
        value = rmse_mean,
        threshold = continuous_tolerance,
        detail = "RMSE of observed-vs-natural mean trajectory"
      )],
      fill = TRUE
    )
  }

  if (nrow(categorical_summary) > 0L) {
    flags <- rbind(
      flags,
      categorical_summary[, .(
        check = paste0("categorical_tvc_max_abs_diff:", variable),
        ok = is.finite(max_abs_prop_diff) & max_abs_prop_diff <= categorical_tolerance,
        value = max_abs_prop_diff,
        threshold = categorical_tolerance,
        detail = "Maximum absolute observed-vs-natural category proportion gap"
      )],
      fill = TRUE
    )
  }

  list(
    selected_natural = selected_natural,
    risk_trajectory = risk_trajectory,
    survival_trajectory = survival_trajectory,
    continuous_tvc_by_time = continuous_by_time,
    continuous_tvc_summary = continuous_summary,
    categorical_tvc_by_time = categorical_by_time,
    categorical_tvc_summary = categorical_summary,
    flags = flags
  )
}







