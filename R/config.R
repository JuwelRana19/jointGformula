make_gformula_config <- function(id,
                                 time_in = NULL,
                                 time_name = NULL,
                                 time_out = NULL,
                                 outcome = NULL,
                                 outcome_name = NULL,
                                 exposures = NULL,
                                 intvars = NULL,
                                 intervention_exposures,
                                 time_varying_covariates = NULL,
                                 covnames = NULL,
                                 covtypes = NULL,
                                 covranges = NULL,
                                 time_fixed_covariates = NULL,
                                 basecovs = NULL,
                                 factor_vars = NULL,
                                 exposure_lags = NULL,
                                 exposure_types = NULL,
                                 exposure_ranges = NULL,
                                 tvc_lags = NULL,
                                 exposure_order = NULL,
                                 outcome_formula = NULL,
                                 ymodel = NULL,
                                 tvc_formulas = NULL,
                                 covmodels = NULL,
                                 exposure_formulas = NULL,
                                 intervention = "natural_vs_intervention",
                                 outcome_type = "survival",
                                 natural_course = "modeled",
                                 panel_mode = "full_panel") {
  normalize_named <- function(x, keys, default = NULL) {
    if (is.null(x)) {
      x <- rep(default, length(keys))
      names(x) <- keys
      return(x)
    }
    if (is.null(names(x))) {
      if (length(x) != length(keys)) {
        stop("Unnamed type vector must have the same length as its variables.")
      }
      names(x) <- keys
      return(x)
    }
    if (!all(keys %in% names(x))) {
      stop("Missing named entries for: ", paste(setdiff(keys, names(x)), collapse = ", "))
    }
    x[keys]
  }

  if (is.null(time_in)) time_in <- time_name
  if (is.null(outcome)) outcome <- outcome_name
  if (is.null(exposures)) exposures <- intvars
  if (is.null(time_varying_covariates)) time_varying_covariates <- covnames
  if (is.null(time_fixed_covariates)) time_fixed_covariates <- basecovs
  if (is.null(outcome_formula)) outcome_formula <- ymodel
  if (is.null(tvc_formulas)) tvc_formulas <- covmodels

  if (is.null(time_in) || is.null(outcome) || is.null(exposures) ||
      is.null(time_varying_covariates) || is.null(time_fixed_covariates) ||
      is.null(outcome_formula) || is.null(tvc_formulas) || is.null(exposure_formulas)) {
    stop("Config is missing required pieces. Supply id, time_in/time_name, outcome/outcome_name, exposures/intvars, time_varying_covariates/covnames, time_fixed_covariates/basecovs, outcome_formula/ymodel, tvc_formulas/covmodels, and exposure_formulas.")
  }

  if (is.null(time_out)) {
    stop("Please supply time_out explicitly.")
  }

  if (is.null(exposure_lags)) {
    exposure_lags <- paste0(exposures, "_lag1")
  }
  if (is.null(tvc_lags)) {
    tvc_lags <- paste0(time_varying_covariates, "_lag1")
  }
  if (is.null(exposure_order)) {
    exposure_order <- exposures
  }

  covtypes <- normalize_named(covtypes, time_varying_covariates, default = "categorical")
  exposure_types <- normalize_named(exposure_types, exposures, default = "lognormal")

  stopifnot(length(exposures) == length(intervention_exposures))
  stopifnot(length(exposures) == length(exposure_lags))
  stopifnot(length(time_varying_covariates) == length(tvc_lags))
  stopifnot(all(exposure_order %in% exposures))
  stopifnot(all(names(tvc_formulas) %in% time_varying_covariates))
  stopifnot(all(names(exposure_formulas) %in% exposures))

  list(
    id = id,
    time_in = time_in,
    time_name = time_in,
    time_out = time_out,
    outcome = outcome,
    outcome_name = outcome,
    exposures = exposures,
    intvars = exposures,
    intervention_exposures = intervention_exposures,
    exposure_lags = exposure_lags,
    exposure_types = exposure_types,
    exposure_ranges = exposure_ranges,
    time_varying_covariates = time_varying_covariates,
    covnames = time_varying_covariates,
    covtypes = covtypes,
    covranges = covranges,
    time_fixed_covariates = time_fixed_covariates,
    basecovs = time_fixed_covariates,
    baseline_covariates = time_fixed_covariates,
    tvc = time_varying_covariates,
    tvc_lags = tvc_lags,
    factor_vars = factor_vars,
    exposure_order = exposure_order,
    outcome_formula = outcome_formula,
    ymodel = outcome_formula,
    tvc_formulas = tvc_formulas,
    covmodels = tvc_formulas,
    exposure_formulas = exposure_formulas,
    intervention = intervention,
    outcome_type = outcome_type,
    natural_course = natural_course,
    panel_mode = panel_mode
  )
}


make_toronto_mixture_config <- function(id = "UniqID",
                                        time_in = "TimeInn",
                                        time_out = "TimeOut",
                                        outcome = "status",
                                        exposures = c("NO2_RF", "BC_RF", "PM1_RF"),
                                        intervention_exposures = c("NO2_RF_NoTrucks", "BC_RF_NoTrucks", "PM1_RF_NoTrucks"),
                                        time_varying_covariates = c("CSize", "Urban_form", "ethniconcentration", "deprivation", "dependency", "instability"),
                                        time_fixed_covariates = c(
                                          "age", "sex", "visible_minority", "IndigenousIdentity",
                                          "landed_migration", "marsth", "Education", "employment",
                                          "Occupation", "income_inadequacy"
                                        )) {
  make_gformula_config(
    id = id,
    time_in = time_in,
    time_out = time_out,
    outcome = outcome,
    exposures = exposures,
    intervention_exposures = intervention_exposures,
    time_varying_covariates = time_varying_covariates,
    covtypes = stats::setNames(rep("categorical", length(time_varying_covariates)), time_varying_covariates),
    time_fixed_covariates = time_fixed_covariates,
    exposure_lags = paste0(exposures, "_lag1"),
    exposure_types = stats::setNames(rep("lognormal", length(exposures)), exposures),
    tvc_lags = paste0(time_varying_covariates, "_lag1"),
    exposure_order = c("BC_RF", "PM1_RF", "NO2_RF"),
    outcome_formula = status ~
      NO2_RF + BC_RF + PM1_RF +
      NO2_RF:TimeOut + BC_RF:TimeOut + PM1_RF:TimeOut +
      TimeOut + I(TimeOut^2) +
      age + sex + visible_minority + IndigenousIdentity + landed_migration +
      marsth + Education + employment + Occupation + income_inadequacy +
      CSize + Urban_form + ethniconcentration + deprivation + dependency + instability,
    tvc_formulas = list(
      CSize = CSize ~ TimeOut + I(TimeOut^2) +
        age + sex + visible_minority + IndigenousIdentity + landed_migration +
        marsth + Education + employment + Occupation + income_inadequacy +
        CSize_lag1 + deprivation_lag1 + dependency_lag1 + instability_lag1 + ethniconcentration_lag1,
      Urban_form = Urban_form ~ TimeOut + I(TimeOut^2) +
        age + sex + visible_minority + IndigenousIdentity + landed_migration +
        marsth + Education + employment + Occupation + income_inadequacy +
        Urban_form_lag1 + CSize_lag1 + CSize +
        deprivation_lag1 + dependency_lag1 + instability_lag1 + ethniconcentration_lag1,
      ethniconcentration = ethniconcentration ~ TimeOut + I(TimeOut^2) +
        age + sex + visible_minority + IndigenousIdentity + landed_migration +
        marsth + Education + employment + Occupation + income_inadequacy +
        ethniconcentration_lag1 + CSize + Urban_form +
        deprivation_lag1 + dependency_lag1 + instability_lag1,
      deprivation = deprivation ~ TimeOut + I(TimeOut^2) +
        age + sex + visible_minority + IndigenousIdentity + landed_migration +
        marsth + Education + employment + Occupation + income_inadequacy +
        deprivation_lag1 + CSize + Urban_form + ethniconcentration +
        dependency_lag1 + instability_lag1,
      dependency = dependency ~ TimeOut + I(TimeOut^2) +
        age + sex + visible_minority + IndigenousIdentity + landed_migration +
        marsth + Education + employment + Occupation + income_inadequacy +
        dependency_lag1 + CSize + Urban_form + ethniconcentration + deprivation +
        instability_lag1,
      instability = instability ~ TimeOut + I(TimeOut^2) +
        age + sex + visible_minority + IndigenousIdentity + landed_migration +
        marsth + Education + employment + Occupation + income_inadequacy +
        instability_lag1 + CSize + Urban_form + ethniconcentration + deprivation + dependency
    ),
    exposure_formulas = list(
      BC_RF = log(BC_RF) ~ BC_RF_lag1 + BC_RF_lag1:TimeOut +
        TimeOut + I(TimeOut^2) +
        age + sex + visible_minority + IndigenousIdentity + landed_migration +
        marsth + Education + employment + Occupation + income_inadequacy +
        CSize + Urban_form + ethniconcentration + deprivation + dependency + instability,
      PM1_RF = log(PM1_RF) ~ PM1_RF_lag1 + PM1_RF_lag1:TimeOut +
        BC_RF + BC_RF_lag1 +
        TimeOut + I(TimeOut^2) +
        age + sex + visible_minority + IndigenousIdentity + landed_migration +
        marsth + Education + employment + Occupation + income_inadequacy +
        CSize + Urban_form + ethniconcentration + deprivation + dependency + instability,
      NO2_RF = log(NO2_RF) ~ NO2_RF_lag1 + NO2_RF_lag1:TimeOut +
        BC_RF + PM1_RF + BC_RF_lag1 + PM1_RF_lag1 +
        TimeOut + I(TimeOut^2) +
        age + sex + visible_minority + IndigenousIdentity + landed_migration +
        marsth + Education + employment + Occupation + income_inadequacy +
        CSize + Urban_form + ethniconcentration + deprivation + dependency + instability
    )
  )
}
