# jointGformula

Small internal R package for config-driven parametric g-formula analyses.

## What it gives you

- One reusable engine for joint exposure g-formula work
- Dataset-style configuration similar in spirit to `gfoRmula`
- Flexible per-variable model types for time-varying covariates and exposures
- A quick wrapper with `bootstrap = TRUE/FALSE`
- Stage-level functions so you can start at any point in the workflow
- A Toronto mixture helper as a preset, not a requirement

## Main user-facing functions

- `make_gformula_config()`
- `make_toronto_mixture_config()`
- `validate_gformula_data()`
- `make_lags()`
- `create_mc_skeleton()`
- `fit_joint_models()`
- `simulate_joint_gformula()`
- `assign_joint_intervention()`
- `predict_joint_outcome()`
- `estimate_joint_effects()`
- `run_joint_gformula()`
- `run_joint_pipeline()`
- `bootstrap_joint_gformula()`
- `summarize_joint_bootstrap()`

## Build offline

From a machine with R and build tools:

```r
system("R CMD build jointGformula")
```

Then move the resulting tarball to RDC or another offline machine and install with:

```r
install.packages("jointGformula_0.1.0.tar.gz", repos = NULL, type = "source")
```

## Generic configuration

Prepare factors and transformations in the dataset before calling the package, then define the analysis structure in the config.

```r
library(jointGformula)

cfg <- make_gformula_config(
  id = "UniqID",
  time_name = "TimeInn",
  time_out = "TimeOut",
  outcome_name = "status",
  intvars = c("NO2_RF", "BC_RF", "PM1_RF"),
  intervention_exposures = c("NO2_RF_NoTrucks", "BC_RF_NoTrucks", "PM1_RF_NoTrucks"),
  covnames = c("CSize", "Urban_form", "ethniconcentration", "deprivation", "dependency", "instability"),
  covtypes = c(
    CSize = "categorical",
    Urban_form = "categorical",
    ethniconcentration = "categorical",
    deprivation = "categorical",
    dependency = "categorical",
    instability = "categorical"
  ),
  basecovs = c("age", "sex", "visible_minority", "IndigenousIdentity", "landed_migration", "marsth", "Education", "employment", "Occupation", "income_inadequacy"),
  exposure_types = c(NO2_RF = "lognormal", BC_RF = "lognormal", PM1_RF = "lognormal"),
  ymodel = status ~ NO2_RF + BC_RF + PM1_RF + NO2_RF:TimeOut + BC_RF:TimeOut + PM1_RF:TimeOut + TimeOut + I(TimeOut^2) + age + sex + visible_minority + IndigenousIdentity + landed_migration + marsth + Education + employment + Occupation + income_inadequacy + CSize + Urban_form + ethniconcentration + deprivation + dependency + instability,
  covmodels = list(
    CSize = CSize ~ TimeOut + I(TimeOut^2) + age + sex + visible_minority + IndigenousIdentity + landed_migration + marsth + Education + employment + Occupation + income_inadequacy + CSize_lag1 + deprivation_lag1 + dependency_lag1 + instability_lag1 + ethniconcentration_lag1,
    Urban_form = Urban_form ~ TimeOut + I(TimeOut^2) + age + sex + visible_minority + IndigenousIdentity + landed_migration + marsth + Education + employment + Occupation + income_inadequacy + Urban_form_lag1 + CSize_lag1 + CSize + deprivation_lag1 + dependency_lag1 + instability_lag1 + ethniconcentration_lag1,
    ethniconcentration = ethniconcentration ~ TimeOut + I(TimeOut^2) + age + sex + visible_minority + IndigenousIdentity + landed_migration + marsth + Education + employment + Occupation + income_inadequacy + ethniconcentration_lag1 + CSize + Urban_form + deprivation_lag1 + dependency_lag1 + instability_lag1,
    deprivation = deprivation ~ TimeOut + I(TimeOut^2) + age + sex + visible_minority + IndigenousIdentity + landed_migration + marsth + Education + employment + Occupation + income_inadequacy + deprivation_lag1 + CSize + Urban_form + ethniconcentration + dependency_lag1 + instability_lag1,
    dependency = dependency ~ TimeOut + I(TimeOut^2) + age + sex + visible_minority + IndigenousIdentity + landed_migration + marsth + Education + employment + Occupation + income_inadequacy + dependency_lag1 + CSize + Urban_form + ethniconcentration + deprivation + instability_lag1,
    instability = instability ~ TimeOut + I(TimeOut^2) + age + sex + visible_minority + IndigenousIdentity + landed_migration + marsth + Education + employment + Occupation + income_inadequacy + instability_lag1 + CSize + Urban_form + ethniconcentration + deprivation + dependency
  ),
  exposure_formulas = list(
    BC_RF = log(BC_RF) ~ BC_RF_lag1 + BC_RF_lag1:TimeOut + TimeOut + I(TimeOut^2) + age + sex + visible_minority + IndigenousIdentity + landed_migration + marsth + Education + employment + Occupation + income_inadequacy + CSize + Urban_form + ethniconcentration + deprivation + dependency + instability,
    PM1_RF = log(PM1_RF) ~ PM1_RF_lag1 + PM1_RF_lag1:TimeOut + BC_RF + BC_RF_lag1 + TimeOut + I(TimeOut^2) + age + sex + visible_minority + IndigenousIdentity + landed_migration + marsth + Education + employment + Occupation + income_inadequacy + CSize + Urban_form + ethniconcentration + deprivation + dependency + instability,
    NO2_RF = log(NO2_RF) ~ NO2_RF_lag1 + NO2_RF_lag1:TimeOut + BC_RF + PM1_RF + BC_RF_lag1 + PM1_RF_lag1 + TimeOut + I(TimeOut^2) + age + sex + visible_minority + IndigenousIdentity + landed_migration + marsth + Education + employment + Occupation + income_inadequacy + CSize + Urban_form + ethniconcentration + deprivation + dependency + instability
  )
)
```

## Quick run

```r
cfg <- make_toronto_mixture_config()
dat <- readRDS("NO2.RDS")
validate_gformula_data(dat, cfg)

out <- run_joint_pipeline(
  data = dat,
  config = cfg,
  mc_size = 1000,
  seed = 1234,
  bootstrap = FALSE
)

out$estimates
```

## Supported variable types

For `covtypes`, the package currently supports:

- `categorical`
- `binary`
- `normal`
- `bounded normal`

For `exposure_types`, the package currently supports:

- `lognormal`
- `normal`
- `bounded normal`
- `binary`

Use `covranges` or `exposure_ranges` when you want bounded variables clipped to user-specified limits.
## Additions

- `describe_model_plan()` prints the planned outcome, TVC, and exposure models with fitting method and formula.
- Parallel bootstrap now saves checkpoint files after each completed batch.
- `batch_size` controls how many bootstrap replicates are processed before each parallel checkpoint save.

## Marginal cumulative incidence

For survival outcomes, `predict_joint_outcome()` returns:

- `hazard`: predicted conditional probability of the outcome at each time
- `cumrisk`: cumulative incidence built from those hazards within subject over time

Use `summarize_risk_trajectory()` to get the mean risk curve over time from `cumrisk`. The returned natural and intervention datasets keep only the active exposure columns used for prediction.
