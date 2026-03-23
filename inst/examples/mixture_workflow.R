library(jointGformula)

cfg <- make_toronto_mixture_config()
dat <- readRDS("NO2.RDS")
dat <- make_lags(dat, cfg)
validate_gformula_data(dat, cfg)

# Quick run without bootstrap
out <- run_joint_pipeline(
  data = dat,
  config = cfg,
  mc_size = 100,
  seed = 1234,
  bootstrap = FALSE
)

print(out$estimates)

# Start from an intermediate stage if needed
mc_dat <- create_mc_skeleton(dat, cfg, mc_size = 100, seed = 1234)
mods <- fit_joint_models(dat, cfg)
sim_dat <- simulate_joint_gformula(mc_dat, cfg, mods)

natural <- predict_joint_outcome(sim_dat, cfg, mods$outcome_model, use_intervention = FALSE)
# Let predict_joint_outcome() apply the configured intervention directly.
intervention <- predict_joint_outcome(sim_dat, cfg, mods$outcome_model, use_intervention = TRUE)

print(estimate_joint_effects(natural, intervention, cfg))
