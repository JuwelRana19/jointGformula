\name{run_joint_gformula}
\alias{run_joint_gformula}
\alias{create_mc_skeleton}
\alias{fit_joint_models}
\alias{simulate_joint_gformula}
\alias{predict_joint_outcome}
\alias{estimate_joint_effects}
\alias{summarize_risk_trajectory}
\title{Run the jointGformula Pipeline}
\description{
Fit the nuisance models, create the Monte Carlo skeleton, simulate the natural course,
apply the joint intervention, and estimate final cumulative risks and contrasts.
}
\details{
When time-varying covariate models use other covariates from the same time point as predictors,
the order in \code{config$tvc} defines the within-time simulation sequence. Set
\code{config$tvc} so earlier variables in that vector are available before later models
are evaluated.
}
\usage{
run_joint_gformula(
  data,
  config,
  mc_size = 10000L,
  seed = 1234L,
  replace_mc = NULL,
  temp = 1.25,
  natural_course = config$natural_course \%||\% "modeled"
)
}
\arguments{
  \item{data}{Observed longitudinal input data.}
  \item{config}{Configuration object created for \code{jointGformula}.}
  \item{mc_size}{Monte Carlo sample size.}
  \item{seed}{Random seed for Monte Carlo resampling and simulation.}
  \item{replace_mc}{Whether to sample subjects with replacement when constructing the Monte Carlo panel.}
  \item{temp}{Sampling temperature for categorical simulation draws.}
  \item{natural_course}{Character string; \code{"modeled"} uses the simulated natural course,
  \code{"observed"} predicts the natural course on the observed data, and
  \code{c("observed", "modeled")} returns both while setting \code{natural}
  to the observed version.}
}
\value{
A list containing \code{mc_data}, \code{sim_data}, \code{natural},
\code{natural_observed}, \code{natural_modeled}, \code{intervention},
\code{estimates}, the fitted models, and the configuration.
}
