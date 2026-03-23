\name{diagnose_positivity_joint_gformula}
\alias{diagnose_positivity_joint_gformula}
\title{Diagnose Positivity and Extrapolation in jointGformula}
\usage{
diagnose_positivity_joint_gformula(
  observed_data = NULL,
  config,
  natural = NULL,
  intervention = NULL,
  gformula_obj = NULL,
  model_fit = NULL
)
}
\arguments{
  \item{observed_data}{Observed longitudinal data.}
  \item{config}{A configuration object created by \code{make_gformula_config()}.}
  \item{natural}{Natural-course simulated data.}
  \item{intervention}{Intervention simulated data.}
  \item{gformula_obj}{Optional output from \code{run_joint_gformula()}, used as a shortcut to supply observed, natural, intervention, and model fit objects.}
  \item{model_fit}{Optional fitted model object from \code{fit_joint_models()} or \code{run_joint_gformula()}. If provided, link-scale summaries are returned.}
}
\value{
A list with overall support checks, time-specific support checks, mean exposure comparisons by time, an overlap-classifier summary, and optional link-scale summaries.
}
\description{
Provides lightweight diagnostics for positivity and extrapolation in the joint g-formula workflow. The helper compares observed, natural-course, and intervention exposure distributions overall and by time, summarizes intervention-versus-natural mean differences, and fits a simple classifier to assess overlap between observed and intervention covariate patterns.
}
\examples{
\dontrun{
fit <- run_joint_gformula(NO2_complete_clean, cfg, mc_size = 100, seed = 1234)
pos_diag <- diagnose_positivity_joint_gformula(gformula_obj = fit, config = cfg)
pos_diag$overall_support
pos_diag$time_support
pos_diag$overlap_summary
}
}
