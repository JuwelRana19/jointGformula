\name{diagnose_internal_validity_joint_gformula}
\alias{diagnose_internal_validity_joint_gformula}
\title{Diagnose Internal Validity for jointGformula}
\description{
Compare the observed risk trajectory with the natural-course risk trajectory and summarize
observed-vs-natural agreement for time-varying covariates.
}
\usage{
diagnose_internal_validity_joint_gformula(
  observed_data = NULL,
  config,
  natural = NULL,
  gformula_obj = NULL,
  use_natural = c("auto", "observed", "modeled", "default"),
  risk_tolerance = 0.01,
  continuous_tolerance = 0.1,
  categorical_tolerance = 0.05
)
}
\arguments{
  \item{observed_data}{Observed longitudinal data. Optional when \code{gformula_obj} is supplied.}
  \item{config}{Configuration object created for \code{jointGformula}.}
  \item{natural}{Natural-course prediction dataset. Optional when \code{gformula_obj} is supplied.}
  \item{gformula_obj}{Result from \code{run_joint_gformula()}. If supplied, the helper can select the natural-course output automatically.}
  \item{use_natural}{Which natural-course output to use when \code{gformula_obj} is supplied.}
  \item{risk_tolerance}{Threshold used to flag observed-vs-natural cumulative risk gaps.}
  \item{continuous_tolerance}{Threshold used to flag RMSE in continuous or binary TVC mean trajectories.}
  \item{categorical_tolerance}{Threshold used to flag maximum absolute differences in categorical TVC proportions.}
}
\value{
A list containing the observed-vs-natural risk trajectory, survival trajectory,
continuous and categorical TVC comparison tables, and a set of diagnostic flags.
}
