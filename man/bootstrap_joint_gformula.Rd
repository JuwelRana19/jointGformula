\name{bootstrap_joint_gformula}
\alias{bootstrap_joint_gformula}
\alias{summarize_joint_bootstrap}
\alias{diagnose_joint_bootstrap}
\title{Bootstrap jointGformula Fits}
\description{
Run a subject-level cluster bootstrap for \code{jointGformula}. In regular mode the function returns
replicate-level risk estimates. In audit mode it also stores the per-replicate bootstrap data,
Monte Carlo data, simulated data, natural-course output, intervention output, and estimated effects.
}
\usage{
bootstrap_joint_gformula(
  data,
  config,
  n_boot = 200L,
  mc_size = 10000L,
  seed = 1234L,
  replace_mc = NULL,
  temp = 1.25,
  natural_course = config$natural_course \%||\% "modeled",
  parallel = FALSE,
  n_workers = 2L,
  batch_size = NULL,
  checkpoint_file = NULL,
  verbose = TRUE,
  audit_mode = FALSE,
  keep_boot_objects = FALSE
)

summarize_joint_bootstrap(boot_res)

diagnose_joint_bootstrap(boot_obj, config = NULL, tolerance = 1e-12)
}
\arguments{
  \item{data, config, mc_size, seed, replace_mc, temp, natural_course}{See \code{run_joint_gformula()}.}
  \item{n_boot}{Number of bootstrap replicates.}
  \item{parallel}{Logical; whether to evaluate bootstrap replicates in parallel.}
  \item{n_workers}{Number of workers when \code{parallel = TRUE}.}
  \item{batch_size}{Number of replicates per parallel batch.}
  \item{checkpoint_file}{Optional RDS file to update during the bootstrap run.}
  \item{verbose}{Logical; print bootstrap progress.}
  \item{audit_mode}{Logical; if \code{TRUE}, keep replicate objects for auditing.}
  \item{keep_boot_objects}{Logical; explicit switch to retain audit objects.}
  \item{boot_res}{Result from \code{bootstrap_joint_gformula()}.}
  \item{boot_obj}{Bootstrap result or audit object to diagnose.}
  \item{tolerance}{Numerical tolerance for monotonicity checks.}
}
\value{
\code{bootstrap_joint_gformula()} returns either a data table of bootstrap risk estimates or,
in audit mode, a list with \code{results}, \code{audit_objects}, and \code{audit_mode}.
\code{summarize_joint_bootstrap()} returns a summary table of bootstrap means and percentile intervals.
\code{diagnose_joint_bootstrap()} returns summary checks for bootstrap success and, when audit objects are available,
replicate-level structural validation checks.
}
