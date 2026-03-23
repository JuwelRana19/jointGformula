\name{run_joint_pipeline}
\alias{run_joint_pipeline}
\title{Convenience Wrapper for jointGformula}
\description{
Run either the main \code{jointGformula} pipeline or its bootstrap wrapper using a single entry point.
}
\usage{
run_joint_pipeline(
  data,
  config,
  mc_size = 10000L,
  seed = 1234L,
  replace_mc = NULL,
  temp = 1.25,
  bootstrap = FALSE,
  natural_course = config$natural_course \%||\% "modeled",
  n_boot = 200L,
  parallel = FALSE,
  n_workers = 2L,
  batch_size = NULL,
  checkpoint_file = NULL,
  verbose = TRUE
)
}
\arguments{
  \item{data, config, mc_size, seed, replace_mc, temp, natural_course}{See \code{run_joint_gformula()}.}
  \item{bootstrap}{Logical; if \code{TRUE}, run the bootstrap workflow instead of the single-run workflow.}
  \item{n_boot, parallel, n_workers, batch_size, checkpoint_file, verbose}{See \code{bootstrap_joint_gformula()}.}
}
\value{
If \code{bootstrap = FALSE}, returns the result of \code{run_joint_gformula()}.
If \code{bootstrap = TRUE}, returns a list with bootstrap results, bootstrap summary, and the configuration.
}
