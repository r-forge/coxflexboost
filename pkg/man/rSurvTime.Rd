\name{rSurvTime}
\alias{rSurvTime}

\title{ Simulate Arbitrary Survival Times }
\description{
  Simulates arbitrary (Cox-type) survival times using the inversion
  method. The survival times do \bold{\emph{not} need to follow a known distribution},
  as e.g., the exponential or Weibull distribution. Furthermore,
  time-varying effects can be sampled, i.e., the data does not
  neccessarily need to statisfy the proportinal hazards assumption.
}
\usage{
rSurvTime(lambda, x, cens_fct, upper = 1000, ..., file = NULL)
}

\arguments{
  \item{lambda}{ function. Baseline hazard \eqn{\lambda(t, x)}{lambda(t,x)}
    where time must be first argument. }
  \item{x}{ matrix. (Sampled) values for covariates (without time) }
  \item{cens_fct}{ function. Function to compute (random) censoring }
  \item{upper}{ upper boundary of the interval the random survival times fall into }
  \item{\dots}{ further arguments to be passed to \code{lambda} or \code{cens_fct}. }
  \item{file}{ character. name of the data file the generated data set should be stored into
               (e.g., "survtimes.RData") or NULL if the dataset should directly be
               returned in R }
}
\details{
  By specifying an (arbitrary) hazard rate and a censoring function it
  is possible to sample survival times according to this scheme. The
  concept of Bender et al. (2005) is expanded to even allow for
  time-varying effects in the hazard rate.

  The function \code{lambda} needs to be vectorized in its first
  argument, the time. To achieve this even if no the hazard rate is
  specified as time-constant a nice trick is to spefify \code{0 * time}
  explizitly in lambda (see example for more details).
  Another important remark is that \code{lambda} expects \code{x} to be
  a vector, i.e., we just pass one row of the matrix of coefficients to
  lambda (see example below).

  The censoring function \code{cens_fct} takes as first argument the
  \bold{uncensored} (sampled) survival times. Further arguments can be
  specified (and accessed via the \code{\dots} argument in
  \code{rSurvTime}). The function must return a matrix consisting of the
  \bold{observed} survival times and the non-censoring indicator (see
  example below).
}
\value{
  A \code{data.frame} consisting of the observed survival time
  (\code{time}), the non-censoring indicator (\code{event}) and further
  covariates \code{x} is returned. If \code{file} is specified, the
  data.frame is additionally stored on the disc.
}
\references{
  Ralph Bender and Thomas Augustin and Maria Blettner (2005),
  Generating Survival Times to Simulate Cox Proportional Hazards
  Models.
  \emph{Statistics in Medicine}, \bold{24}, 1713--1723.
}
\author{ Benjamin Hofner }
\note{ Note that the covariate values should not become too large as
  otherwise numerical problems may occure. A linear predictor in the
  range from -3 to 3 has shown to be sufficiently small.
}
\examples{
set.seed(1234)

## sample covariates first
X <- matrix(NA, nrow=400, ncol=3)
X[,1] <- runif(400, -1, 1)
X[,2] <- runif(400, -1, 1)
X[,3] <- runif(400, -1, 1)

## time-constant hazard rate (needs to be a vectororized in the first argument time)
lambda <- function(time, x){
  exp(0 * time + 0.7 * x[1] + x[2]^2 + 0 * x[3])
}

## time-dependent hazard rate
lambda <- function(time, x){
  exp(log(time) + 0.7 * x[1] + x[2]^2 + 0 * x[3])
}

## specify censoring function
cens_fct <- function(time, mean_cens){
  ## censoring times are independent exponentially distributed
  censor_time <- rexp(n = length(time), rate = 1/mean_cens)
  event <- (time <= censor_time)
  t_obs <- apply(cbind(time, censor_time), 1, min)
  ## return matrix of observed survival times and event indicator
  return(cbind(t_obs, event))
}

data <- rSurvTime(lambda, X, cens_fct, mean_cens = 5)

## it is also possible to sample time-varying effects
lambda <-  function(time, x){
  exp(log(time) + 0.7 * x[1] + x[2]^2 + (x[3] * time)/3)
}

data <- rSurvTime(lambda, X, cens_fct, mean_cens = 5)
}

\keyword{ survival }
\keyword{ distribution }

