\name{cfboost}
\alias{cfboost}
\alias{cfboost.formula}
\alias{cfboost_fit}
\title{ Likelihood Based Boosting for Flexible, Structured Survival
  Models }
\description{
  Likelihood-based boosting algorithm to fit flexible,
  structured survival models with component-wise linear or P-spline
  base-learners. Variable selection and model choice are built in
  features.
}
\usage{
cfboost(x, \dots)
\method{cfboost}{formula}(formula, data = list(), weights = NULL,
   na.action = na.omit, control = boost_control(), \dots)
}

\arguments{
  \item{x}{ an object to be pased to \code{cfboost.formula}. }
  \item{formula}{ \code{\link{formula}} of the model to be fit.
    The response \code{y} must be of class \code{\link[survival]{Surv}}. The model
    is specified using the base-learner \code{\link{bols}} and
    \code{\link{bbs}}.}
  \item{data}{ data frame. Contains the data of the model that is to be fitted. }
  \item{weights}{ integer vector. Optional weights of the observations in
    \code{data}. Currently only weights 0 or 1 are supported. See Details. }
  \item{na.action}{ a function that defines how to handle \code{NA}s in
    the data. }
  \item{control}{ an object of class \code{\link{boost_control}}. }
  \item{\dots}{ further arguments to be passed to subfunctions. }
}
\details{
  A structured, flexible Cox-type survival model is fitted to the data
  by likelihood-based boosting. The (component-wise) base-learners are
  specified via the \code{formula}. Examples for the specification of
  base-learners can be found either in \code{\link{bols}} or in
  \code{\link{bbs}}.

  The \code{weights} (a.t.m.) can only be used to specify a learning
  sample which consists of observations with \code{weights == 1} and
  and an out-of-bag sample with \code{weights == 0}. The latter can for
  example be used to determine the appropriate stopping iteration of the
  algorithm.
}
\value{
  An object of class \code{cfboost} is returned:
  \item{data}{data object after data pre-processing (a list of class \code{boost_data}). }
  \item{ensemble}{vector of selected base-learners. }
  \item{ensembless}{list of coefficient estimates for each iteration.
  Only the coefficient of the selected base-learner is stored; i.e.,
  \code{ensembless[[i]]} is the coefficient vector of base-learner
  \code{ensemble[i]} (in iteration \code{i}). }
  \item{fit}{vector of fitted values (in the last iteration).}
  \item{offset}{offset value.}
  \item{control}{control parameters as specified.}
  \item{response}{the specified response variable (an object of class \code{\link[survival]{Surv}}).}
  \item{risk}{vector of the risk (negative log likelihood) computed on
  the training or validation sample (see \code{\link{boost_control}}).}
  \item{weights}{weights used for model fitting.}
  \item{df}{matrix of estimated degrees of freedom for smooth base-learners.}
  \item{coefs}{(list of) estimated coefficients (in the final boosting iteration).}
  \item{predict}{function for prediction (to be accessed via
    \code{\link[=predict.cfboost]{predict}}).}
  \item{call}{the model call.}
}
\references{
  ~~ PUT TECH.REP. HERE ~~
}
\author{ Benjamin Hofner }
\seealso{ \code{\link{boost_control}} for the control arguments for
  \code{cfboost} and \code{\link{bols}} or \code{\link{bbs}} for the
  available base-learners.
  Methods to process the output can be found in
  \code{\link[CoxFlexBoost]{methods}}. A complete example can be found in
  \code{\link{CoxFlexBoost-package}}.
}
\examples{
## a simple example (without time-varying effects)

set.seed(1234)

## sample covariates first
X <- matrix(NA, nrow=400, ncol=3)
X[,1] <- runif(400, -1, 1)
X[,2] <- runif(400, -1, 1)
X[,3] <- runif(400, -1, 1)

## time-dependent hazard rate
lambda <- function(time, x){
  exp(0 * time + 0.7 * x[1] + x[2]^2)
}

## specify censoring function
cens_fct <- function(time, mean_cens){
  censor_time <- rexp(n = length(time), rate = 1/mean_cens)
  event <- (time <= censor_time)
  t_obs <- apply(cbind(time, censor_time), 1, min)
  return(cbind(t_obs, event))
}
data <- rSurvTime(lambda, X, cens_fct, mean_cens = 5)

ctrl <- boost_control( mstop = 100, risk="oobag")
weights <- sample(c(0,1), 400, replace=TRUE, prob=c(0.25, 0.75))

## fit (a simple) model
model <- cfboost(Surv(time, event) ~ bbs(x.1) + bbs(x.2) + bbs(x.3),
                 control = ctrl, data = data, weights = weights)

model <- model[mstop(model)]
summary(model)

## fit (a simple) model with model choice
## i.e., with decomposition of base-learners
model_2 <- cfboost(Surv(time, event) ~ bols(x.1) + bbs(x.1, df=1, center=TRUE)
                   + bols(x.2) + bbs(x.2, df=1, center=TRUE)
                   + bols(x.3) + bbs(x.3, df=1, center=TRUE),
                   control = ctrl, data = data, weights = weights)

model_2 <- model_2[mstop(model_2)]
## only bols(x.1) and bbs(x.2, ...)  are selected (as desired)
summary(model_2)

}

\keyword{ survival }
\keyword{ smooth }
\keyword{ models }

