\name{CoxFlexBoost-package}
\alias{CoxFlexBoost-package}
\alias{CoxFlexBoost}
\alias{coxflexboost}
\docType{package}
\title{
CoxFlexBoost: Boosting Flexible Cox Models (with Time-Varying Effects)
}
\description{
Likelihood-based boosting approach to fit flexible,
structured survival models with component-wise linear or P-spline
base-learners. Variable selection and model choice are built in
features.
}
\details{
\tabular{ll}{
Package: \tab CoxFlexBoost\cr
Type: \tab Package\cr
Version: \tab 0.7-0\cr
Date: \tab 2009-03-XX\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}
}
\author{
  Benjamin Hofner <firstname.lastname@imbe.med.uni-erlangen.de>
}

\seealso{
  Important functions include: \code{\link{rSurvTime}} for simulations
  of survival times,
  \code{\link{cfboost}} for the boosting algorithm, \code{\link{bols}}
  and \code{\link{bbs}} for the base-learners and
  \code{\link[CoxFlexBoost]{methods}} for methods to be used with
  results.
}

\references{
  Benjamin Hofner, Torsten Hothorn and Thomas Kneib (2008),
  Variable Selection and Model Choice in Structured Survival Models.
  \emph{Department of Statistics, Technical Report No. 43}.
  \url{http://epub.ub.uni-muenchen.de/7901/}
}

\examples{
## a simple example

#################
# simulate data #
#################

set.seed(6102)

## sample covariates first
X <- matrix(NA, nrow=400, ncol=3)
X[,1] <- runif(400, -1, 1)
X[,2] <- runif(400, -1, 1)
X[,3] <- runif(400, -1, 1)

## time-constant hazard rate
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

##################
# estimate model #
##################

ctrl <- boost_control( mstop = 100, risk="oobag")
weights <- c(rep(1, 300), rep(0, 100))

## fit (a simple) model
model <- cfboost(Surv(time, event) ~ bbs(x.1) + bbs(x.2) + bbs(x.3),
                 control = ctrl, data = data, weights = weights)

#########################
# processing the output #
#########################

## determine appropriate stopping iteration
mstop(model)

## summary (with mstop)
summary(model)
## or (almost) same result
model

## plot the risk
plot(risk(model))

## proceed with optimal model
model_opt <- model[mstop(model)]

## extract coefficients of optimal model
## (not very meaningfull for P-spline base-learner)
coef(model_opt)

## plot optimal model
plot(model_opt)

## extract fitted values
fitted(model_opt)

## make prediction for the first out-of-bag observation
predict(model_opt, newdata = data[301,])

## look at the selection frequencies for optimal model)
freq.sel(model_opt)

###################################
# model with time-varying effects #
###################################

\dontrun{
## time-dependent hazard rate
lambda <- function(time, x){
  exp(log(time) + 0.7 * x[1] + x[2]^2)
}

data <- rSurvTime(lambda, X, cens_fct, mean_cens = 5)

## fit model (quite time-consuming)
model <- cfboost(Surv(time, event) ~ bbsTime(time) + bbs(x.1) + bbs(x.2) + bbs(x.3),
                 control = ctrl, data = data, weights = weights)
}

##################################################
# estimate model with model choice decomposition #
##################################################

## the next model estimation is very time-consuming

\dontrun{
model_mod_choice <- cfboost(Surv(time, event) ~ bolsTime(time)
    + bbsTime(time, df = 1, center = TRUE)
    + bbsTime(time, df = 1, center = TRUE)
    + bols(x.1) + bbs(x.1, df = 1, center = TRUE)
    + bolsTime(time, z = x.1)
    + bbsTime(time, z = x.1, df = 1, center = TRUE)
    + bols(x.2) + bbs(x.2, df = 1, center = TRUE)
    + bolsTime(time, z = x.2)
    + bbsTime(time, z = x.2, df = 1, center = TRUE)
    + bols(x.3) + bbs(x.3, df = 1, center = TRUE)
    + bolsTime(time, z = x.3)
    + bbsTime(time, z = x.3, df = 1, center = TRUE),
  data = data, control = ctrl, weights = weights)
}

## proceed as above
}

\keyword{ package }
\keyword{ survival }
\keyword{ smooth }
\keyword{ models }
