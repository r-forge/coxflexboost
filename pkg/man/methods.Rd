\name{methods}
\alias{print.cfboost}
\alias{summary.cfboost}
\alias{plot.cfboost}
\alias{coef.cfboost}
\alias{Extract}
\alias{Subscript}
\alias{[}
\alias{[.cfboost}

\title{ (Generic) Functions for cfboost Objects }
\description{
  Functions to print or plot objects from \code{\link{cfboost}}. Further
  functions to extract the coefficients are specified. The extract
  function allows to reduce the model to an earlier boosting iteration,
  as e.g., the optimal stopping iteration computed with \code{\link{mstop}}.
}
\usage{
\method{print}{cfboost}(x, ...)
\method{summary}{cfboost}(object, ...)
\method{plot}{cfboost}(x, which = NULL, ask = TRUE && dev.interactive(),
     type = "b", ylab = expression(f[partial]), add_rug = TRUE,
     color.palette = c("heat.colors", "terrain.colors", "topo.colors",
     "cm.colors", "rainbow", "none"), ...)
\method{coef}{cfboost}(object, ...)
\method{[}{cfboost}(object, i, ...)
}
\arguments{
  \item{x}{ an object of class \code{cfboost} }
  \item{object}{ an object of class \code{cfboost} }
  \item{which}{ integer vector. Only the plots for the given
       base-learners are plotted. Per default, only the selected
       base-learners are plotted. }
  \item{ask}{ logic. If true a user input is required between each plot,
       see \code{\link[=par]{par(ask=.)}}. }
  \item{type}{ type of plot to be drawn. See \code{\link{plot}} for more
       details.}
  \item{ylab}{ A title for the y axis. }
  \item{add_rug}{ logic. Determines if \code{\link{rug}}s are added. }
  \item{color.palette}{ character. Determines how time-varying effects
       of non-binary covariates should be ploted. A color palette of
       "heat.colors", "terrain.colors", "topo.colors", "cm.colors",
       "rainbow" is specified (as character). See \code{\link{rainbow}}
       for details. If "none" is specified, all observations are
       printed in black. }
  \item{i}{ integer. Index specifying the model to extract. See example
       for more details.}
  \item{\dots}{ additional arguments (not used a.t.m.) }
}
\details{
  The function \code{print} prints basic information about the model
  specified by \code{object} and returns the argument \code{object}
  invisible. \code{summary} is a wrapper to the \code{print}
  function put gives additional information on the frequency of
  selection of the base-learners (see \code{\link{freq.sel}}).

  \code{plot} gives a simple plotting interface for the estimated
  effects of the model \code{object}.

  The function \code{coef} extracts the coefficients from the model
  \code{object}.

  \code{object[i]} extracts the model from the i-th boosting iteration.
}
\value{
  \code{object[i]} returns the \code{cfboost} model-object from the i-th
  iteration.
}
\seealso{ \code{\link{cfboost}} for model fitting.
  \code{\link{freq.sel}} for the selection frequencies of base-learners.}
\examples{
## a simple example

#################
# simulate data #
#################

set.seed(4321)

## sample covariates first
X <- matrix(NA, nrow=400, ncol=3)
X[,1] <- runif(400, -1, 1)
X[,2] <- runif(400, -1, 1)
X[,3] <- runif(400, -1, 1)

## time-constant hazard rate
lambda <- function(time, x){
  exp(0.1 * time - 0.7 * x[1] + 0.5 * x[2])
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

ctrl <- boost_control( mstop = 80, risk="oobag")
weights <- c(rep(1, 300), rep(0, 100))

## non-censored observations (in-bag)
sum(data$event[weights==1])/nrow(data[weights==1,])

## fit (a simple) model
model <- cfboost(Surv(time, event) ~ bolsTime(time) + bols(x.1) + bols(x.2) + bols(x.3),
                 control = ctrl, data = data, weights = weights)

#########################
# processing the output #
#########################

## estimate optimal mstop
(stop.opt <- mstop(model))
model_opt <- model[stop.opt]

model
## summary for mstop = 100
summary(model)
## summary (with "optimal" mstop)
summary(model_opt)

## plot of baseline hazard
plot(model_opt, which=1)
## plot of x.1 und x.2
par(mfrow=c(1,2))
plot(model_opt, which = c(2,3), ask = FALSE)

## extract coefficients
coef(model_opt)
## almost correct estimations
}

\keyword{ misc }
\keyword{ methods }
