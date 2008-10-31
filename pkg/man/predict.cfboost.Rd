\name{predict.cfboost}
\alias{fitted.cfboost}
\alias{predict.cfboost}

\title{ Predictions (for New Observations) }
\description{
  Functions to extract the fitted values of models fitted via
  \code{cfboost} and make predictions for new observations.
}
\usage{
\method{fitted}{cfboost}(object, type = c("hazard", "log-hazard"), ...)
\method{predict}{cfboost}(object, newdata = NULL, type = c("hazard", "log-hazard"), trace = TRUE, ...)
}

\arguments{
  \item{object}{ object of class \code{cfboost} }
  \item{newdata}{ optional. a data frame in which to look for variables with
          which to predict. }
  \item{type}{ character. Specifies the type of the returned prediction
    (either the hazard rate or the log hazard rate). }
  \item{trace}{ ??? What does this option do ??? }
  \item{\dots}{ Further arguments to be passed to subsequent functions. }
}
\details{
  Note that the data frame newdata needs to have the same variable names
  as the original data (or as specified in the base-learners via the
  \code{xname} and \code{zname} options).
}
\value{
  \code{fitted} returns a matrix of the fitted values in the last
  boosting iteration of \code{object}.

  \code{predict} can be used to make predictions on new data sets and
  returns a vector of the predicted (log) hazard rates.
}
\author{ Benjamin Hofner }
\seealso{ \code{\link{cfboost}} for the boosting procedure and
  \code{\link{bbs}} or \code{\link{bols}} for the specification of the base-learners.}
\examples{
## a simple example
set.seed(1234)

## sample covariates first
X <- matrix(NA, nrow=400, ncol=3)
X[,1] <- runif(400, -1, 1)
X[,2] <- runif(400, -1, 1)
X[,3] <- runif(400, -1, 1)

## time-dependent hazard rate
lambda <- function(time, x){
  exp(log(time) + 0.7 * x[1] + x[2]^2)
}

## specify censoring function
cens_fct <- function(time, mean_cens){
  censor_time <- rexp(n = length(time), rate = 1/mean_cens)
  event <- (time <= censor_time)
  t_obs <- apply(cbind(time, censor_time), 1, min)
  return(cbind(t_obs, event))
}
data <- rSurvTime(lambda, X, cens_fct, mean_cens = 5)

ctrl <- boost_control(risk="oobag")
weights <- c(rep(1, 300), rep(0, 100))

## fit (a simple) model
model <- cfboost(Surv(time, event) ~ bbs(x.1) + bbs(x.2) + bbs(x.3),
                 control = ctrl, data = data, weights = weights)

## extract fitted values
fitted(model)

## make prediction for the first out-of-bag observation
predict(model, newdata = data[301,])

}
\keyword{ models }
\keyword{ methods }
\keyword{ survival }

