\name{cfboost}
\alias{cfboost}
\alias{cfboost.formula}
\alias{cfboost_fit}
\title{ Likelihood Based Boosting for Flexible, Structured Survival
  Models }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cfboost(x, \dots)
\method{cfboost}{formula}(formula, data = list(), weights = NULL,
   na.action = na.omit, control = boost_control(), \dots)
}

\arguments{
  \item{x}{ an object to be pased to cfboost.formula }

  \item{formula}{ formula of the model to be fit }
  \item{data}{ data frame. Contains the data of the model that is to be fitted }
  \item{weights}{ integer vector. Optional weights of the observations in
    \code{data}. Currently only weights 0 or 1 are supported. See Details. }
  \item{na.action}{ a function that defines how to handle \code{NA}s in
    the data }
  \item{control}{ an object of class \code{\link{boost_control}} }
  \item{\dots}{ further arguments to be passed to subfunctions }
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
  \item{data}{data object after data pre-processing}
  \item{ensemble}{vector of selected base-learners}
  \item{ensembless}{list of coefficient estimates in each iteration}
  \item{fit}{vector of fitted values}
  \item{offset}{offset value}
  \item{control}{control parameters as specified}
  \item{risk}{the negative maximum log likelihood}
  \item{weights}{weights used for fitting}
  \item{df}{estimated degrees of freedom for smooth base-learners}
  \item{coefs}{estimated coefficients (in the final boosting iteration)}
  \item{predict}{function for prediction (to be accessed via
        \code{\link[=predict.cfboost]{predict}})}
}
\references{
  ~~ PUT TECH.REP. HERE ~~
}
\author{ Benjamin Hofner }
\seealso{ \code{\link{boost_control}} for the control arguments for
  \code{cfboost} and \code{\link{bols}} or \code{\link{bbs}} for the
  available base-learners.
  Methods to process the output can be found in \code{\link[CoxFlexBoost]{methods}}
}
\examples{
## put an example here
}

\keyword{ survival }
\keyword{ smooth }
\keyword{ models }

