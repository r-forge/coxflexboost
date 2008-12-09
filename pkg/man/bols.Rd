\name{bols}
\alias{bols}
\alias{bolsTime}

\title{ Linear Base-Learner for CoxFlexBoost }
\description{
  Defines a linear base-learner for boosting flexible, structured
  survival models. Both, time-constant base-learners and time-varying
  effects (and thus the baseline hazard) can be specified.
}
\usage{
bols(x, z = NULL, xname = NULL, zname = NULL, timedep = FALSE)
bolsTime(\dots)
}

\arguments{
  \item{x}{ factor or numeric. A vector containing data. }
  \item{z}{ factor or numeric. A vector containing data. }
  \item{xname}{ optional. Name of the variable given in \code{x}. }
  \item{zname}{ optional. Name of the variable given in \code{z}. }
  \item{timedep}{ logic. If \code{timedep = TRUE} the base-learner is
    taken to be time-dependent, i.e., \code{x} must be the time. For more
    details see below. }
  \item{\dots}{ further arguments passed to \code{bols}. }
}
\details{
  The function \code{bolsTime(\dots)} is a wrapper to
  \code{bols(\dots, timedep=TRUE)}.

  Time-dependent base-learners can be utilzed to specify the (log-)
  baseline hazard rate or time-varying effects. In the first case,
  the base-learner is specified as \code{bolsTime(time, z = NULL, \dots)}
  and in the later case, the covariate with a (potential) time-dependent
  effect is specified as \code{z}.
}
\value{
  An object of class \code{baselearner} (and class \code{bols}) is
  returned. It consits of design matrix of the base-learner. Further
  elements are returned as attributes (see \code{attr}) of the object.
  The attributes are mainly for internal use and are, e.g.,
  the current coefficient estimates (\code{coefs}) and a logical
  (\code{timedep}) indicating whether the base-learner specifies
  a time-varying effect or not.
}
\references{
  Benjamin Hofner, Torsten Hothorn and Thomas Kneib (2008),
  Variable Selection and Model Choice in Structured Survival Models.
  \emph{Department of Statistics, Technical Report No. 43}.
  \url{http://epub.ub.uni-muenchen.de/7901/}

  Thomas Kneib, Torsten Hothorn and Gerhard Tutz (2008),
  Variable selection and model choice in geoadditive regression models.
  \emph{Biometrics}. To appear.
  \url{http://epub.ub.uni-muenchen.de/2063/}
}
\seealso{ \code{\link{bbs}} for P-spline base-learners and
  \code{\link{cfboost}} for the boosting algorithm. }
\examples{
x1 <- rnorm(100)
time <- rexp(100,1)

## example for time-constant linear base-learner
str(bols(x1))

## example for linear time-varying effect for x1
str(bolsTime(time, z = x1))

## for more examples and effect decompositions see ?bbs
}

\keyword{ survival }
\keyword{ models }
