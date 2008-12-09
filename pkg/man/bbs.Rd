\name{bbs}
\alias{bbs}
\alias{bbsTime}

\title{ P-Spline Base-Learner for CoxFlexBoost }
\description{
  Defines a P-spline base-learner for boosting flexible, structured
  survival models. Both, time-constant base-learners and time-varying
  effects (and thus the baseline hazard) can be specified.
}
\usage{
bbs(x, z = NULL, knots = 20, degree = 3, differences = 2, df = 4,
    center = FALSE, xname = NULL, zname = NULL, timedep = FALSE)
bbsTime(\dots)
}

\arguments{
  \item{x}{ factor or numeric. A vector containing data. }
  \item{z}{ factor or numeric. A vector containing data. }
  \item{knots}{ a single numerical value defines the number of
    equidistant interior knots of the P-spline, a vector defines the
    knots itself, or \code{NULL}, in which case the number of interior
    knots is choosen to be \eqn{min(\code{length(x)},40)}. }
  \item{degree}{ degree of the P-spline. }
  \item{differences}{ order of the difference penalty used for the
    P-spline. Accepted values are 1,2 and 3. }
  \item{df}{ degrees of freedom. Used to specify the smoothness of the
    base-learner. Lower values for \code{df} correspond to smoother
    (i.e., \dQuote{weaker}) base-learners. For details see \cite{Hofner et
      al.}}
  \item{center}{ logic. If \code{center=TRUE}, the corresponding effect is
    re-parameterized such that the unpenalized part of the fit is substracted and
    only the deviation from the parametric, unpenalized part is fitted.
    The unpenalized part(s) have to be included as separate base-learners
    using \code{\link{bols}} or \code{\link{bolsTime}}. See \cite{Kneib
      et al.} for details. }
  \item{xname}{ optional. Name of the variable given in \code{x}. }
  \item{zname}{ optional. Name of the variable given in \code{z}. }
  \item{timedep}{ logic. If \code{timedep = TRUE} the base-learner is
    taken to be time-dependent, i.e., \code{x} must be a time. For more
    details see below. }
  \item{\dots}{ further arguments passed to \code{bbs}. }
}
\details{
  The function \code{bbsTime(\dots)} is a wrapper to \code{bbs(\dots, timedep=TRUE)}.

  Time-dependent base-learners can be utilzed to specify the (log-)
  baseline hazard rate or time-varying effects. In the first case,
  the base-learner is specified as \code{bbsTime(time, z = NULL, \dots)}
  and in the later case, the covariate with a (potential) time-dependent
  effect is specified as \code{z}.

  In most settings the number of \code{knots} should not exceed 20,
  to keep the computational burden in a reasonable dimension. Thus,
  a \code{\link{warning}} is issued if more knots are used.
}
\value{
  An object of class \code{baselearner} (and class \code{bbs}) is
  returned. It consits of the design matrix of the base-learner. Further
  elements are returned as attributes (see \code{\link{attr}}) of the object.

  The attributes are mainly for internal use. They consist, for example,
  of the degrees of freedom \code{df}. The corresponding value for
  \code{lambda} and the penalty (\code{pen}) are stored as further
  attributes. The current coefficient estimates are stored in
  the base-learner as \code{attr(\dots,"coefs")}.
}
\note{
  The degrees of freedom \code{df} are used to compute the corresponding
  smoothing parameter for the used P-spline. For more details see
  \cite{Hofner et al.}.
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
\seealso{ \code{\link{bols}} for linear base-learners and
  \code{\link{cfboost}} for the boosting algorithm. }
\examples{
x1 <- rnorm(100)
time <- rexp(100,1)
event <- sample(c(0,1), 100, replace=TRUE)
y <- Surv(time,event)

## example for time-constant base-learner
str(bbs(x1, knots = 20, df = 4))

## decomposition of a base-learner
decomp <- y ~ bols(x1) + bbs(x1,knots = 20, df = 1, center = TRUE)

## P-spline base-learner for the (log) baseline hazard
str(bbsTime(time)) # the same as bbs(time, timedep=TRUE)

## example for time-varying effect for x1
str(bbsTime(time, z = x1))

## decomposition of a time-varying effect
decomp <- y ~ bolsTime(time, z = x1) + bbsTime(time, z = x1, df = 1, center = TRUE)
}

\keyword{ survival }
\keyword{ smooth }
\keyword{ models }
