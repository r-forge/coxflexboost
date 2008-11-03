\name{mstop}
\alias{mstop}
\alias{mstop.cfboost}

\title{ Function to Extract (Optimal) Stopping Iteration }
\description{
  Extracts the number of boosting iterations from a model that where
  performed or estimates the optimal number of boosting iterations.
}
\usage{
mstop(object, ...)
\method{mstop}{cfboost}(object, opt = TRUE, ...)
}

\arguments{
  \item{object}{ object of class \code{cfboost}. }
  \item{opt}{ logic. If \code{opt = FALSE} the pre-defined stopping
    iteration is returned, if \code{opt = TRUE} the optimal stopping
    iteration (in the given range of iterations) is returned. }
  \item{\dots}{ (not used a.t.m.) }
}
\details{
  The calculated risk is negative maximum log likelihood for each
  boosting step. The function \code{mstop} can be used to extract the
  number of boosting iterations that were performed (\code{opt = FALSE})
  or the optimal number of boosting iterations \code{opt = TRUE} with
  respect to the in-bag risk or out-of-bag risk. In the first case, the
  risk is computed on the learning sample, in the latter case it is
  computed on the validation sample. The samples can be specified by
  \code{weights} that are given to \code{\link{cfboost}}. The risk that
  is to be computed is specified with the call to \code{\link{cfboost}},
  to be more precise by the \code{\link{boost_control}} function. If an
  out-of-bag sample is specified, the risk that is to be computed can be
  set to both alternatives, if not, only the inbag risk is appropriate.
}
\value{
  The (optimal) number of boosting iterations is returned.
}
\author{ Benjamin Hofner }
\seealso{ \code{\link{boost_control}} for the specification of
  risk-type and \code{\link{cfboost}} for the specification of
  validation samples. }
\examples{
## see for example ?cfboost for usage of mstop()
}

\keyword{ methods }
\keyword{ misc }
