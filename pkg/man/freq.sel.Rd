\name{freq.sel}
\alias{freq.sel}
\alias{print.fs}

\title{ Extracts the Frequency of Selection for Base-Learners }
\description{
  Functions to extract and print the selection frequencies of the
  base-learners from a \code{\link{cfboost}} object.
}
\usage{
freq.sel(object)
\method{print}{fs}(x, \dots)
}
\arguments{
  \item{object}{ an object of class \code{\link{cfboost}}. }
  \item{x}{ an object of class \code{fs} as returned by \code{freq.sel}.}
  \item{\dots}{ additional arguments (not implemented). }
}
\value{
  \code{freq.sel} returns a data.frame of class \code{fs}.
}
\author{ Benjamin Hofner }
\seealso{ \code{\link{cfboost}} for model fitting and \code{\link[CoxFlexBoost]{methods}} for other methods for \code{cfboost} objects. }
\examples{
\dontrun{
## model needs to be estimated in advance
freq.sel(cfb_model)
}
}

\keyword{ misc }
\keyword{ methods }
