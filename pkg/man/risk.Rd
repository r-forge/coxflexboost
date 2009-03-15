\name{risk}
\alias{risk}
\alias{risk.cfboost}
\alias{plot.oobag}
\alias{plot.inbag}
\alias{print.oobag}
\alias{print.inbag}

\title{ Methods to Extract and Plot the (Minimal) Risk }
\description{
  Methods to extract the risk from models fitted via \code{cfboost} and
  to print and show the results.
}
\usage{
\method{risk}{cfboost}(object, ...)
\method{plot}{oobag}(x, xlab = "iteration",
     ylab = "risk (= - log likelihood) in validation sample",
     type = "l", mstop = TRUE, xlim = NULL, ...)
\method{plot}{inbag}(x, xlab = "iteration",
     ylab = "risk (= - log likelihood) in learning sample",
     type = "l", mstop = TRUE, xlim = NULL, ...)
\method{print}{oobag}(x, print = TRUE, ...)
\method{print}{inbag}(x, print = TRUE, ...)
}

\arguments{
  \item{object}{ an object of class \code{cfboost}. }
  \item{x}{ an object of either class \code{oobag} or \code{inbag}. }
  \item{xlab}{  A title for the x axis. }
  \item{ylab}{  A title for the y axis. }
  \item{type}{ type of plot to be drawn. See \code{\link{plot}} for more
       details. }
  \item{mstop}{logic. If \code{mstop = TRUE} the optimal stopping
    iteration is given in the x axis and a vertical line is added to the plot. }
  \item{xlim}{ the limits of the x axis, see \code{\link{plot.default}}.
    Note: reversed axis are not allowed.}
  \item{print}{ logic. If \code{print = TRUE} the result is written to
    the screen. Otherwise the results are only returned invisible. }
  \item{\dots}{ further arguments to be pased on. }

}
\details{
  The function \code{risk} allows to extract the risk (vector) from an
  \code{object} (of class \code{cfboost}).

  Depending on the type of \code{risk} (\code{oobag} or \code{inbag})
  different plotting facilities exist. Note however, that
  \code{plot.inbag} is just a wrapper for \code{plot.oobag} with
  different defaults.

  In the same manner a \code{print} method is defined for
  \code{oobag} and \code{inbag} risk.
}
\value{
  \code{risk} returns the empirical inbag or out-of-bag risk from a
  \code{cfboost} object.

  \code{print} returns the given \code{object} invvisble.
}
\seealso{ for other methods for \code{cfboost} objects see
  \code{\link[CoxFlexBoost]{methods}}.
  See \code{\link{cv}} for cross-validated risk.}
\examples{
## see for example ?CoxFlexBoost-package
## for the usage of risk() and plot(risk()) etc.

}

\keyword{ methods }
\keyword{ misc }