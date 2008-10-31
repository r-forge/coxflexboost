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
     type = "b", ylab = expression(f[partial]), add_rug = TRUE, ...)
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
  \item{i}{ integer. Index specifying the model to extract. See example
       for more details.}
  \item{\dots}{ additional arguments (not used a.t.m.) }
}
\details{
  The function \code{print} prints basic information about the model
  specified by \code{object} and returns the argument \code{object}
  invisible. \code{summary} is just a wrapper to the \code{print}
  function.

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
\author{ Benjamin Hofner }
\seealso{ \code{\link{cfboost}} for model fitting. }
\examples{
## add more examples here
\dontrun{
## estimate model first
(stop.opt <- mstop(cfb_model))
cfb_model <- cfb_model[op.opt]
}
}

\keyword{ misc }
\keyword{ methods }
