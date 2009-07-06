\name{boost_control}
\alias{boost_control}
\title{ Specify Controls for the Boosting Algorithm }
\description{
  Specify control and tuning parameters of the boosting algorithm.
}
\usage{
boost_control(mstop = 100, nu = 0.1, maxit = 30000, risk = c("inbag","oobag", "none"),
    which.offset = c("mle", "zero"), savedata = TRUE, trace = TRUE,  parallel = require("multicore"), hardStop = TRUE)
}
\arguments{
  \item{mstop}{ integer. Initial number of boosting iterations }
  \item{nu}{ double. Step-size between 0 and 1 }
  \item{maxit}{ integer. Controls the maximum number of iterations that
    are used for numerical optimazation via \code{\link[=optim]{optim(\dots, method = "BFGS")}}. }
  \item{risk}{ character. Determines how the empirical risk should be
    computed. It can take the values \code{risk = "inbag"}, i.e., the risk is
    computed for the learning sample, \code{risk = "oobag"}, i.e., the risk is
    computed for the validation sample and \code{risk = "none"}.}
  \item{which.offset}{ character. Indictating the choice of offset. This
    can be either the maximum likelihood estimator of the null model
    (\code{which.offset = "mle"}) or zero (\code{which.offset = "zero"}).}
  \item{savedata}{ logic. Indicating whether the data should be saved in
    the returned \code{cfboost} object. }
  \item{trace}{ logic. Determining if status information should be
    printed during estimation of the model. } 
\item{parallel}{logic. enable parallelization using \pkg{multicore} technology
    for crossvalidation (in function \code{\link{cv}}). If \pkg{multicore} is available
    parallelization is used per default but can be switched off. }
  \item{hardStop}{ logic. Indicates if the initial value of
    \code{mstop} is the maximal number of iterations (\code{mstop =
      TRUE}) or if it should be increased if the algorithm did not
    converge until \code{mstop}. The latter is only possible with
    \code{risk = "oobag"}. }
}
\value{
  An object of class \code{boost_control} is returned (as a list).
}
\seealso{ \code{\link{cfboost}} }

\keyword{ misc }

