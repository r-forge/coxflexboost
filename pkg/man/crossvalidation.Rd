\name{crossvalidation}
\alias{cv}
\alias{cv.cfboost}
\alias{print.cv}
\alias{plot.cv}

\title{ Cross-Validation}
\description{
  Function for the estimation of \code{\link{mstop}} via
  cross-validation and (generic) functions to print or plot the resuls.
}
\usage{

\method{cv}{cfboost}(object, folds,
   grid = c(1:mstop(object, opt=FALSE)), ...)
\method{print}{cv}(x, ...)
\method{plot}{cv}(x, ylab = attr(x, "risk"), ylim = range(x),
     main = attr(x, "call"), ...)
}
\arguments{
  \item{object}{ an object of class \code{\link{cfboost}}.}
  \item{folds}{ a weight matrix with number of rows equal to the number
                of observations. The number of columns corresponds to
                the number of cross-validation runs.}
  \item{grid}{ a vector of iterations the empirical risk
               is to be evaluated for. Per default the empirical
               risks for all iterations \code{1:mstop} are computed and returned.}
  \item{x}{ an object of class \code{cv} }
  \item{ylab}{ A title for the y axis. }
  \item{ylim}{ the y limits of the plot. }
  \item{main}{ the main title of the plot. }
  \item{\dots}{ additional arguments to be passed to \code{\link{plot}}. }
}
\details{
  The number of boosting iterations is a hyper-parameter of the
  boosting algorithms. Cross-validated estimates of the empirical risk
  for different values of \code{mstop} (as given by
  \code{\link{grid}}) are computed, which are used to choose the
  appropriate number of boosting iterations to be applied.

  Different forms of cross-validation can be applied, for example,
  5-fold or 10-fold cross-validation. Bootstrapping is not implemented
  so far. The \code{weights} are defined via the \code{folds} matrix.
  A.t.m. they can only be used to specify a learning
  sample which consists of observations with \code{weights == 1} and
  and an out-of-bag sample with \code{weights == 0}. The latter
  is used to determine the empirical risk (negative log likelihood).
}
\value{
  \code{cv} returns an object of class \code{cv}, which consists of
  a matrix of empirical risks and some further attributes.
}
\seealso{ \code{\link{cfboost}} for model fitting.
   See \code{\link{risk}} for methods to extract the inbag and
   out-of-bag risk and \code{\link{mstop}} for functions to
   extract the (optimal) stopping iteration (based on cross-validation
   or on the inbag and out-of-bag risk).}
\examples{
## fit a model with all observations first

\dontrun{
## (as this takes some minutes)
set.seed(1234)
## sample covariates first
X <- matrix(NA, nrow=400, ncol=3)
X[,1] <- runif(400, -1, 1)
X[,2] <- runif(400, -1, 1)
X[,3] <- runif(400, -1, 1)

## time-dependent hazard rate
lambda <- function(time, x){
   exp(0 * time + 0.7 * x[1] + x[2]^2)
}

## specify censoring function
cens_fct <- function(time, mean_cens){
  censor_time <- rexp(n = length(time), rate = 1/mean_cens)
  event <- (time <= censor_time)
  t_obs <- apply(cbind(time, censor_time), 1, min)
  return(cbind(t_obs, event))
}
daten <- rSurvTime(lambda, X, cens_fct, mean_cens = 5)

ctrl <- boost_control( mstop = 100, risk="none")
## fit (a simple) model
model <- cfboost(Surv(time, event) ~ bbs(x.1) + bbs(x.2) + bbs(x.3),
                 control = ctrl, data = daten)
}

## 5 -fold cross-validation

\dontrun{
## (as this takes some minutes)
n <- nrow(daten)
k <- 5
ntest <- floor(n / k)
cv5f <- matrix(c(rep(c(rep(0, ntest), rep(1, n)), k - 1),
                       rep(0, n * k - (k - 1) * (n + ntest))), nrow = n)
cvm <- cv(model, folds = cv5f)
print(cvm)
plot(cvm)
mstop(cvm)
}
}

\keyword{ misc }
\keyword{ methods }
