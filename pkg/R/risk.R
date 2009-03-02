##############################################
## Methods to Extract and Plot Minimal Risk ##
##############################################

### methods: extract risk for cfboost objects
risk <- function(object, ...)
    UseMethod("risk")

risk.cfboost <- function(object, ...){
    return(object$risk)
}


## methods: plot methods for risk objects
plot.oobag <- function(x, xlab = "iteration", ylab = "risk (= - log likelihood) in validation sample", type = "l",
                       mstop = TRUE, xlim = NULL, ...){
    if (!is.null(xlim)){
        if (length(xlim) != 2)
            stop(sQuote("xlim"), " must be a vector of length 2")
        if (xlim[1] >= xlim[2])
            stop("xlim[1] >= xlim[2]")
        if (!all(xlim %in% 1:length(x)))
            stop (sQuote("xlim"), " must be (an integer valued) element of 1:mstop")
        index <- seq(from = xlim[1], to = xlim[2])
    } else {
        index <- 1:length(x)
    }

    if (mstop)
        xlab <- paste(xlab, " (mstop = ", print(x, print = FALSE)$iteration , ")", sep="")
    plot.default(index, x[index], xlab = xlab, ylab = ylab, type = type, ...)
    if (mstop)
        abline(v = print(x, print = FALSE)$iteration , lty = "dashed")
}

## wrapper to method for oobag-risk
plot.inbag <- function(x, xlab = "iteration", ylab = "risk (= - log likelihood) in learning sample", type = "l",
                       mstop = TRUE, xlim = NULL, ...){
    plot.oobag(x, xlab, ylab, type, mstop, xlim, ...)
}


### methods: extract and print minimal risk

## calculate out of bag minimal risk
print.oobag <- function(x, print = TRUE, ...){
    minrisk <- min(x)
    iteration <- which.min(x)
    RET <- list(minrisk = minrisk, iteration = iteration)
    class(RET) <- "oobag"
    if (print)
         cat("minimal risk (validation sample)", minrisk ,"in iteration: ", iteration, "\n")
    invisible(RET)
}

## calculate in-bag minimal risk
print.inbag <- function(x, print = TRUE, ...){
    minrisk <- min(x)
    iteration <- which.min(x)
    RET <- list(minrisk = minrisk, iteration = iteration)
    class(RET) <- "oobag"
    if (print)
        cat("minimal risk (learning sample)", minrisk ,"in iteration: ", iteration, "\n")
    invisible(RET)
}
