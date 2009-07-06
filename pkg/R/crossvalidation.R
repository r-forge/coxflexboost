cv <- function(object, ...)
    UseMethod("cv")

cv.cfboost <- function(object, folds, grid = c(1:mstop(object, opt=FALSE)), ...){

    oobrisk <- matrix(0, nrow = ncol(folds), ncol = length(grid))
    ctrl <- object$control
    ctrl$risk <- "oobag"
    # fehlt noch: ctrl$savedata <- FALSE
    # fehlt noch: ctrl$saveensss <- FALSE

    if (is.null(object$data))
        stop(sQuote("object"), " does not contain data. Estimate model with option ", sQuote("savedata = TRUE"))

    call <- deparse(object$call)
    data <- object$data$data
    formula <- object$data$formula

    myapply <- lapply
    if (ctrl$parallel && require("multicore")) {
        if (!multicore:::isChild()) {
            myapply <- mclapply
            if (ctrl$trace) {
                ctrl$trace <- FALSE
                cat("\n Running in parallel with `trace = FALSE'\n")
            }
        }
    }

    ## free memory
    rm("object")

    i <- 0

    dummyfct <- function(weights, control, data, formula, grid){
        i <<- i + 1
        if (ctrl$trace) cat("\n>>> Fold ", i, "\n\n")
        model <- cfboost(formula, data = data, control = control, weights = weights)
        ret <- risk(model)[grid]
        rm("model")
        ret
    }

    oobrisk <- myapply(1:ncol(folds), function(i)
                       dummyfct(folds[,i], control = ctrl, data = data, formula = formula, grid = grid),
                       ...)
    oobrisk <- t(as.data.frame(oobrisk))
    oobrisk <- oobrisk/colSums(folds == 0)
    colnames(oobrisk) <- grid
    rownames(oobrisk) <- 1:nrow(oobrisk)
    attr(oobrisk, "call") <- call
    attr(oobrisk, "mstop") <- grid
    attr(oobrisk, "risk") <- "empirical risk (neg. log likelihood)"
    class(oobrisk) <- "cv"
    oobrisk
}


print.cv <- function(x, ...) {
    cat("\n\t Cross-validated risk \n\t Call:",
              attr(x, "call"), "\n\n")
    print(colMeans(x))
    cat("\n\t Optimal number of boosting iterations:", mstop(x), "\n")
    return(invisible(x))
}

plot.cv <- function(x, ylab = attr(x, "risk"), ylim = range(x),
                        main = attr(x, "call"), ...) {

    cm <- colMeans(x)
    plot(1:ncol(x), cm, ylab = ylab, ylim = ylim,
         type = "n", lwd = 2,
         xlab = "Number of boosting iterations",
         main = main, ...)
    out <- apply(x, 1, function(y) lines(1:ncol(x),y, col = "lightgrey"))
    rm(out)
    ms <- which.min(cm)
    lines(c(ms, ms), c(min(c(0, ylim[1] * ifelse(ylim[1] < 0, 2, 0.5))), cm[ms]), lty = 2)
    lines(1:ncol(x), cm, type = "l")
}
