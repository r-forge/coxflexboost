##########################################
## Helper Functions Needed in cfboost() ##
##########################################

## convenience function for prediction
predict.baselearner <- function(object, newdata = NULL, newcoefs = NULL){
    ###
    # object        base-learner used for prediction
    # x             new x-values (default: NULL)
    # z             new z-values (default: NULL)
    if (!inherits(object,"baselearner"))
        stop(sQuote("object"), " must be of class ", sQuote("baselearner"))
    if (is.null(newcoefs)){
        return(as.vector(attr(object, "predict")(coefs = attr(object, "coefs"), newdata)))
    } else {
        if (length(newcoefs) != length(attr(object, "coefs")))
            stop(sQuote("newcoefs"), "must have length ", length(attr(object, "coefs")))
        return(as.vector(attr(object, "predict")(coefs = newcoefs, newdata)))
    }
}

## prediction for time-dependent variables on time_grid, using z-values from base-learner (if present)
predict_td <- function(object, ...)
    UseMethod("predict_td")

predict_td.baselearner <- function(object, time_grid){
    ###
    # object        base-learner for which to predict
    # time_grid     a list of time grids
    if (!inherits(object,"baselearner")) stop(sQuote("object"), " must be of class ", sQuote("baselearner"))

    x <- unlist(time_grid)
    xname <- get("xname", environment(attr(object, "predict")))
    z <- get("z", environment(attr(object, "predict")))
    if (!is.null(z)){
        zname <- get("zname", environment(attr(object, "predict")))
        z <- rep(z, each = length(time_grid[[1]]))
        newdata <- data.frame(cbind(x, z))
        names(newdata) <- c(xname, zname)
    } else {
        newdata <- data.frame(x)
        names(newdata) <- c(xname)
    }
    return(predict(object, newdata = newdata))
}

## convenience function for updating coefficients of base-learner 'object'
updatecoefs <- function(object, newcoefs){
    ###
    # object        base-learner that should be updated
    # newcoefs      new coefficients (added to the current coefficients)
    #
    # returned:     updated object
    #               to update a current base-learner you need to make a call like this:
    #               x <- updatecoefs(x, c(...))
    if (!inherits(object,"baselearner")) stop(sQuote("object"), " must be of class ", sQuote("baselearner"))
    attr(object, "coefs") <- attr(object, "coefs") + newcoefs
    return(object)
}

## convenience function to extract coefficients from base-learner (and multiply them with nu)
getcoefs <-  function(object, nu = 1){
    ###
    # object        base-learner which holds the coefficients you want to extract
    # nu            (optional) Specify nu to get the reduced coefficients from boosting algorithm
    attr(object, "coefs") * nu
}

## function for calculation of the initial fit (i.e. offset)
getoffset <- function(y, which.offset = "mle"){
#    if (which.offset == "mle"){
#        logLH <- function(f, y){
#            ## risk for constant fit "f"
#            sum((y[,2] * f - exp(f) * y[,1]))
#        }
#        if (plot) {
#            f <- seq(interval[1],interval[2], length = 300)
#            llh <- rep(NA, length = length(f))
#            for (i in 1:length(f)) llh[i] <- logLH(f[i],y)
#            plot(f, llh, type = "l")
#        }
#        os <- optimize(logLH, interval = interval, y = y, maximum = TRUE)$maximum
#        return(os)
#    }
    if (which.offset == "mle"){
        return(log(sum(y[,2]) / sum(y[,1])))
    }
    if (which.offset == "zero"){
        return(0)
    }
}


## controls for boosting algorithm
## (adapted version from mboost)
boost_control <- function(mstop = 100, nu = 0.1, maxit = 30000, risk = c("inbag", "oobag", "none"),
                          which.offset = c("mle", "zero"), savedata = TRUE,
                          center = FALSE, trace = TRUE, hardStop = TRUE) {

    which.offset <- match.arg(which.offset)
    risk <- match.arg(risk)
    RET <- list(mstop = mstop, nu = nu, maxit = maxit,
                risk = risk, which.offset = which.offset,
                savedata = savedata, center = center, trace = trace,
                hardStop = hardStop)
    class(RET) <- c("boost_control")
    RET
}

## data preprocessing
## (adapted version from mboost)
boost_dpp <- function(formula, data, weights = NULL, na.action = na.omit, ...) {
    if (is.null(weights)) {
        env <- ModelEnvFormula(formula, data, na.action = na.action, ...)

        y <- env@get("response")
        if (length(y) != 1)
            stop("cannot deal with multivariate response variables")
        y <- y[[1]]
        x <- env@get("designMatrix")

        RET <- list(x = x, y = y)
        RET$formula <- formula
        RET$menv <- env
        class(RET) <- "boost_data"
        return(RET)
    } else {
        if (NROW(data) != length(weights))
            stop(sQuote("weights"), " is not of length ", NROW(data))
        if (!all(weights == 1 | weights == 0))
            stop(sQuote("weights"), "must be all 1 (learning data) or 0 (test data)")

        ##
        # make two environments (in-bag, out-of-bag)
        data_inbag <- data[as.logical(weights),]
        data_oob <- data[as.logical(1 - weights),]

        env_inbag <- ModelEnvFormula(formula, data_inbag, na.action = na.action,...)
        env_oob <- ModelEnvFormula(formula, data_oob, na.action = na.action,...)

        ## build up appropriate objects (inbag)
        y_inbag <- env_inbag@get("response")
        if (length(y_inbag) != 1)
            stop("cannot deal with multivariate response variables")
        y_inbag <- y_inbag[[1]]
        x_inbag<- env_inbag@get("designMatrix")

        RET_inbag <- list(x = x_inbag, y = y_inbag)
        RET_inbag$formula <- formula
        RET_inbag$menv <- env_inbag
        class(RET_inbag) <- "boost_data"

        ## same for oobag
        y_oob <- env_oob@get("response")
        if (length(y_oob) != 1)
            stop("cannot deal with multivariate response variables")
        y_oob <- y_oob[[1]]
        x_oob<- env_oob@get("designMatrix")

        RET_oob <- list(x = x_oob, y = y_oob)
        RET_oob$formula <- formula
        RET_oob$menv <- env_oob

        RET_inbag$oob <- RET_oob

        return(RET_inbag)
    }
}
