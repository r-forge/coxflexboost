##########################################
## Helper Functions Needed in cfboost() ##
##########################################

## convenience function for prediction
predict.baselearner <- function(object, newdata = NULL, newcoefs = NULL, ...){
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
#predict_td <- function(object, ...)
#    UseMethod("predict_td")

#predict_td.baselearner <- function(object, time_grid){
predict_td <- function(object, time_grid){
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
                          trace = TRUE, hardStop = TRUE) {

    which.offset <- match.arg(which.offset)
    risk <- match.arg(risk)
    RET <- list(mstop = mstop, nu = nu, maxit = maxit,
                risk = risk, which.offset = which.offset,
                savedata = savedata, trace = trace, hardStop = hardStop)
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

## helper function for df2lambda as used in bbs()
helper_fct <- function(y, x, offset, pen, subdivisions = 100){
    time <- y[,1]
    delta <- y[,2]
    coefs <- rep(0, ncol(x))

    if(attr(x, "timedep")){
        ## make grid and compute grid-width if there is ANY time-dependent base-learner
        sub = subdivisions      ## subdivisions of time
        n = length(time)        ## number of observations
        grid <- function(upper, length){
            ## helper function to compute grid
            seq(from = 0, to = upper, length = length)
        }
        ## make grid
        grid <- lapply(time, grid, length = sub)
        trapezoid_width <- rep(NA, n)
        for (i in 1:n)
            trapezoid_width[i] <- grid[[i]][2]  # = second element, as first element == 0 and the grid equidistant for every i
    } else {
        grid = NULL
        trapezoid_width = NULL
    }
    ## build design matrix for base-learner
    if (attr(x, "timedep")){
        xd <- unlist(grid)
        xname <- get("xname", environment(attr(x, "predict")))
        zd <- get("z", environment(attr(x, "predict")))
        if (!is.null(zd)){
            zname <- get("zname", environment(attr(x, "predict")))
            zd <- rep(zd, each = length(grid[[1]]))
            newdata <- data.frame(cbind(xd, zd))
            names(newdata) <- c(xname, zname)
        } else {
            newdata <- data.frame(xd)
            names(newdata) <- c(xname)
        }
        desMat <- attr(x,"designMat")(newdata = newdata)
    } else {
        desMat <- attr(x,"designMat")()
    }
    F_pen <- function(lambda, coefficients=NULL){
        if(!is.null(coefficients))
            coefs <- coefficients
        Fisher_mat <- integr_fisher(x, coefs, desMat,
                                    predictions = list(offset = exp(offset), tconst = 1, td = 1),
                                    controls = list(grid = grid, trapezoid_width = trapezoid_width, upper = time, nu = 0.1, which = "fisher"))
        if (is.null(pen)) pen <- 0
        return(list(F = Fisher_mat, F_pen = Fisher_mat + lambda * pen))
    }
    return(F_pen)
}

penalty <- function(obj){
    dummy <- attr(obj, "pen")
    if (!is.null(dummy))
        dummy <- dummy()
    return(dummy)
}

## (adapted version from mboost)
complete_cases <- function(x, y = NULL, z = NULL) {

    tmp <- list(x = x, y = y, z = z)
    tmp <- tmp[!sapply(tmp, is.null)]
    rowSums(sapply(tmp, is.na)) == 0
}
