########################
## Boosting Algorithm ##
########################

### generic method for likelihood-based boosting with component-wise P-splines
### for fitting Cox-type additive models (flexible Cox boosting)
cfboost <- function(x, ...)
    UseMethod("cfboost")

### formula interface
cfboost.formula <- function(formula, data = list(), weights = NULL, na.action = na.omit,  control = boost_control(), ...) {
    ## construct design matrix etc.
    object <- boost_dpp(formula, data, weights, na.action)
    ## fit the ensemble
    object$input <- object$menv@get("input")
    if (!is.null(weights))
        object$oob$input <- object$oob$menv@get("input")
    RET <- cfboost_fit(object, control = control, data = data, weights = weights, ...)
    RET$call <- match.call()
    return(RET)
}

### fitting-function
cfboost_fit <- function(object, control = boost_control(), data, weights = NULL, ...) {

    ## data and base-learner
    x <- object$input
    class(x) <- "list"

    oob <- list(x = object$oob$input , y = object$oob$y)
    if (!is.null(oob$x)) class(oob$x) <- "list"

    y <- object$y
    if (!inherits(y, "Surv")) stop("response is not an object of class ", sQuote("Surv"))

    ## hyper parameters
    mstop <- control$mstop
    risk <- control$risk
    nu <- control$nu
    trace <- control$trace
    tracestep <- options("width")$width / 2
    maxit <- control$maxit # maximum iterations in optim (see PMLE())
    which.offset <- control$which.offset
    hardStop <- control$hardStop # option: "continue boosting if minimum not reached" or "stop"

    ## check if (enough) oob-data is present if risk="oobag"
    if (risk == "oobag"){
        checksum <- sum(weights==0)
        if(checksum == 0)
            stop("All observations are used for estimation.\nSpecify some weights equal to zero for ", sQuote("risk = oobag"))
        if (checksum < length(weights)/10)
            warning("Less than 1/10 of the data used for out-of-bag risk.\n", sQuote("object$risk"), " might be misleading.")
    }

    ## the ensemble
    ens <- rep(NA, mstop)
    ensss <- vector(mode = "list", length = mstop)

    ## vector of empirical risks for all boosting iterations
    mrisk <- numeric(mstop)
    mrisk[1:mstop] <- NA

    maxll <- numeric(length(x))
    coefs <- logLH <- vector(mode = "list", length = length(x))

    fit <- fit_oob <- offset <- getoffset(y, which.offset)
    if (trace)
        cat("Offset: ", offset, "\n")


    mstart <- 1
    hSi <- 1     # number of iterations in the repeat loop
    df_est <- matrix(NA, nrow = mstop, ncol = length(x)) # matrix of estimated degrees of freedom

    ## compute df2lambda which depends on the offset and on y
    for (i in 1:length(x)){
        if (!is.null( attr(x[[i]], "lambda"))){
            attr(x[[i]],"lambda") <-  attr(x[[i]], "lambda")(y, offset)
            attr(x[[i]],"pen") <- attr(x[[i]],"pen")(attr(x[[i]],"lambda"))
        }
    }

    ##################################
    #### start boosting iteration ####
    ##################################
    repeat{
      for (m in mstart:mstop) {
        if (trace)
          cat("Step ", m, "\n")

        ## fit MLE component-wise
        for (i in 1:length(x)) {
            maxll[i] <- NA
            dummy_ens <- ens[1:m]   # get the first m-1 selected base-learners
            dummy_ens[m] <- i       # and set the m-th base-learner temporarily to i
            ## try to compute the (component-wise) penalized MLE
            dummy <- try(PMLE(y, x, offset, fit, dummy_ens, nu, maxit))
            if (inherits(dummy, "try-error")) next
            coefs[[i]] <- dummy$par
            maxll[i] <- dummy$maxll
            logLH[[i]] <- dummy$logLH
            if (!is.null(dummy$df)) df_est[m,i] <- dummy$df
        }
        if (all(is.na(maxll)))
            stop("could not fit base learner in boosting iteration ", m)

        ## select base-learner
        xselect <- which.max(maxll)

        ## output for debugging
        if (trace)
            cat("\tSelected: ", xselect, "   with log LH ", maxll[xselect],"\n", sep = "")
        ## update step
        fit <- fit + nu * x[[xselect]] %*% coefs[[xselect]]

        ## save the model, i.e., the selected coefficient and base-learner
        ens[m] <- xselect
        ensss[[m]] <- coefs[[xselect]]

        ## save updated parameters in x[[xselect]]
        x[[xselect]] <- updatecoefs(x[[xselect]], coefs[[xselect]])

        if (risk == "inbag"){
            mrisk[m] <- - logLH[[xselect]](coefs[[xselect]] * nu)
            if (trace)
                cat("\trisk (inbag) = ", mrisk[m], "\n")
        }
        if (risk == "oobag"){
            ## make a call to PMLE function to built matrices and get the logLH function (without estimation of coefficients)
            dummy <- PMLE(oob$y, oob$x, offset, fit_oob, ens[1:m], nu, maxit, estimation = FALSE)
            mrisk[m] <- - dummy$logLH(coefs[[xselect]] * nu)
            if (trace)
                cat("\trisk (oobag) = ", mrisk[m], "\n")

            ## update step for oob data
            oob$x[[xselect]] <- updatecoefs(oob$x[[xselect]], coefs[[xselect]])
            fit_oob <- fit_oob + nu *  oob$x[[xselect]] %*% coefs[[xselect]]
        }
      }
      if (hardStop | risk != "oobag" | which.min(mrisk) < (mstop - mstop * 0.2/hSi) | hSi == 5) break

      ## else if minimum is at the end of the boosting algorithm,
      ## don't stop but proceed with boosting: Therefore, increase mstop
      ## and the length of the arrays and vectors for storage of results
      warning("risk in test sample seems not to be minimal. ", sQuote("mstop"), " increased.")

      hSi <- hSi + 1
      increase <- 100
      ens <-  c(ens, rep(NA, increase))
      dummy <- vector("list", length = (mstop + increase))
      dummy[1:mstop] <- ensss
      ensss <- dummy
      mrisk <- c(mrisk, rep(NA, increase))
      df_est <- rbind(df_est, matrix(NA, nrow = increase, ncol = ncol(df_est)))# matrix of estimated degrees of freedom
      mstart <- mstop + 1
      mstop <- control$mstop <- mstop + increase
    }

    class(mrisk) <- risk

    RET <- list(data = object,          ### original object
                ensemble = ens,         ### selected base-learners
                ensembless = ensss,     ### list of coefficients in each iteration
                fit = fit,              ### vector of fitted values
                offset = offset,        ### offset
                control = control,      ### control parameters
                response = y,           ### the response variable
                risk = mrisk,           ### the negative maximum log LH
                weights = weights,      ### weights used for fitting
                df = df_est,            ### estimated degrees of freedom for smooth base-learners
                coefs = lapply(x[1:length(x)], getcoefs, nu = nu)  ### coefficients
    )

    RET$predict <- function(newdata = NULL, mstop = mstop, ...) {
        if (!is.null(newdata)) {
            if (is.null(colnames(newdata)))
                stop("missing column names for ", sQuote("newdata"))
            if (is.matrix(newdata)) newdata <- as.data.frame(newdata)
        }
        lp <- offset
        for (m in 1:mstop)
            lp <- lp + nu * predict(x[[ens[m]]], newdata = newdata, newcoefs = ensss[[m]])
        return(lp)
    }

    class(RET) <- c("cfboost")
    return(RET)
}
