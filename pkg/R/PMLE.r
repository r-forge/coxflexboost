###################
## Penalized MLE ##
###################

## Function for penalized ML estimation
PMLE <- function(y, x, offset, fit, ens, nu, maxit, subdivisions = 100, estimation = TRUE, trace){

    time <- y[,1]
    delta <- y[,2]

    ## currently added base-learner
    added_bl <- ens[length(ens)]
    ## number of coefficients to estimate
    coefs <- rep(0, ncol(x[[added_bl]]))
    ## get penalty for currently added base-learner
    pen <- attr(x[[added_bl]], "pen")
    ## get selected base-learners (uniquely)
    ens_uni <- unique(ens)

    if(any(sapply(x[ens_uni], attr, which = "timedep"))){
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

    ## remove the most recently added base-learner
    ens_uni <- ens_uni[ens_uni != added_bl]

    ## initialization (needed if either no time-dependent or no time-constant base-learner present)
    predictions_tconst = 0
    predictions_td = 0

    ## data pre-processing
    if(length(ens_uni) > 0){
        ## extract time-dependency
        timedep <- sapply(x[ens_uni], attr, which = "timedep")
        if (any(!timedep)){
            ## predictions for time-constant base-learners
            predictions_tconst <- sapply(x[ens_uni[!timedep]], predict)
            predictions_tconst <- apply(predictions_tconst, 1, sum)
        }
        if (any(timedep)){
            ## prediction for time-dependent base-learners
            predictions_td <- sapply(x[ens_uni[timedep]], predict_td, time_grid = grid)
            predictions_td <- apply(predictions_td, 1, sum)
            predictions_td <- matrix(predictions_td, nrow = n, ncol = sub, byrow = TRUE)
        }
    }

    ## build design matrix for currently added base-learner
    if (attr(x[[added_bl]], "timedep")){
        xd <- unlist(grid)
        xname <- get("xname", environment(attr(x[[added_bl]], "predict")))
        zd <- get("z", environment(attr(x[[added_bl]], "predict")))
        if (!is.null(zd)){
            zname <- get("zname", environment(attr(x[[added_bl]], "predict")))
            zd <- rep(zd, each = length(grid[[1]]))
            newdata <- data.frame(cbind(xd, zd))
            names(newdata) <- c(xname, zname)
        } else {
            newdata <- data.frame(xd)
            names(newdata) <- c(xname)
        }
        desMat <- attr(x[[added_bl]],"designMat")(newdata = newdata)
    } else {
        desMat <- attr(x[[added_bl]],"designMat")()
    }

    exp_offset <- exp(offset)
    exp_pred_tconst <- exp(predictions_tconst * nu)
    exp_pred_td <- exp(predictions_td * nu)

    logLH_pen <- function(coefs){
        log_lik <- sum(delta * (fit + x[[added_bl]] %*% coefs)
                       - integr(x[[added_bl]], coefs, desMat,
                                predictions = list(offset = exp_offset, tconst = exp_pred_tconst, td = exp_pred_td),
                                controls = list(grid = grid, trapezoid_width = trapezoid_width, upper = time, nu = nu)))
        if (is.null(pen)) pen <- 0 else pen <- 0.5 * (coefs %*% pen %*% coefs)
        return(log_lik - pen)
    }

    s_pen <- function(coefs){
        int <- integr_score(x[[added_bl]], coefs, desMat,
                            predictions = list(offset = exp_offset, tconst = exp_pred_tconst, td = exp_pred_td),
                            controls = list(grid = grid, trapezoid_width = trapezoid_width, upper = time, nu = nu))
        score_vec <- t(delta) %*% x[[added_bl]] - int
        if (is.null(pen)) pen <- 0 else pen <- (pen %*% coefs)
        return(score_vec - as.vector(pen))
    }

    F_pen <- function(coefs){
        Fisher_mat <- integr_fisher(x[[added_bl]], coefs, desMat,
                                    predictions = list(offset = exp_offset, tconst = exp_pred_tconst, td = exp_pred_td),
                                    controls = list(grid = grid, trapezoid_width = trapezoid_width, upper = time, nu = nu, which = "fisher"))
        if (is.null(pen)) pen <- 0
        return(list(F = Fisher_mat, F_pen = Fisher_mat + pen))
    }

    if (estimation){
        repeat{
            ## optimize (penalized) log-likelihood
            result <- optim(par = coefs, fn = logLH_pen, method = "BFGS", control = list(fnscale = -1, maxit = maxit))

            if (result$convergence == 0){
                if (!is.null(pen)){
                    ## compute (penalized and unpenalized) _observed_ Fisher matrix
                    result$Fisher <- F_pen(result$par)
                    if (is.null(pen)) pen <- 0
                    ## compute estimated degrees of freedom
                    m <- result$Fisher$F  %*% solve(result$Fisher$F_pen)
                    result$df <- sum(diag(m))  # trace of m
                    ## set penalty matrix NULL
                    pen <- NULL
                }
                ## calculate criterion value (i.e., unpenalized log-likelihood)
                result$maxll <- logLH_pen(result$par)
                ## return function for (unpenalized) log-likelihood
                result$logLH <- logLH_pen
                return(result)
            }
            maxit <- maxit * 10
            warning(sQuote("optim"), " did not converge. ",
                    sQuote("maxit"), "was increased. For faster computation please increase ", sQuote("boost_control(maxit)"))
            if (maxit > 5000000) stop(sQuote("optim"), " did not converge.")
        }
    } else {
        pen <- NULL
        ## return function for (unpenalized) log-likelihood
        return(list(logLH = logLH_pen))
    }
}

