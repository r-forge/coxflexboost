###########################
## P-Spline Base-Learner ##
###########################

bbsTime <- function(...){
    bbs(..., timedep=TRUE)
}

bbs <- function(x, z = NULL, knots = 20, degree = 3, differences = 2, df = 4,
                   center = FALSE, xname = NULL, zname = NULL, lambda = NULL, timedep = FALSE) {
    ## df determines lambda corresponding to the degrees of freedom in the first iteration

    cc <- mboost:::complete_cases(x = x, z = z)

    if (is.null(xname)) xname <- deparse(substitute(x))
    if (is.null(zname)) zname <- deparse(substitute(z))

    if (is.factor(x) || (df <= 2 && !center))
        return(bols(x = x, z = z, xname = xname, zname = zname))

    if (!differences %in% 1:3)
        stop(sQuote("differences"), " are not in 1:3")
    if ((!center) && (df < differences))
        stop(sQuote("df"), " is less than ", sQuote("differences"))
    if(center && (degree < (differences-1)))
        stop(sQuote("degree"), " is less than ", sQuote("differences"), "-1")
    if (length(unique(x)) < 6)
        stop(sQuote(xname), " has less than 6 unique values")

    n.kn <- function(n) {
        ## Number of inner knots
        if(n < 40) n
        else 40
    }

    if (is.null(knots)) {
        n <- length(x)
        nk <- n.kn(n)
        ### ADDED: maximal 20 knots (to reduce computational burden)
        if (nk > 20){
            warning("Maximal 20 knots allowed. Number of knots set to 20")
            nk <- 20
        }
        knots <- seq(from = min(x, na.rm = TRUE),
                     to = max(x, na.rm = TRUE), length = nk)
        knots <- knots[2:(length(knots) - 1)]
    } else {
        if (length(unique(diff(knots))) > 1)
            warning("non-equidistant ", sQuote("knots"),
                    " might be inappropriate")
    }

    if (length(knots) == 1) {
        ### ADDED: maximal 20 knots (to reduce computational burden)
        if (knots > 20){
            warning("Maximal 20 knots allowed. Number of knots set to 20")
            knots <- 20
        }
        knots <- seq(from = min(x, na.rm = TRUE),
                     to = max(x, na.rm = TRUE), length = knots+2)
        knots <- knots[2:(length(knots) - 1)]
    }
    boundary.knots <- range(x, na.rm = TRUE)

    newX <- function(x, z = NULL, na.rm = TRUE) {
        if (na.rm) {
            x <- x[cc]
            if (!is.null(z))
                z <- z[cc]
        }
        X <- bs(x, knots = knots, degree = degree, intercept = TRUE,
                Boundary.knots = boundary.knots)
        if (!is.null(z))
            X <- X * z
        if (center) {
            K <- diff(diag(ncol(X)), differences = differences)
            X <- tcrossprod(X, K) %*% solve(tcrossprod(K))
        }
        return(X)
    }

    X <- newX(x, z)
    Xna <- X
    if (any(!cc))
        Xna <- newX(x, z, na.rm = FALSE)

    if (center) {
        K <- diag(ncol(X))
    } else {
        K <- diff(diag(ncol(X)), differences = differences)
        K <- crossprod(K, K)
    }

    predictfun <- function(coefs, newdata = NULL) {
        if (is.null(newdata)) return(Xna %*% coefs)
        nX <- newX(x = newdata[[xname]], z = newdata[[zname]], na.rm = FALSE)
        nX %*% coefs
    }

    designMat <- function(newdata = NULL){
        if (is.null(newdata)) return(Xna)
        newX(x = newdata[[xname]], z = newdata[[zname]], na.rm = FALSE)
    }

    attr(X, "designMat") <- designMat
    attr(X, "df") <- df
    attr(X, "lambda") <- lambda
    if (is.null(lambda)) lambda <- 1
    attr(X, "pen") <- lambda * K
    ### ATTENTION: pen is assumed to be lambda * K
    ### so one needs make a call to df2lambda and store this in lambda and pen
    attr(X, "timedep") <- timedep
    attr(X, "coefs") <- rep(0, ncol(X))
    attr(X, "predict") <- predictfun
    class(X) <- c("baselearner", "bbs")
    return(X)
}
