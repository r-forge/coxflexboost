#########################
## Linear Base-Learner ##
#########################

bolsTime <- function(...){
    bols(..., timedep=TRUE)
}

bols <- function(x, z = NULL, xname = NULL, zname = NULL, timedep=FALSE) {

    if (is.null(xname)) xname = deparse(substitute(x))
    if (is.null(zname)) zname = deparse(substitute(z))

    cc <- mboost:::complete_cases(x = x, z = z)

    newX <- function(x, z = NULL, na.rm = TRUE) {
        if (na.rm) {
            x <- x[cc]
            if (!is.null(z))
                z <- z[cc]
        }
        X <- model.matrix(~ x)
        if (any(!cc) & !na.rm) {
            Xtmp <- matrix(NA, ncol = ncol(X), nrow = length(cc))
            Xtmp[cc,] <- X
            X <- Xtmp
        }
        if (!is.null(z)) X <- X * z
        X
    }
    X <- newX(x, z)
    Xna <- X
    if (any(!cc))
        Xna <- newX(x, z, na.rm = FALSE)

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
    attr(X, "timedep") <- timedep
    attr(X, "coefs") <- rep(0, ncol(X))
    attr(X, "predict") <- predictfun
    class(X) <- c("baselearner", "bols")
    return(X)
}


