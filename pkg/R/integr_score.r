############################################
## Integration Using the Trapezoidal Rule ##
## Adapted for Score Function             ##
############################################

integr_score <- function(x, coefs, desMat, predictions = list(), controls = list(), ...){
    ###
    # x             currently added base-learner
    # coefs         the current coefficients
    # predictions   list of 'offset', 'tconst', 'td'
    #                 offset           is the offset
    #                 tconst      is the prediction for time-constant base-learners
    #                 td          is the prediction for time-dependent base-learners
    # controls      list of 'grid' , 'trapezoid_width', 'upper', 'nu', 'which'
    #                 grid             a (time) grid
    #                 trapezoid_width  the distances in grid
    #                 nu               fraction of fit to include
    #                 upper            vector (here: of times t[i])

    if (length(predictions$offset) != 1) stop(sQuote("offset"), " must be a single constant")

    coefs <- coefs + controls$nu * attr(x, "coefs")
    foo <- desMat %*% coefs
    if(attr(x, "timedep")){
        foo <- matrix(foo, nrow = length(controls$grid), ncol = length(controls$grid[[1]]), byrow = TRUE)
        ## component-wise multiplication
        predictions$td <- predictions$td * exp(foo)
    } else {
        ## component-wise multiplication
        predictions$tconst <- predictions$tconst * exp(foo)
    }

    if(attr(x, "timedep")){
        dummy <- vector("list", ncol(desMat))
        for (i in 1:ncol(desMat)){
            dummy[[i]] <- matrix(desMat[,i], nrow = length(controls$grid), ncol = length(controls$grid[[1]]), byrow = TRUE)
            dummy[[i]] <- dummy[[i]] * predictions$td  # component-wise multiplication of matrices
        }
        predictions$td <- dummy
    } else {
        predictions$tconst <- desMat * as.vector(predictions$tconst)
    }

    ###
    ## integrating over time for time dependent part by applying the trapezoid rule
    ## if controls$grid != NULL, i.e. time-dependent base-learners present
    ## returning exp(0) * time as integral if no time-dependent base-learners are present

    ## returned objects:
    ## 1.) time constant base-learner currently added --> automatically returned as matrix with
    ##     nrow = # observations and ncol = # parameters used for base-learner
    ## 2.) time dependent base-learner currently added --> by sapply we get a matrix with same
    ##     specifications as 1.)

    if (!is.null(controls$grid)){
        sub <- length(controls$grid[[1]])
        if (!is.list(predictions$td)){
            predictions$td <- controls$trapezoid_width * (0.5 * (predictions$td[,1] + predictions$td[,sub]) + apply(predictions$td[,2:(sub-1)], 1, sum))
        } else {
            for (i in 1:length(predictions$td))
                predictions$td[[i]] <- controls$trapezoid_width * (0.5 * (predictions$td[[i]][,1] + predictions$td[[i]][,sub]) + apply(predictions$td[[i]][,2:(sub-1)], 1, sum))
            predictions$td <- sapply(predictions$td, cbind)
        }
    } else {
        predictions$td <- controls$upper
    }

    ## multiplying time dependent integral with constants
    int <- predictions$offset * predictions$tconst * predictions$td
    return(apply(int, 2, sum))
}
