############################################
## Integration Using the Trapezoidal Rule ##
############################################

integr <- function(x, coefs, desMat, predictions = list(), controls = list(), ...){
    ###
    # x             currently added base-learner
    # coefs         the current coefficients
    # predictions   list of 'offset', 'tconst', 'td'
    #                 offset           is the exp(offset)
    #                 tconst           is the exp(prediction for time-constant base-learners)
    #                 td               is the exp(prediction for time-dependent base-learners)
    # controls      list of 'grid' , 'trapezoid_width', 'upper', 'nu', 'which'
    #                 grid             a (time) grid
    #                 trapezoid_width  the distances in grid
    #                 nu               fraction of fit to include
    #                 upper            vector (here: of times t[i])
    #                 which            can be "NULL", "score", "fisher" returning the integral needed
    #                                  for the log likelihood ("NULL"), score vector or fisher matrix

    if (length(predictions$offset) != 1) stop(sQuote("offset"), " must be a single constant")

    coefs <- coefs + controls$nu * attr(x, "coefs")
    foo <- desMat %*% coefs
    if(attr(x, "timedep")){
        foo <- matrix(foo, nrow = length(controls$grid), ncol = length(controls$grid[[1]]), byrow = TRUE)
        predictions$td <- predictions$td * exp(foo)
    } else {
        predictions$tconst <- predictions$tconst * exp(foo)
    }

    ###
    # integrating over time for time dependent part by applying the trapezoid rule
    # if controls$grid != NULL, i.e. time-dependent base-learners present
    # returning exp(0) * time as integral if no time-dependent base-learners are present
    if (!is.null(controls$grid)){
        sub <- length(controls$grid[[1]])
        predictions$td <- controls$trapezoid_width * (0.5 * (predictions$td[,1] + predictions$td[,sub]) + apply(predictions$td[,2:(sub-1)], 1, sum))
    } else {
        predictions$td <- controls$upper
    }

    ###
    # multiplying time dependent integral with constants
    return(predictions$offset * predictions$tconst * predictions$td)
}
