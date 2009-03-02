############################################
## Integration Using the Trapezoidal Rule ##
## Adapted for Fisher Matrix              ##
############################################

integr_fisher <- function(x, coefs, desMat, predictions = list(), controls = list(), ...){
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
        predictions$td <- predictions$td * exp(foo)
    } else {
        predictions$tconst <- predictions$tconst * exp(foo)
    }

    if(attr(x, "timedep")){
        ####
        ## trapezoidal rule if the _current_ base-learner is time-dependent
        ####
        sub <- length(controls$grid[[1]])
        dummy <- vector("list", length(controls$grid))
        ## for all observations
        for (i in 1:length(controls$grid)){
            ## for all times in time-grid
            forward <- sub * (i-1) # shift after the first iteration from 1 --> 101 etc.
            dummy[[i]] <- matrix(0, ncol = ncol(desMat), nrow = ncol(desMat))
            for (j in 1:sub){
                MAT <- desMat[j + forward,] %o% desMat[j + forward,] * predictions$td[i,j]
                # trapezoidal rule applied if current base-learner is time-dependent
                if (j == 1 || j == sub){
                    dummy[[i]] <- dummy[[i]] + 0.5 * MAT
                } else {
                    dummy[[i]] <- dummy[[i]] + MAT
                }
            }
            dummy[[i]] <- dummy[[i]] * controls$trapezoid_width[i]
        }
    } else {
        ####
        ## trapezoidal rule if the _current_ base-learner is _not_ time-dependent
        ####
        dummy <- vector("list", nrow(desMat))
        for (i in 1:nrow(desMat)){
            MAT <- desMat[i,] %o% desMat[i,]
            dummy[[i]] <- MAT * predictions$tconst[i,]
        }
        # integrating over time for time dependent part by applying the trapezoid rule
        # if controls$grid != NULL, i.e. time-dependent base-learners present
        # returning exp(0) * time as integral if no time-dependent base-learners are present
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
    }

    ## multiplying time dependent integral with constants
    if(attr(x, "timedep")){
        output <- matrix(0, ncol = ncol(desMat), nrow = ncol(desMat))
        predictions$tconst <- predictions$offset * predictions$tconst * rep(1, length(dummy)) # last part only to insure that predictions$tconst is a vector of length dummy
        for (i in 1:length(dummy)){
            output <- output + predictions$tconst[i] * dummy[[i]]
        }

    } else {
        output <- matrix(0, ncol = ncol(desMat), nrow = ncol(desMat))
        predictions$td <- predictions$offset * predictions$td
        for (i in 1:length(dummy)){
             output <- output + predictions$td[i] * dummy[[i]]
        }
    }

    return(output)
}
