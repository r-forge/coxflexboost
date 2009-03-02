###########################################
## Function to (Randomly) Sample Data    ##
## According to an Arbitrary Hazard Rate ##
###########################################

rSurvTime <- function(lambda, x, cens_fct, upper = 1000, ..., file = NULL){
    ###
    # lambda        function. Baseline hazard \lambda(t, x) (time must be first argument)
    # x             matrix. (sampled) values for covariates
    # cens_fct      function. Function to compute (random) censoring
    # upper         upper boundary of the interval the random survival times fall into
    # file          name of the data file the generated data set should be stored into
    #               (e.g., "survtimes.RData") or NULL if the dataset should directly be
    #               returned in R

    if (!is.matrix(x)) x <- cbind(x)
    time <- rep(NA, nrow(x))

    Lambda <- function(lambda, x, time){
        integrate(lambda, 0, time, x = x)$value
    }

    InvLambda <- function(Lambda, lambda, x){
        negLogU <- - log(runif(1, 0, 1))
        rootfct <- function(time) {negLogU - Lambda(lambda, x, time)}
        return(uniroot(rootfct, interval = c(0, upper))$root)
    }

    for (i in 1:nrow(x)){
        time[i] = InvLambda(Lambda, lambda, x[i,])
    }

    time_event = cens_fct(time, ...)

    data = data.frame(time = time_event[,1], event = time_event[,2], x = x)
    if (!is.null(file)){
        save(data, file = file)
        invisible(data)
    } else {
        return(data)
    }
}
