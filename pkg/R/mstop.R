######################################################
## Function to Extract (Optimal) Stopping Iteration ##
######################################################

mstop <- function(object, ...)
    UseMethod("mstop")

mstop.cfboost <- function(object, opt = TRUE, ...){
    if(!opt){
        #warning("returned value is not the optimal mstop")
        return(NROW(object$ensemble))
    } else {
        mr <- print(risk(object), print = FALSE)
        if (class(mr) == "inbag")
            warning("computing ", sQuote("mstop"), " on learning sample may lead to overfitting")
        return(mr$iteration)
    }
}
