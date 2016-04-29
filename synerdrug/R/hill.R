#' Hill function
#'
#' @param x dose to be evaluated
#' @param Emax Emax parameter of hill function
#' @param E0 \eqn{E_0} parameter of hill function
#' @param n slope parameter of hill function
#' @param EC50 EC50 parameter of hill function

#'
#' @return numeric
#' @export
#'
#' @examples
#' hill(0.5, 100, 0, 2, 0.5)
hill <- function(x, Emax, E0, n, EC50){
    res <- (Emax - E0) / (1 + (EC50 / x) ^ n) + E0
    return(res)
}

#' Inverse function of Hill equation
#'
#' @param E effect
#' @param Emax Emax parameter of hill function
#' @param E0 \eqn{E_0} parameter of hill function
#' @param n slope parameter of hill function
#' @param EC50 EC50 parameter of hill function
#'
#' @return numeric
#' @export
#'
#' @examples
#' revhill(50, 100, 0, 2, 0.5)
revhill <- function(E, Emax, E0, n, EC50){
    d <- (Emax - E0) / (E - E0) - 1
    d <- d ^ (1 / n)
    x <- EC50 / d
    return(x)
}


EDpred <- function(E, resnls){
    revhill(E,  coefficients(resnls)[1], coefficients(resnls)[2], coefficients(resnls)[4], coefficients(resnls)[3])
}


estimateHill <- function(d){
    mod <- nlsLM(value ~ E0 + (Emax - E0) / (1 + (EC50 / conc) ^ n),
                start = list(Emax = 1, E0 = 0, EC50 = median(d[,"conc"]), n = 1),
                data = d,lower=c(0,0,0,-Inf))
    return(mod)
}






