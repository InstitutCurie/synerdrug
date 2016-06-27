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


EDpred <- function(E, coeffs){
    revhill(E,  coeffs[1], coeffs[2], coeffs[4], coeffs[3])
}


estimateHill <- function(d, nbParam = c(4, 3)){
    nbParam <- checkTypeHill(nbParam)
    if (nbParam == 4) mod <- nlsLM(value ~ E0 + (Emax - E0) / (1 + (EC50 / conc) ^ n),
                start = list(Emax = 1, E0 = 0, EC50 = median(d[,"conc"]), n = 1),
                data = d, lower = c(0, 0, 0, -Inf))
    else mod <- nlsLM(value ~ E0 + (Emax - E0) / (1 + (EC50 / conc) ^ n),
                      start = list(Emax = 1, E0 = 0, EC50 = median(d[,"conc"]), n = 1),
                      data = d, lower = c(0, 0, 0, -Inf), upper = c(+Inf, 0, +Inf, +Inf))
    coeff <- coefficients(mod)
    ## check E0 and Emax values
    if (coeff["Emax"] < coeff["E0"]) {
        coeff[c("Emax", "E0")] <- coeff[c("E0", "Emax")]
        coeff["n"] <- -coeff["n"]
    }
    return(coeff)
}






