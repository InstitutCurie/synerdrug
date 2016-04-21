## Loewe simple, resolution de equation a/A+b/B=1
resLoewe <- function(a, b, EmaxA, EmaxB, E0A, E0B, nA, nB, EC50A, EC50B){
    fun2min <- function(x, a, b, EmaxA, EmaxB, E0A, E0B, nA, nB, EC50A, EC50B){
        a / (revhill(x, EmaxA, E0A, nA, EC50A)) + b / (revhill(x, EmaxB, E0B, nB, EC50B)) - 1
    }
    res <- uniroot(fun2min, sort(c(E0A, E0B, EmaxA, EmaxB))[2:3],
                   a = a, b = b, EmaxA = EmaxA, EmaxB = EmaxB, E0A = E0A,
                   E0B = E0B, nA = nA, nB = nB, EC50A = EC50A, EC50B = EC50B,
                   extendInt = "yes")$root
    return(res)
}

#' Compute the predicted effect under Loewe model
#'
#' @param a dose of drug A
#' @param b dose of drug B
#' @param modA hill model for drug A
#' @param modB hill model for drug B
#'
#' @return expected effect of combination under Loewe additive model
#' @export
#'
#' @examples NULL
loewe <- function(a, b, modA, modB){
   resLoewe(a, b, EmaxA = coefficients(modA)["Emax"], EmaxB = coefficients(modB)["Emax"],
             E0A = coefficients(modA)["E0"], E0B = coefficients(modB)["E0"],
             nA = coefficients(modA)["n"], nB = coefficients(modB)["n"],
             EC50A = coefficients(modA)["EC50"], EC50B = coefficients(modB)["EC50"])
}


tryLoewe <- function(a, b, modA, modB){
    t <- try(loewe(a, b, modA, modB), silent = TRUE)
    if (class(t) == "try-error") res <- NA
    else res <- t
    return(res)
}
