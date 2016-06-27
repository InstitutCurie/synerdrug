## Loewe simple, resolution de equation a/A+b/B=1
resLoewe <- function(a, b, EmaxA, EmaxB, E0A, E0B, nA, nB, EC50A, EC50B){
    fun2min <- function(x, a, b, EmaxA, EmaxB, E0A, E0B, nA, nB, EC50A, EC50B){
        A <- a / (revhill(x, EmaxA, E0A, nA, EC50A))
        if (!is.finite(A)) A <- 0
        B <- b / (revhill(x, EmaxB, E0B, nB, EC50B))
        if (!is.finite(B)) B <- 0
        return(A + B - 1)
    }
    res <- uniroot(fun2min, sort(c(E0A, E0B, EmaxA, EmaxB))[c(2, 3)] + c(1e-9, -1e-9),  ## add eps to range to compute f at extremas
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
   resLoewe(a, b, EmaxA = modA["Emax"], EmaxB = modB["Emax"],
             E0A = modA["E0"], E0B = modB["E0"],
             nA = modA["n"], nB = modB["n"],
             EC50A = modA["EC50"], EC50B = modB["EC50"])
}


tryLoewe <- function(a, b, modA, modB){
    t <- try(loewe(a, b, modA, modB), silent = TRUE)
    if (class(t) == "try-error") res <- NA
    else res <- t
    return(res)
}
