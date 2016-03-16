#' Hill function
#'
#' @param x
#' @param Emax
#' @param E0
#' @param n
#' @param EC50
#'
#' @return numeric
#' @export
#'
#' @examples
#' hill(0.5, 100, 0, 2, 0.5)
hill <- function(x, Emax, E0, n, EC50){
    res <- (Emax-E0)/(1+(EC50/x)^n)+E0
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
    d <- (Emax - E0)/(E - E0) - 1
    d <- d^(1/n)
    x <- EC50/d
    return(x)
}
