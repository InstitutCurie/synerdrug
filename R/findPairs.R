## find paisrs of doses to obtain given effect
findPairs <- function(y, d, val, rangeRoot){
    fun2min <- function(x, d, y, val){
        interp(d$PosB, d$PosA, d$value, duplicate = "mean", xo = x, yo = y, linear = TRUE)$z - val
    }
    uniroot(fun2min, rangeRoot, d = d, y = y, val = val)$root
}

tryFindPairs <- function(y, d, val, rangeRoot){
    t <- try(findPairs(y, d, val, rangeRoot))
    if (class(t) == "try-error") res <- NA
    else res <- t
    return(res)
}
