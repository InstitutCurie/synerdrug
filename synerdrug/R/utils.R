#' Convert character to numeric
#'
#' @param x string to convert
#'
#' @return numeric
c2num <- function(x){
    as.numeric(gsub(",", ".", x, fixed = TRUE))
}

#' Automatic detection of decimal character in file
#'
#' @param file input file
#'
#' @return character "." or ","
#'
detectDec <- function(file){
    d <- readLines(file, n = 3)
    npoint <- length(grep("\\.", d))
    ncomma <- length(grep(",", d))
    return(ifelse(npoint > ncomma, ".", ","))
}


#' Generate sequence of doses
#'
#' @param doses data doses
#' @param n number of points in each doses interval (inclusive)
#'
#' @return sequence of doses
#'
mySeq <- function(doses, n){
    d <- data.frame(from = doses[-length(doses)], to = doses[-1])
    return(unique(as.vector(apply(d, 1, function(x, n) seq(x[1], x[2], length.out = n), n))))
}

