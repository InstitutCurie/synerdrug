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


## generate sequence of doses
mySeq <- function(doses, n){
    d <- data.frame(from = doses[-length(doses)], to = doses[-1])
    return(unique(as.vector(apply(d, 1, function(x, n) seq(x[1], x[2], length.out = n), n))))
}



## convert mod object to hill function
mod2ll4 <- function(coef){
    ##coef <- coefficients(mod)
    f <- function(x) hill(x, coef[1], coef[2], coef[4], coef[3])
    return(f)
}




## interpolate values on surface, inspired by fields::interp.surface
##
## @param x
## @param y
## @param z
## @param xo
## @param yo
interpol <- function(x, y, z, xo, yo){
    nx <- length(x)
    ny <- length(y)
    lx <- approx(x, 1:nx, xo)$y
    ly <- approx(y, 1:ny, yo)$y
    lx1 <- floor(lx)
    ly1 <- floor(ly)
    ex <- lx - lx1
    ey <- ly - ly1
    ex[lx1 == nx] <- 1
    ey[ly1 == ny] <- 1
    lx1[lx1 == nx] <- nx - 1
    ly1[ly1 == ny] <- ny - 1
    interp <- z[cbind(lx1, ly1)] * (1 - ex) * (1 - ey) +
        z[cbind(lx1 + 1, ly1)] * ex * (1 - ey) +
        z[cbind(lx1, ly1 + 1)] * (1 - ex) * ey +
        z[cbind(lx1 + 1, ly1 + 1)] * ex * ey
    return(interp)
}


## check typeHill value
checkTypeHill <- function(x){
    int <- intersect(c(4, 3), x)
    if (length(int) == 0) stop("typeHill must be 3 or 4")
    return(int[1])
}

## trim extreme values
trimvalues <- function(x, limits) {
    x <- ifelse(x < limits[1], limits[1], x)
    x <- ifelse(x > limits[2], limits[2], x)
    return(x)
}

