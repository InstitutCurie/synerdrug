fun2Optim <- function(x, d, y, val){
    abs(interp(d$PosB, d$PosA, d$value, duplicate="mean", xo=y, yo=x, linear=TRUE)$z - val)
}
