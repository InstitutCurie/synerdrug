fun2Optim <- function(x, d, y, effect){
    abs(interpol(x = 1:max(d[, "Pos1"]), y = 1:max(d[, "Pos2"]), z = matrix(d[, "value"], ncol = max(d[, "Pos2"])), xo = y, yo = x) - effect)
}

