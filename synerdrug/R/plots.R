themeNobg <- function(){
    theme(plot.background = element_blank(),
          legend.background = element_blank())
}

theme4Heatmap <- function(){
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none"
    ) + themeNobg()
}

matrixPlot <- function(d, co, ro, addVal = FALSE, what = c("value", "HSA", "Chou", "Bliss", "Loewe", "LoeweExcess")){
    d[, "valueFill"] <- d[, "value"]
    if (what %in% c("value", "Loewe")) {
        d[, "valueFill"] <- trimvalues(d[, "valueFill"], c(0, 1))
        g <- ggplot(d, aes_string(x = co, y = ro, z = "value"))
                g <- g + geom_tile(aes_string(fill = "valueFill"), color = "gray40")
        g <- g + scale_fill_gradientn(colours = rev(brewer.pal(11, "Spectral")), breaks = seq(0, 1, by = 0.2), limits = c(0, 1))
        if (addVal) g <- g + geom_text(aes_string(x = co, y = ro, label = "round(value, digits=2)"))

    }
    if (what == "LoeweExcess") {
        d[, "valueFill"] <- trimvalues(d[, "valueFill"], c(-0.3, 0.3))
        g <- ggplot(d, aes_string(x = co, y = ro, z = "value"))
        g <- g + geom_tile(aes_string(fill = "valueFill"), color = "gray40")
        g <- g + scale_fill_gradientn(colours = rev(brewer.pal(11, "Spectral")), breaks = seq(-0.3, 0.3, by = 0.1), limits = c(-0.3, 0.3))
        if (addVal) g <- g + geom_text(aes_string(x = co, y = ro, label = "round(value, digits=2)"))
    }
    if (what %in% c("HSA", "Chou", "Bliss")) {
        d[, "valueFill"] <- trimvalues(d[, "valueFill"], c(0.5, 1.5))
        g <- ggplot(d, aes_string(x = co, y = ro, z = "value"))
        g <- g + geom_tile(aes_string(fill = "valueFill"), color = "gray40")
        g <- g + scale_fill_gradient(low = "#132B43", high = "#56B1F7", breaks = seq(0.5, 1.5, by = 0.1), limits = c(0.5, 1.5))
        if (addVal) g <- g + geom_text(aes_string(x = co, y = ro, label = "round(value, digits=2)"), color = "white")

    }
    return(g)
}

plotHeatmap <- function(object, what = c("value", "HSA", "Chou", "Bliss", "Loewe", "LoeweExcess"), addVal = TRUE, main = ""){
    what <- match.arg(what, c("value", "HSA", "Chou", "Bliss", "Loewe", "LoeweExcess"))
    data <- meanData(object)
    d <- data[, c("Pos1", "Pos2", what)]
    colnames(d)[3] <- "value"
    g <- matrixPlot(d, "Pos1", "Pos2", addVal = addVal, what = what)
    drugs <- drugNames(object)
    g <- g + xlab(drugs[1]) + ylab(drugs[2])
    g <- g + ggtitle(main)
    g <- g + theme4Heatmap()
    doses <- doses(object)
    g <- g + coord_cartesian(expand = FALSE) + scale_x_continuous(breaks=1:length(doses[[1]]), labels=doses[[1]]) + scale_y_continuous(breaks = 1:length(doses[[2]]), labels=doses[[2]])
    return(g)
}

plotSurface <- function(object, ...) {
    data <- meanData(object)
    doses <- doses(object)
    drugs <- drugNames(object)
    xo <- mySeq(1:max(data[, "Pos1"]), n = 4)
    yo <- mySeq(1:max(data[, "Pos2"]), n = 4)
    pos <- expand.grid(xo, yo)
    z <- interpol(x = 1:max(data[, "Pos1"]), y = 1:max(data[, "Pos2"]), z = matrix(data[, "value"], ncol = length(doses[[2]])), xo = pos$Var1, yo = pos$Var2)
    z <- matrix(z, ncol = length(yo))
    nbcol <- 11
    color <- rev(brewer.pal(nbcol, "Spectral"))
    nrz <- nrow(z)
    ncz <- ncol(z)
    zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
    facetcol <- cut(zfacet, nbcol)
    surf <- list(x = xo / max(xo), y = yo/max(yo), z = z)
    persp(surf, col = color[facetcol], box = TRUE, axes = TRUE, border = "gray30", scale = FALSE, xlab = drugs[1], ylab = drugs[2], zlab = content(object), ...)
}

## parallel plot
parPlot <- function(object, ref, xaxis = c("pos", "dose")){
    xaxis <- match.arg(xaxis, c("pos", "dose"))
    d <- meanData(object)
    drugs <- drugNames(object)
    notRef <- setdiff(drugs, ref)
    d[, notRef] <- as.factor(d[, notRef])
    if (xaxis == "dose") g <- ggplot(d, aes_string(x = ref, y = "value", color = notRef))
    else {
        g <- ggplot(d, aes_string(x = paste0("Pos", match(ref, drugs)), y = "value", color = notRef))
    }
    g <- g + geom_point() + geom_line()
    g <- g + ylab(content(object)) + xlab(ref)
    return(g)
}



## courbe reponse individuelle
respPlot <- function(object, drug, xlog = TRUE){
    if (!(drug %in% drugNames(object))) stop("Wrong drug name")
    d <- expData(object)
    d <- d[d[, setdiff(colnames(d), c("value", "rep", drug))] == 0, ]
    g <- ggplot(d, aes_string(x = drug, y = "value"))
    if (any(colnames(d) == "rep")) {
        g <- g + geom_point(aes(shape = factor(rep))) + scale_shape_discrete(name = "Replicate")
    } else g <- g + geom_point()
    if (xlog) g <- g + scale_x_log10()
    if (length(respInd(object))) {
        mod <- mod2ll4(respInd(object)[[drug]])
        g <- g + stat_function(fun = mod, color = "red")
    }
    return(g)
}


## median effect plot
medianPlot <- function(object, ref){
    d <- meanData(object)
    d[, "fa"] <- log(d[, "x"] / (1 - d[, "x"]))
    ggplot(d, aes(x = A, y = fa, color=factor(B))) + geom_point() + geom_line() + scale_x_log10()
    }



isobologram <- function(object, effect = 0.5, mode = c("linear", "SP")){
    d <- meanData(object)
    mode <- match.arg(mode, c("linear", "SP"))
    drugs <- drugNames(object)
    doses <- doses(object)
    respInd <- respInd(object)
    ## convert matrix to dose space
    scaleA <- smooth.spline(1:length(doses[[1]]), doses[[1]])
    scaleB <- smooth.spline(1:length(doses[[2]]), doses[[2]])
    invscaleA <- smooth.spline(doses[[1]], 1:length(doses[[1]]))

    ## isobole additive
    ED1 <- EDpred(effect, respInd[[drugs[1]]])
    ED2 <- EDpred(effect, respInd[[drugs[2]]])
    dlin <- data.frame(x = c(0,ED2), y = c(ED1, 0), xend = c(ED2, 0), yend = c(0, ED1))
    if (mode == "linear") g <- ggplot() + geom_segment(aes_string(x = "x", y = "y", xend = "xend", yend = "yend"), data = dlin, lty = 2)
    else {
        g <- ggplot()
        ## mode 2A
        Apos <- seq(0, ED1, length.out = 30)
        Aeffect <- mod2ll4(respInd[[1]])(Apos)
        Beq <- EDpred(Aeffect, respInd[[2]])
        f <- function(x, Beq) abs(effect - mod2ll4(respInd[[2]])(Beq + x))
        Bpos <- lapply(Beq, function(B) optimize(f, interval = range(doses[[2]]), Beq = B))
        Bpos <- sapply(Bpos, "[[", 1)
        g <- g + geom_line(aes(x = x, y = y, linetype = "II"), data = data.frame(x = Bpos, y = Apos))
        ## mode 2B
        Bpos <- seq(0, ED2, length.out = 30)
        Beffect <- mod2ll4(respInd[[2]])(Bpos)
        Aeq <- EDpred(Beffect, respInd[[1]])
        f <- function(x, Aeq) abs(effect - mod2ll4(respInd[[1]])(Aeq + x))
        Apos <- lapply(Aeq, function(A) optimize(f, interval = range(doses[[1]]), Aeq = A))
        Apos <- sapply(Apos, "[[", 1)
        g <- g + geom_line(aes(x = x, y = y, linetype = "II"), data = data.frame(x = Bpos, y = Apos))
        ## mode 1
        Apos <- seq(0, ED1, length.out = 30)
        Aeffect <- mod2ll4(respInd[[1]])(Apos)
        deltaEffect <- effect - Aeffect
        Bpos <- EDpred(deltaEffect, respInd[[2]])
        g <- g + geom_line(aes(x = x, y = y, linetype = "I"), data = data.frame(x = Bpos, y = Apos))
        g <- g + scale_linetype_manual(name = "mode", values = c("I" = "dotted", "II" = "dashed"))
    }

    ## observed data
    ypos <- seq(1, predict(invscaleA, ED1)$y, length.out = 20)
    resOptim <- lapply(ypos, function(y, d, effect, doses) optimise(fun2Optim, interval = 1:length(doses[[1]]), d = d, y = y, effect = effect), d = d, effect = effect, doses = doses)

    id <- which(sapply(resOptim, "[", 2) < 0.01)
    xpos <- unlist(sapply(resOptim, "[", 1)[id])
    ypos <- ypos[id]

    xdose <- predict(scaleB, xpos)$y
    ydose <- predict(scaleA, ypos)$y
    dobs <- data.frame(x = xdose, y = ydose)
    g <- g + geom_line(aes(x = x, y = y), data = dobs) + geom_point(aes(x = x, y = y), data = dobs)

    g <- g + xlab(drugNames(object)[2]) + ylab(drugNames(object)[1])
    return(g)
}



