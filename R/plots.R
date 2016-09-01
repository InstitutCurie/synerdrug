## utils ----

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

## heatmap ----

#' Plot drugSyn object as heatmap
#'
#' @param object a DrugSyn object
#' @param what character containing the name of the variable to plot
#' @param addVal logical indicating if values should be displayed on heatmap
#' @param main character containing title of the plot
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#' data(exDrugSyn)
#' ## Experimental data
#' plotHeatmap(exDrugSyn)
#' ## Loewe predicted effect
#' plotHeatmap(exDrugSyn, what = "Loewe")
#' ## Chou CI
#' plotHeatmap(exDrugSyn, what = "Chou")
plotHeatmap <- function(object, what = c("value", "HSA", "Chou", "Bliss", "Loewe", "LoeweExcess"), addVal = TRUE, main = ""){
    what <- match.arg(what, c("value", "HSA", "Chou", "Bliss", "Loewe", "LoeweExcess"))
    data <- meanData(object)
    d <- data[, c("Pos1", "Pos2", what)]
    colnames(d)[3] <- "value"
    g <- matrixPlot(d, "Pos2", "Pos1", addVal = addVal, what = what)
    drugs <- drugNames(object)
    g <- g + xlab(drugs[2]) + ylab(drugs[1])
    g <- g + ggtitle(main)
    g <- g + theme4Heatmap()
    doses <- doses(object)
    g <- g + coord_cartesian(expand = FALSE) + scale_x_continuous(breaks = 1:length(doses[[2]]), labels = doses[[2]]) + scale_y_continuous(breaks = 1:length(doses[[1]]), labels = doses[[1]])
    return(g)
}

## surface ----

#' Plot drugSyn object as surface
#'
#' @param object a \code{DrugSyn} object
#' @param what data to plot
#' @param gridDensity density of the grid for interpolation
#' @param ... other arguments passed to \code{\link{persp}}
#'
#' @return a plot
#' @export
#'
#' @examples
#' #' data(exDrugSyn)
#' ## Experimental data
#' plotSurface(exDrugSyn)
#' ## change azimuth
#' plotSurface(exDrugSyn, theta = -80)
#' ## Loewe excess
#' plotSurface(exDrugSyn, what = "LoeweExcess")
plotSurface <- function(object, what = c("value", "Loewe", "LoeweExcess", "HSA", "Bliss", "Chou"), gridDensity = 4, ...) {
    what <- match.arg(what, c("value", "Loewe", "LoeweExcess", "HSA", "Bliss", "Chou"))
    data <- meanData(object)
    doses <- doses(object)
    drugs <- drugNames(object)
    yo <- mySeq(1:max(data[, "Pos1"]), n = gridDensity)
    xo <- mySeq(1:max(data[, "Pos2"]), n = gridDensity)
    pos <- expand.grid(xo, yo)
    z <- interpol(x = 1:max(data[, "Pos2"]), y = 1:max(data[, "Pos1"]), z = matrix(data[, what], ncol = length(doses[[1]])), xo = pos$Var1, yo = pos$Var2)
    z <- matrix(z, ncol = length(yo))
    nbcol <- 11
    if (what %in% c("value", "Loewe", "LoeweExcess")) {
        color <- rev(brewer.pal(nbcol, "Spectral"))
    } else color <- scales::seq_gradient_pal("#132B43", "#56B1F7")(seq(0,1, length.out = nbcol))
    nrz <- nrow(z)
    ncz <- ncol(z)
    zfacet <- (z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]) / 4
    if (what %in% c("value", "Loewe")) {
        zfacet <- trimvalues(zfacet, c(0, 1))
        facetcol <- cut(zfacet, seq(0, 1, length.out = 12))
    } else if (what == "LoeweExcess") {
        zfacet <- trimvalues(zfacet, c(-0.3, 0.3))
        facetcol <- cut(zfacet, seq(-0.3, 0.3, length.out = 12))
    } else {
        zfacet <- trimvalues(zfacet, c(0.5, 1.5))
        facetcol <- cut(zfacet, seq(0.5, 1.5, length.out = 12))
    }
    z[!is.finite(z)] <- NA
    surf <- list(x = xo / max(xo), y = yo / max(yo), z = z)

    ## overwrite default args by ...
    args <- list(x = surf, col = color[facetcol], box = TRUE, axes = TRUE, border = "gray30", scale = FALSE, xlab = drugs[2], ylab = drugs[1], zlab = content(object))
    args <- modifyList(args, list(...))
    #persp(surf, col = color[facetcol], box = TRUE, axes = TRUE, border = "gray30", scale = FALSE, xlab = drugs[2], ylab = drugs[1], zlab = content(object), ...)
    do.call(persp, args)
}


## parallel lines ----

#' Plot DrugSyn object as parallel lines
#'
#' @param object a \code{DrugSyn} object
#' @param ref character, name of the drug plotted as x-axis
#' @param xaxis either "pos" or "dose", indicates if the scale of x-axis is in doses of positions in data matrix
#'
#' @return a \code{ggplot} object
#' @export
#'
#' @examples
#' data(exDrugSyn)
#' parPlot(exDrugSyn, ref = "drugA")
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

## individual responses ----

#' Plot individual drug response
#'
#' @param object a \code{DrugSyn} object
#' @param drug character, name of the drug
#' @param xlog boolean, describe if th x-axis is in log-scale
#'
#' @return a \code{ggplot} object
#' @export
#'
#' @examples
#' data(exDrugSyn)
#'respPlot(exDrugSyn, "drugA")
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

## median effect ----

## median effect plot
medianPlot <- function(object, ref){
    d <- meanData(object)
    d[, "fa"] <- log(d[, "x"] / (1 - d[, "x"]))
    B <- "B"
    ggplot(d, aes_string(x = "A", y = "fa", color = factor(B))) + geom_point() + geom_line() + scale_x_log10()
    }

## isobologram ----

#' Plot an isobologram
#'
#' @param object a \code{DrugSyn} object
#' @param effect numeric, effect at which the isobologram is computed
#' @param mode type of the isobologram, see Details
#'
#' @return a plot
#' @details The \code{mode} argument defines the isobologram type. Use \code{linear} for classical isobologram where additivity is represented as a line; use \code{SP} to plot isobologram as defined by Steel and Peckham where additivity is represented as an enveloppe.
#' @export
#'
#' @examples
#' data(exDrugSyn)
#' ## linear isobologram
#' isobologram(exDrugSyn)
#' ## Steel and Peckham isobologram
#' isobologram(exDrugSyn, mode = "SP")
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
        f2A <- function(x, Beq) abs(effect - mod2ll4(respInd[[2]])(Beq + x))
        Bpos <- lapply(Beq, function(B) optimize(f2A, interval = range(doses[[2]]), Beq = B))
        Bpos <- sapply(Bpos, "[[", 1)
        Bpos[!is.finite(Beq)] <- NA
        g <- g + geom_line(aes_string(x = "x", y = "y", linetype = "type"), data = data.frame(x = Bpos, y = Apos, type = "II"))
        ## mode 2B
        Bpos <- seq(0, ED2, length.out = 30)
        Beffect <- mod2ll4(respInd[[2]])(Bpos)
        Aeq <- EDpred(Beffect, respInd[[1]])
        f2B <- function(x, Aeq) abs(effect - mod2ll4(respInd[[1]])(Aeq + x))
        Apos <- lapply(Aeq, function(A) optimize(f2B, interval = range(doses[[1]]), Aeq = A))
        Apos <- sapply(Apos, "[[", 1)
        Apos[!is.finite(Aeq)] <- NA

        g <- g + geom_line(aes_string(x = "x", y = "y", linetype = "type"), data = data.frame(x = Bpos, y = Apos, type = "II"))
        ## mode 1
        Apos <- seq(0, ED1, length.out = 30)
        Aeffect <- mod2ll4(respInd[[1]])(Apos)
        deltaEffect <- effect - Aeffect
        Bpos <- EDpred(deltaEffect, respInd[[2]])
        g <- g + geom_line(aes_string(x = "x", y = "y", linetype = "type"), data = data.frame(x = Bpos, y = Apos, type = "I"))
        g <- g + scale_linetype_manual(name = "mode", values = c("I" = "dotted", "II" = "dashed"))
    }

    ## observed data
    ##ypos <- seq(1, predict(invscaleA, ED1)$y, length.out = 20)
    ypos <- predict(invscaleA, seq(0, ED1, length.out = 20))$y

    fun2Optim <- function(x, d, y, effect){
        abs(interpol(x = 1:max(d[, "Pos1"]), y = 1:max(d[, "Pos2"]), z = matrix(d[, "value"], ncol = max(d[, "Pos2"])), xo = y, yo = x) - effect)
    }

    resOptim <- lapply(ypos, function(y, d, effect, doses) optimise(fun2Optim, interval = 1:length(doses[[1]]), d = d, y = y, effect = effect), d = d, effect = effect, doses = doses)

    id <- which(sapply(resOptim, "[", 2) < 0.01)
    xpos <- unlist(sapply(resOptim, "[", 1)[id])
    ypos <- ypos[id]

    xdose <- predict(scaleB, xpos)$y
    ydose <- predict(scaleA, ypos)$y
    ## restrict isobole values to [0, ED50]
    xdose <- trimvalues(xdose, c(0, ED2))
    ydose <- trimvalues(ydose, c(0, ED1))

    dobs <- data.frame(x = xdose, y = ydose)
    g <- g + geom_line(aes_string(x = "x", y = "y"), data = dobs) + geom_point(aes_string(x = "x", y = "y"), data = dobs)

    g <- g + xlab(drugNames(object)[2]) + ylab(drugNames(object)[1])
    return(g)
}



