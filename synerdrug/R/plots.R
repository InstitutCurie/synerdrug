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

matrixPlot <- function(d, co, ro, addVal=FALSE, contour=FALSE){
    g <- ggplot(d, aes_string(x = co, y = ro, fill = "value", z = "value")) + geom_tile()
    g <- g + scale_fill_gradientn(colours = rev(brewer.pal(11, "Spectral")), breaks = seq(0,1, by = 0.2), limits = c(0,1))
    if (addVal) g <- g + geom_text(aes_string(x = co, y = ro, label = "round(value, digits=2)"))
    if (contour) g <- g + stat_contour(breaks = seq(0,1,0.1), linetype = 2)
    return(g)
}

plotHeatmap <- function(object, type = c("exp", "HSA", "Chou", "Bliss", "Loewe", "LoeweExcess"), addVal = TRUE, contour = FALSE, main = ""){
    data <- meanData(object)
    d <- data[, c("Pos1", "Pos2", type)]
    colnames(d)[3] <- "value"
    g <- matrixPlot(d, "Pos1", "Pos2", addVal = addVal, contour = contour)
    drugs <- drugNames(object)
    g <- g + xlab(drugs[1]) + ylab(drugs[2])
    g <- g + ggtitle(main)
    g <- g + theme4Heatmap()
    ## TODO: gérer les couleurs en fct du type
    return(g)
}

plotSurface <- function(object, ...) {
    data <- meanData(object)
    doses <- doses(object)
    drugs <- drugNames(object)
    xo <- mySeq(1:max(data[, "Pos1"]), n = 10)
    yo <- mySeq(1:max(data[, "Pos2"]), n = 10)
    pos <- expand.grid(xo, yo)
    z <- interpol(x = 1:max(data[, "Pos1"]), y = 1:max(data[, "Pos2"]), z = matrix(data[, "x"], ncol = length(doses[[2]])), xo = pos$Var1, yo = pos$Var2)
    z <- matrix(z, ncol = length(yo))
    nbcol <- 11
    color <- rev(brewer.pal(nbcol, "Spectral"))
    nrz <- nrow(z)
    ncz <- ncol(z)
    zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
    facetcol <- cut(zfacet, nbcol)
    surf <- list(x = xo / max(xo), y = yo/max(yo), z = z)
    persp(surf, col = color[facetcol], box = TRUE, axes = TRUE, border = "gray20", scale = FALSE, xlab = drugs[1], ylab = drugs[2], zlab = content(object), ...)
}


parPlot <- function(object, ref = 1){
    d <- meanData(object)
    drugNames <- drugNames(object)
    B <- drugNames[ref]
    A <- drugNames[-ref]
    d[, B] <- as.factor(d[, B])
    g <- ggplot(d, aes_string(x = A, y = "x", color = B, group = B)) + geom_point() + geom_line() + scale_x_log10()
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







