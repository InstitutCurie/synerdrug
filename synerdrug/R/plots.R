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




