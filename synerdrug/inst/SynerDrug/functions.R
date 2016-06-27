## heatmap of data matrix
matrixPlot <- function(d, co, ro, addVal=FALSE, contour=FALSE){
    g <- ggplot(d, aes_string(x=co, y=ro, fill="value", z="value")) + geom_tile()
    g <- g + scale_fill_gradientn(colours=rev(brewer.pal(11, "Spectral")), breaks=seq(0,1, by=0.2), limits=c(0,1))
    if (addVal) g <- g + geom_text(aes_string(x=co, y=ro, label="round(valueOrig, digits=2)"))
    if (contour) g <- g + stat_contour(breaks=seq(0,1,0.1), linetype=2)
    return(g)
}

## parrallel plot
## TODO couleurs à ajuster
parPlot <- function(d, A, B){
    d[, B] <- as.factor(d[, B])
    g <- ggplot(d, aes_string(x=A, y="value", color=B, group=B)) + geom_point() + geom_line()
    return(g)
}

themeNobg <- theme(plot.background = element_blank(),
                   legend.background=element_blank())

theme4Heatmap <-  theme(panel.grid.minor = element_blank(),
                        panel.grid.major = element_blank(),
                        panel.background = element_blank(),
                        axis.line=element_blank(),
                        ## axis.text.x=element_blank(),
                        ## axis.text.y=element_blank(),
                        axis.ticks=element_blank(),
                        legend.position="none"
                        ) + themeNobg 

## modelisation par hill
hill <- function(x, Emax, E0, n, EC50){
    res <- (Emax-E0)/(1+(EC50/x)^n)+E0
    return(res)
}

## predhill <- function(x, resnls){
##     res <- hill(x, coefficients(resnls)[1], coefficients(resnls)[3], coefficients(resnls)[2])
##     return(res)
## }

revhill <- function(E, Emax, E0, n, EC50){
    d <- (Emax - E0)/(E - E0) - 1
    d <- d^(1/n)
    x <- EC50/d
    return(x)
}

EDpred <- function(E, resnls){
    revhill(E,  coefficients(resnls)[1], coefficients(resnls)[2], coefficients(resnls)[4], coefficients(resnls)[3])
}

## courbe réponse individuelle
respPlot <- function(d, drug){
    g <- ggplot(d, aes_string(x=drug, y="value"))
    if (any(colnames(d)=="rep")){
        g <- g + geom_point(aes(shape=factor(rep))) + scale_shape_discrete(name="Replicate")
    } else g <- g + geom_point() 
    return(g)    
}



## conversion char to num
c2num <- function(x){
    as.numeric(gsub(",", ".", x, fixed = TRUE))
}

## ## calcul isobole
## iso <- function(a, E, EmaxA, EmaxB, nA, nB, EC50A, EC50B){
##     B <- revhill(E, EmaxB, nB, EC50B)
##     ga <- 1 + (EC50A^nA)/(a^nA)
##     ga <- (EmaxB/EmaxA) * ga - 1
##     ga <- ga^(1/nB)
##     b <- B - EC50B/ga    
##     return(b)
## }

## isobole <- function(a, E, modA, modB){
##     return(iso(a, E, coefficients(modA)[1],coefficients(modB)[1], coefficients(modA)[3], coefficients(modB)[3], coefficients(modA)[2], coefficients(modB)[2]))
## }


## Loewe simple, résolution de equation a/A+b/B=1
resLoewe <- function(a, b, EmaxA, EmaxB, E0A, E0B, nA, nB, EC50A, EC50B){
    fun2min <- function(x, a, b, EmaxA, EmaxB, E0A, E0B, nA, nB, EC50A, EC50B){
        A <- a / (revhill(x, EmaxA, E0A, nA, EC50A))
        if (!is.finite(A)) A <- 0
        B <- b / (revhill(x, EmaxB, E0B, nB, EC50B))
        if (!is.finite(B)) B <- 0
        return(A + B - 1)
    }
    res <- uniroot(fun2min, sort(c(E0A, E0B, EmaxA, EmaxB))[2:3] + c(1e-9, -1e-9), a=a, b=b, EmaxA=EmaxA, EmaxB=EmaxB, E0A=E0A, E0B=E0B, nA=nA, nB=nB, EC50A=EC50A, EC50B=EC50B, extendInt="yes")$root
    return(res)
}

Loewe <- function(a,b, modA, modB){
    resLoewe(a, b, EmaxA=coefficients(modA)["Emax"], EmaxB=coefficients(modB)["Emax"], E0A=coefficients(modA)["E0"], E0B=coefficients(modB)["E0"], nA=coefficients(modA)["n"], nB=coefficients(modB)["n"], EC50A=coefficients(modA)["EC50"], EC50B=coefficients(modB)["EC50"])
}

                                   
tryLoewe <- function(a,b, modA, modB){
    t <- try(Loewe(a, b, modA, modB))
    if (class(t) == "try-error") res <- NA
    else res <- t
    return(res)
}

## fct pour trouver couples de A,B qui donnent une réponse attendue
## utilisation de uniroot via predict(gam) - expected
findPairs <- function(y, d, val, rangeRoot){
    fun2min <- function(x, d, y, val){
        interp(d$PosB, d$PosA, d$value, duplicate="mean", xo=x, yo=y, linear=TRUE)$z - val
    }
    uniroot(fun2min, rangeRoot, d=d, y=y, val=val)$root
}

tryFindPairs <- function(y, d, val, rangeRoot){
    t <- try(findPairs(y, d, val, rangeRoot))
    if (class(t) == "try-error") res <- NA
    else res <- t
    return(res)
}


mySeq <- function(doses, n){
    d <- data.frame(from=doses[-length(doses)], to=doses[-1])
    return(unique(as.vector(apply(d, 1, function(x, n) seq(x[1], x[2], length.out=n), n))))
}



fun2Optim <- function(x, d, y, val){
    abs(interp(d$PosB, d$PosA, d$value, duplicate="mean", xo=y, yo=x, linear=TRUE)$z - val)
}


detectDec <- function(file){
    d <- readLines(file, n=3)
    npoint <- length(grep("\\.", d))
    ncomma <- length(grep(",", d))
    return(ifelse(npoint > ncomma, ".", ","))
}
