d <- read.table("../longData.txt", header = TRUE)


d2 <- read.table("~/projets/SynerDrug/simdata2.txt", sep = "\t")
doseB <- d2[1,]
d2 <- melt(d2[-1,], id.vars = "V1")
d2$B <- unlist(doseB[as.character(d2$variable)])
d2$A <- d2$V1
d2$rep <- 1
d2 <- d2[,colnames(d)]

object <- new("DrugSyn", data = d2, doses = list(A = unique(d2$A), B = unique(d2$B)))
object <- makeRespInd(object)



object <- makeDrugSyn(d2, doses = list(A = unique(d2$A), B = unique(d2$B)))
object <- computeDataMean(object)
object <- computeBliss(object)
object <- computeHSA(object)



## réponse ind
respPlot(object, "A")
respPlot(object, "C")

## parallel plot
parPlot(object, 2)

## HSA
## Bliss
## Loewe
## CI
## Isobologram Loewe + Steel
## median effect



## Loewe
dm <- meanData(object)
mods <- respInd(object)
Loewe(dm$A[22], dm$B[22], mods[[1]], mods[[2]])

mapply(tryLoewe, a=dm$A, b=dm$B, MoreArgs = list(modA=mods[[1]], modB=mods[[2]]))

#####################################################################################
## reflexion
## plot : 1 methode à décliner avec un type = c(heatmap, parPlot, respInd) ?
## faire le calcul des indices dès la création de l'objet ?









## interpolation
x <- 1:3
y <- 1:3
z <- matrix(c(10,25,12,43, 32,45,56,78, 100), ncol = 3)
#z <- matrix(1:4, ncol = 2)
xo <- c(1.1, 1.3, 2.5, 2.7, 2.9)
yo <- c(1.5, 2.7, 1.8, 1.8, 2.9)

xo <- seq(1, 3, by=0.1)
yo <- seq(1,3, by=0.1)

interpol <- function(x, y, z, xo, yo){
    ## z matrice de données
    ## x, y coordonnées de la grille
    ## from fields::interp.surface
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


pos <- expand.grid(xo, yo)

my <- interpol(x, y, z, pos$Var1, pos$Var2)
my <- cbind(pos, my)


zz<-melt(z)
ak<-melt(akima::interp(zz$Var1,zz$Var2,zz$value, xo, yo, linear = TRUE)$z)
ak$Var1 <- xo[ak$Var1]
ak$Var2 <- yo[ak$Var2]

a <- merge(ak, my)




