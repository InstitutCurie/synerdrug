require(reshape2)


d <- read.table("../longData.txt", header = TRUE)


d2 <- read.table("~/projets/SynerDrug/simdata2.txt", sep = "\t")
doseB <- d2[1,]
d2 <- melt(d2[-1,], id.vars = "V1")
d2$B <- unlist(doseB[as.character(d2$variable)])
d2$A <- d2$V1
d2$rep <- 1
d2 <- d2[,colnames(d)]

object <- new("DrugSyn", data = d2, doses = list(A = unique(d2$A), B = unique(d2$B)), content = "Death", drugNames=c("A", "B"))
object <- makeRespInd(object)


colnames(d2) <- c("value", "rep", "AA", "BB")

object <- makeDrugSyn(d, doses = list(A = unique(d$A), B = unique(d$B)), content = "Death")
#object <- computeDataMean(object)
#object <- computeBliss(object)
#object <- computeHSA(object)



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
my <- cbind(my, my2)

zz <- melt(z)
ak <- melt(akima::interp(zz$Var1,zz$Var2,zz$value, xo, yo, linear = TRUE)$z)
ak$Var1 <- xo[ak$Var1]
ak$Var2 <- yo[ak$Var2]

a <- merge(ak, my)

## tests sur vrai données
xo <- yo <- mySeq(sort(unique(d2$A)),1)
pos <- expand.grid(xo, yo)
ak <- melt(akima::interp(d2$A, d2$B, d2$value, xo, yo, linear = TRUE)$z)
ak$Var1 <- xo[ak$Var1]
ak$Var2 <- yo[ak$Var2]
my <- interpol(sort(unique(d2$A)), sort(unique(d2$B)), matrix(d2$value, ncol=6), pos$Var1, pos$Var2)
my <- cbind(pos, my)
a <- merge(ak, my)

ggplotly(ggplot(a, aes(x=value, y=my)) + geom_point())
ggplotly(ggplot(a, aes(x=Var1, y=Var2, fill=value-my)) + geom_tile()+scale_x_log10()+scale_y_log10())

a$pos1 <- match(a$Var1, xo)
a$pos2 <- match(a$Var2, yo)

ggplotly(ggplot(a, aes(x=pos1, y=pos2, fill=value-my)) + geom_tile())






D <- read.table("~/Desktop/Shiny DB75-Erlot 5J.txt", sep="\t")
NAlines <-  which(apply(D, 1, function(x) all(is.na(x))))
if (length(NAlines) == 0) NAlines <- nrow(D) + 1
nB <- ncol(D) - 1
nA <- NAlines[1] - 2
doseB <- unlist(D[1,-1])
D <- D[-c(1, NAlines, NAlines+1), ]
doseA <- D[,1]
nRep <- nrow(D)/nA
D <- melt(D, id.vars="V1")
D$concB <- doseB[D$variable]
D$variable <- NULL
D$rep <- rep(rep(1:nRep, each = nA), nB)
D$concA <- D$V1
D$V1 <- NULL
D$PosA <- match(D$concA, sort(unique(D$concA)))
D$PosB <- match(D$concB, sort(unique(D$concB)))
D$A <- D$PosA
D$B <- D$PosB
colnames(D) <- c("value", "concB", "rep", "concA", "PosA", "PosB", "A", "B")
D$valueOrig <- D$value
D$value <- D$value/100


D$A <- D$concA
D$B <- D$concB
#D$value <- 1 - D$value
D <- D[, c("value", "rep", "A", "B")]
object <- makeDrugSyn(D, doses = list(A = unique(D$A), B = unique(D$B)))

modA <- respInd(object)[[1]]
modB <- respInd(object)[[2]]


d <- subset(D, B == 0)
d$conc <- d$A
mod1 <- nlsLM(value ~ E0 + (Emax - E0) / (1 + (EC50 / conc) ^ n),
             start = list(Emax = 1, E0 = 0, EC50 = median(d[,"conc"]), n = 1),
             data = d)
mod2 <- nlsLM(1-value ~ E0 + (Emax - E0) / (1 + (EC50 / conc) ^ n),
              start = list(Emax = 1, E0 = 0, EC50 = median(d[,"conc"]), n = 1),
              data = d)

pred1 <- predict(mod1, list(conc=seq(0,20, length.out = 200)))
pred2 <- predict(mod2, list(conc=seq(0,20, length.out = 200)))

EmaxA = coefficients(modA)["Emax"]; EmaxB = coefficients(modB)["Emax"]
E0A = coefficients(modA)["E0"]; E0B = coefficients(modB)["E0"]
nA = coefficients(modA)["n"]; nB = coefficients(modB)["n"]
EC50A = coefficients(modA)["EC50"]; EC50B = coefficients(modB)["EC50"]

res <- uniroot(fun2min, c(0,1),
               a = a, b = b, EmaxA = EmaxA, EmaxB = EmaxB, E0A = E0A,
               E0B = E0B, nA = nA, nB = nB, EC50A = EC50A, EC50B = EC50B,
               extendInt = "yes")$root

ggplotly(plot(object, what="heatmap", type="x", addVal=TRUE))


