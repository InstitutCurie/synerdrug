## ---- echo = FALSE-------------------------------------------------------
require(knitr)
opts_chunk$set(echo = TRUE, fig.width = 5, fig.height = 4, prompt = TRUE, warning = FALSE)

## ------------------------------------------------------------------------
require(synerdrug)
data(exDrugSyn)

## ------------------------------------------------------------------------
exDrugSyn
drugNames(exDrugSyn)
doses(exDrugSyn)

## ------------------------------------------------------------------------
plot(exDrugSyn)

## ------------------------------------------------------------------------
plot(exDrugSyn, type = "surface")
plot(exDrugSyn, type = "surface", theta = -45, phi = 45)

## ---- fig.height = 3, fig.width = 6--------------------------------------
plot(exDrugSyn, type = "parallel", ref = "drugA")
plot(exDrugSyn, type = "parallel", ref = "drugB")


## ---- fig.height = 3, fig.width = 6--------------------------------------
plot(exDrugSyn, type = "ind", drug = "drugA")
plot(exDrugSyn, type = "ind", drug = "drugB")

## ------------------------------------------------------------------------
respInd(exDrugSyn)

## ------------------------------------------------------------------------

plot(exDrugSyn, type = "heatmap", what = "Bliss")
plot(exDrugSyn, type = "surface", what = "LoeweExcess")


## ------------------------------------------------------------------------
isobologram(exDrugSyn, effect = 0.5)
isobologram(exDrugSyn, effect = 0.5, mode = "SP")

