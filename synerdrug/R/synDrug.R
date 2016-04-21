setClass("DrugSyn", representation(data = "data.frame", doses = "list", respInd = "list", drugNames = "vector", dataMean = "data.frame", content = "character"))

setValidity("DrugSyn", function(object){

    return(TRUE)
})

setMethod("show", signature(object = "DrugSyn"),
          function(object){
              print(dim(expData(object)))
              print(doses(object))
          })


## doses
setMethod("doses", signature(object = "DrugSyn"), function(object){
    object@doses
})
setReplaceMethod("doses", signature(object = "DrugSyn", value = "list"), function(object, value){
    object@doses <- value
    validObject(object)
    return(object)
})


setMethod("expData", signature(object = "DrugSyn"), function(object){
    object@data
})
setReplaceMethod("expData", signature(object = "DrugSyn", value = "data.frame"), function(object, value){
    object@data <- value
    return(object)
})

setMethod("meanData", signature(object = "DrugSyn"), function(object){
    object@dataMean
})
setReplaceMethod("meanData", signature(object = "DrugSyn", value = "data.frame"), function(object, value){
    object@dataMean <- value
    return(object)
})

setMethod("respInd", signature(object = "DrugSyn"), function(object){
    object@respInd
})
setReplaceMethod("respInd", signature(object = "DrugSyn", value = "list"), function(object, value){
    object@respInd <- value
    validObject(object)
    return(object)
})

setMethod("drugNames", signature(object = "DrugSyn"), function(object){
    object@drugNames
})
setReplaceMethod("drugNames", signature(object = "DrugSyn", value = "vector"), function(object, value){
    object@drugNames <- value
    validObject(object)
    return(object)
})

########################################################
## !!!! handle drug names properly
########################################################
setMethod("makeRespInd", signature(object = "DrugSyn"), function(object){
    data <- expData(object)
    dataA <- subset(data, B == 0)
    dataA$B <- NULL
    dataB <- subset(data, A == 0)
    dataB$A <- NULL
    colnames(dataA)[3] <- colnames(dataB)[3] <- "conc"
    respA <- estimateHill(dataA)
    respB <- estimateHill(dataB)
    respInd(object) <- list(A = respA, B = respB)
    validObject(object)
    return(object)
})

makeDrugSyn <- function(data, doses, content = c("Death", "Survival")){
    drugNames <- names(doses)
    content <- match.arg(content, c("Death", "Survival"))
    data[, "Pos1"] <- match(data[, drugNames[1]], sort(doses[[1]]))
    data[, "Pos2"] <- match(data[, drugNames[2]], sort(doses[[2]]))
    object <- new("DrugSyn", data = data, doses = doses, drugNames = names(doses))
    content(object) <- content
    object <- makeRespInd(object)
    object <- computeDataMean(object)
    object <- computeBliss(object)
    object <- computeLoewe(object)
    object <- computeChou(object)
    return(object)
}

setMethod("computeDataMean", signature(object = "DrugSyn"), function(object){
    data <- expData(object)
    drugNames <- drugNames(object)
    dataMean <- aggregate(data$value, by = list(A = data$A, B = data$B), mean)
    dataMean[, "Pos1"] <- match(dataMean[, "A"], sort(unique(data[, "A"])))
    dataMean[, "Pos2"] <- match(dataMean[, "B"], sort(unique(data[, "B"])))
    ## predict values for individual responses
    respInd <- respInd(object)
    if (is.null(respInd)) stop("Please compute individual responses")
    dataMean$predA <- predict(respInd$A, list(conc = dataMean$A))
    dataMean$predB <- predict(respInd$B, list(conc = dataMean$B))
    meanData(object) <- dataMean
    validObject(object)
    return(object)
})

setMethod("computeBliss", signature(object = "DrugSyn"), function(object){
    data <- meanData(object)
    if (is.null(data)) {
        stop("Please compute meanData before computing Bliss")
    }
    data$Bliss <- (data$predA + data$predB - data$predA * data$predB) / data$x
    meanData(object) <- data
    validObject(object)
    return(object)
})

setMethod("computeHSA", signature(object = "DrugSyn"), function(object){
    data <- meanData(object)
    if (is.null(data)) {
        stop("Please compute meanData before computing Bliss")
    }
    data$HSA <- pmax(data$predA, data$predB) / data$x
    meanData(object) <- data
    validObject(object)
    return(object)
})


setMethod("computeLoewe", signature(object = "DrugSyn"), function(object){
    data <- meanData(object)
    mods <- respInd(object)
    content <- content(object)
    if (content == "Death") pminmax <- pmax
    else pminmax <- pmin
    if (is.null(data)) {
        stop("Please compute meanData")
    }
    if (is.null(mods)) {
        stop("Please compute individual responses")
    }
    data$Loewe <- mapply(tryLoewe, a = data$A, b = data$B, MoreArgs = list(modA = mods[[1]], modB = mods[[2]]))
    data$Loewe <- ifelse(is.na(data$Loewe), pminmax(data$predA, data$predB), data$Loewe)
    data$LoeweExcess <- data$x - data$Loewe
    meanData(object) <- data
    validObject(object)
    return(object)
    })


setMethod("computeChou", signature(object = "DrugSyn"), function(object){
    data <- meanData(object)
    mods <- respInd(object)
    if (is.null(data)) {
        stop("Please compute meanData")
    }
    if (is.null(mods)) {
        stop("Please compute individual responses")
    }
    data$Chou <- data$A / EDpred(data$A, mods[[1]]) + data$B / EDpred(data$B, mods[[2]])
    meanData(object) <- data
    validObject(object)
    return(object)
})

## get HSA and Bliss values
setMethod("HSA", signature(object = "DrugSyn"), function(object){
    data <- meanData(object) ## what if not computed ?
    data[, c("A", "B", "predA", "predB", "x", "HSA")]
})
setMethod("Bliss", signature(object = "DrugSyn"), function(object){
    data <- meanData(object) ## what if not computed ?
    data[, c("A", "B", "predA", "predB", "x", "Bliss")]
})
setMethod("Loewe", signature(object = "DrugSyn"), function(object){
    data <- meanData(object)
    data[, c("A", "B", "predA", "predB", "x", "Loewe", "LoeweExcess")]
})


setMethod("content", signature(object = "DrugSyn"), function(object){
    object@content
    })
setReplaceMethod("content", signature(object = "DrugSyn", value = "character"), function(object, value = c("Death", "Survival")){
    value <- match.arg(value, c("Death", "Survival"))
    object@content <- value
    validObject(object)
    return(object)
})


#############################################################"
setMethod("plot", signature(x = "DrugSyn"), function(x){
    NULL
    })
