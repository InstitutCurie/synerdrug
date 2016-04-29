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

makeDrugSyn <- function(data, doses, content = c("Death", "Survival")){
    drugNames <- names(doses)
    content <- match.arg(content, c("Death", "Survival"))
    data[, "Pos1"] <- match(data[, drugNames[1]], sort(doses[[1]]))
    data[, "Pos2"] <- match(data[, drugNames[2]], sort(doses[[2]]))
    object <- new("DrugSyn", data = data, doses = doses, drugNames = names(doses))
    content(object) <- content
    drugNames(object) <- drugNames
    object <- makeRespInd(object)
    object <- computeDataMean(object)
    object <- computeHSA(object)
    object <- computeBliss(object)
    object <- computeLoewe(object)
    object <- computeChou(object)
    return(object)
}

setMethod("makeRespInd", signature(object = "DrugSyn"), function(object){
    data <- expData(object)
    drugs <- drugNames(object)
    dataA <- data[data[drugs[2]] == 0, ]
    dataA[, drugs[2]] <- NULL
    dataB <- data[data[drugs[1]] == 0, ]
    dataB[, drugs[1]] <- NULL
    colnames(dataA)[3] <- colnames(dataB)[3] <- "conc"
    respA <- estimateHill(dataA)
    respB <- estimateHill(dataB)
    resp <- list(respA, respB)
    names(resp) <- drugs
    respInd(object) <- resp
    validObject(object)
    return(object)
})




setMethod("computeDataMean", signature(object = "DrugSyn"), function(object){
    data <- expData(object)
    drugs <- drugNames(object)
    listToAg <- setNames(list(data[, drugs[1]], data[, drugs[2]]), drugs)
    dataMean <- aggregate(data$value, by = listToAg, mean)
    dataMean[, "Pos1"] <- match(dataMean[, drugs[1]], sort(unique(data[, drugs[1]])))
    dataMean[, "Pos2"] <- match(dataMean[, drugs[2]], sort(unique(data[, drugs[2]])))
    ## predict values for individual responses
    respInd <- respInd(object)
    if (is.null(respInd)) stop("Please compute individual responses")
    dataMean[, paste0("pred", drugs[1])] <- predict(respInd[[drugs[1]]], list(conc = dataMean[, drugs[1]]))
    dataMean[, paste0("pred", drugs[2])] <- predict(respInd[[drugs[2]]], list(conc = dataMean[, drugs[2]]))
    meanData(object) <- dataMean
    validObject(object)
    return(object)
})

setMethod("computeBliss", signature(object = "DrugSyn"), function(object){
    data <- meanData(object)
    if (is.null(data)) {
        stop("Please compute meanData before computing Bliss")
    }
    drugs <- drugNames(object)
    pred1 <- data[, paste0("pred", drugs[1])]
    pred2 <- data[, paste0("pred", drugs[2])]
    data$Bliss <- (pred1 + pred2 - pred1 * pred2) / data$x
    meanData(object) <- data
    validObject(object)
    return(object)
})

setMethod("computeHSA", signature(object = "DrugSyn"), function(object){
    data <- meanData(object)
    content <- content(object)
    pminmax <- ifelse(content == "Death", pmax, pmin)
    drugs <- drugNames(object)
    if (is.null(data)) {
        stop("Please compute meanData before computing Bliss")
    }
    data$HSA <- pminmax(data[, paste0("pred", drugs[1])], data[, paste0("pred", drugs[2])]) / data$x
    meanData(object) <- data
    validObject(object)
    return(object)
})


setMethod("computeLoewe", signature(object = "DrugSyn"), function(object){
    data <- meanData(object)
    mods <- respInd(object)
    content <- content(object)
    drugs <- drugNames(object)
    if (content == "Death") pminmax <- pmax
    else pminmax <- pmin
    if (is.null(data)) {
        stop("Please compute meanData")
    }
    if (is.null(mods)) {
        stop("Please compute individual responses")
    }
    data$Loewe <- mapply(tryLoewe, a = data[, drugs[1]], b = data[, drugs[2]], MoreArgs = list(modA = mods[[1]], modB = mods[[2]]))
    data$Loewe <- ifelse(is.na(data$Loewe), pminmax(data[, paste0("pred", drugs[1])], data[, paste0("pred", drugs[2])]), data$Loewe)
    data$LoeweExcess <- data$x - data$Loewe
    meanData(object) <- data
    validObject(object)
    return(object)
    })


setMethod("computeChou", signature(object = "DrugSyn"), function(object){
    data <- meanData(object)
    mods <- respInd(object)
    drugs <- drugNames(object)
    if (is.null(data)) {
        stop("Please compute meanData")
    }
    if (is.null(mods)) {
        stop("Please compute individual responses")
    }
    dose1 <- data[, drugs[1]]
    dose2 <- data[, drugs[2]]
    data$Chou <- dose1 / EDpred(dose1, mods[[1]]) + dose2 / EDpred(dose2, mods[[2]])
    meanData(object) <- data
    validObject(object)
    return(object)
})

## get Loewe, HSA and Bliss values
setMethod("HSA", signature(object = "DrugSyn"), function(object){
    extractData(object, "HSA")
})
setMethod("Bliss", signature(object = "DrugSyn"), function(object){
    extractData(object, "Bliss")
})
setMethod("Loewe", signature(object = "DrugSyn"), function(object){
    extractData(object, c("Loewe", "LoeweExcess"))
})
setMethod("Chou", signature(object = "DrugSyn"), function(object){
    extractData(object, "Chou")
})

extractData <- function(object, what){
    data <- meanData(object)
    drugs <- drugNames(object)
    d <- data[, c(drugs, paste0("pred", drugs), "x", what)]
    return(d)
}

setMethod("content", signature(object = "DrugSyn"), function(object){
    object@content
    })
setReplaceMethod("content", signature(object = "DrugSyn", value = "character"), function(object, value = c("Death", "Survival")){
    value <- match.arg(value, c("Death", "Survival"))
    object@content <- value
    validObject(object)
    return(object)
})


#############################################################
setMethod("plot", signature(x = "DrugSyn"), function(x, what = c("heatmap", "parallel", "ind", "surface"), ...){
    what <- match.arg(what, c("heatmap", "parallel", "ind", "surface"))
    if (what == "surface") {
        plotSurface(x, ...)
        return(invisible())
    }
    if (what == "heatmap") g <- plotHeatmap(x, ...)
    if (what == "parallel") g <- parPlot(x, ...)
    if (what == "ind") g <- respPlot(x, ...)
    return(g)
    })
