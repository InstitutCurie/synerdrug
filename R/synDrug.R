## class definition ----

#' DrugSyn class
#'
#' S4 class to store drug combination data
#'
#' @slot data data.frame
#' @slot doses list
#' @slot respInd list
#' @slot drugNames vector
#' @slot dataMean data.frame
#' @slot content character
#' @slot typeHill numeric
#' @slot range character
#'
#' @return
#' @export
#' @exportClass DrugSyn
#'
#' @examples
setClass("DrugSyn", representation(data = "data.frame", doses = "list", respInd = "list", drugNames = "vector", dataMean = "data.frame", content = "character", typeHill = "numeric", range = "character"))

setValidity("DrugSyn", function(object){

    return(TRUE)
})

## show ----
setMethod("show", signature(object = "DrugSyn"),
          function(object){
              drugs <- drugNames(object)
              doses <- doses(object)
              cat("DrugSyn object\n")
              cat("drugs:", drugs)
              #print(dim(expData(object)))
              cat("\ndoses: \n\t")
              cat(drugNames(object)[1], paste0(doses[[1]], collapse = " "))
              cat("\n\t")
              cat(drugNames(object)[2], paste0(doses[[2]], collapse = " "))
          })


###########################################
## setter / getter ----

## doses ----
#' @describeIn DrugSyn get doses values
#' @param object a DrugSyn
#' @param value replacement value
#' @export
setMethod("doses", signature(object = "DrugSyn"), function(object){
    object@doses
})
#' @describeIn DrugSyn set doses values
#' @export
setReplaceMethod("doses", signature(object = "DrugSyn", value = "list"), function(object, value){
    object@doses <- value
    validObject(object)
    return(object)
})


#' @describeIn DrugSyn extract experimental data
#' @export
setMethod("expData", signature(object = "DrugSyn"), function(object){
    object@data
})
#' @describeIn DrugSyn set experimental data
#' @export
setReplaceMethod("expData", signature(object = "DrugSyn", value = "data.frame"), function(object, value){
    object@data <- value
    return(object)
})

## meanData getter/setter ----
setMethod("meanData", signature(object = "DrugSyn"), function(object){
    object@dataMean
})
setReplaceMethod("meanData", signature(object = "DrugSyn", value = "data.frame"), function(object, value){
    object@dataMean <- value
    return(object)
})

## individual responses getter/setter ----
#' @describeIn DrugSyn extract individual responses
#' @export
setMethod("respInd", signature(object = "DrugSyn"), function(object){
    object@respInd
})
#' @describeIn DrugSyn set individual responses
#' @export
setReplaceMethod("respInd", signature(object = "DrugSyn", value = "list"), function(object, value){
    object@respInd <- value
    validObject(object)
    return(object)
})

## drug names getter / setter ----
#' @describeIn DrugSyn extract drug names
#' @export
setMethod("drugNames", signature(object = "DrugSyn"), function(object){
    object@drugNames
})
#' @describeIn DrugSyn set drug names
#' @export
setReplaceMethod("drugNames", signature(object = "DrugSyn", value = "vector"), function(object, value){
    oldNames <- drugNames(object)
    object@drugNames <- value
    names(object@doses) <- value
    names(object@respInd) <- value
    coln <- colnames(expData(object))
    names(value) <- oldNames
    colnames(expData(object))[coln %in% oldNames] <- value[coln[coln %in% oldNames]]
    coln <- colnames(meanData(object))
    colnames(meanData(object))[coln %in% oldNames] <- value
    names(value) <- paste0("pred", oldNames)
    colnames(meanData(object))[coln %in% paste0("pred", oldNames)] <- paste0("pred", value[coln[coln %in% paste0("pred", oldNames)]])
    validObject(object)
    return(object)
})

## hill getter/setter ----
setMethod("typeHill", signature(object = "DrugSyn"), function(object){
    object@typeHill
    })
setReplaceMethod("typeHill", signature(object = "DrugSyn", value = "numeric"), function(object, value){
    object@typeHill <- value
    validObject(object)
    return(object)
})

## makeDrugSyn ----
#' Create DrugSyn object
#'
#' @param data data.frame of experimental values
#' @param doses list with doses of the two drugs
#' @param content character, type of the experiment
#' @param typeHill numeric, number of parameters for the Hill function
#' @param range "Fraction" or "Percentage", defines the range of the values
#'
#' @return a \code{DrugSyn} object
#' @export
#'
#' @examples
#'
#' ## data.frame of experimental values: d
#' object <- makeDrugSyn(d, doses = list(A = sort(unique(d$A)), B = sort(unique(d$B))),
#'  content = "Death")
makeDrugSyn <- function(data, doses, content = c("Death", "Survival"), typeHill = c(4, 3), range = c("Fraction", "Percentage")){
    drugNames <- names(doses)
    content <- match.arg(content, c("Death", "Survival"))
    typeHill <- checkTypeHill(typeHill)
    range <- match.arg(range, c("Fraction", "Percentage"))
    data <- data[order(data[, drugNames[2]], data[, drugNames[1]]), ]
    data[, "Pos1"] <- match(data[, drugNames[1]], sort(doses[[1]]))
    data[, "Pos2"] <- match(data[, drugNames[2]], sort(doses[[2]]))
    if (range == "Percentage") data[, "value"] <- data[, "value"] / 100
    object <- new("DrugSyn", data = data, doses = doses, drugNames = names(doses), typeHill = typeHill, range = range)
    content(object) <- content
    object <- makeRespInd(object)
    object <- computeDataMean(object)
    object <- computeHSA(object)
    object <- computeBliss(object)
    object <- computeLoewe(object)
    object <- computeChou(object)
    drugNames(object) <- drugNames
    return(object)
}

## makeRespInd ----
setMethod("makeRespInd", signature(object = "DrugSyn"), function(object){
    data <- expData(object)
    drugs <- drugNames(object)
    dataA <- data[data[drugs[2]] == 0, ]
    dataA[, drugs[2]] <- NULL
    dataB <- data[data[drugs[1]] == 0, ]
    dataB[, drugs[1]] <- NULL
    colnames(dataA)[3] <- colnames(dataB)[3] <- "conc"
    nbParam <- typeHill(object)
    respA <- estimateHill(dataA, nbParam)
    respB <- estimateHill(dataB, nbParam)
    resp <- list(respA, respB)
    names(resp) <- drugs
    respInd(object) <- resp
    validObject(object)
    return(object)
})

## computeDataMean ----
setMethod("computeDataMean", signature(object = "DrugSyn"), function(object){
    data <- expData(object)
    drugs <- drugNames(object)
    listToAg <- setNames(list(data[, drugs[1]], data[, drugs[2]]), drugs)
    dataMean <- aggregate(data$value, by = listToAg, mean)
    dataMean[, "Pos1"] <- match(dataMean[, drugs[1]], sort(unique(data[, drugs[1]])))
    dataMean[, "Pos2"] <- match(dataMean[, drugs[2]], sort(unique(data[, drugs[2]])))
    colnames(dataMean)[3] <- "value"
    ## predict values for individual responses
    respInd <- respInd(object)
    if (is.null(respInd)) stop("Please compute individual responses")
    dataMean[, paste0("pred", drugs[1])] <- mod2ll4(respInd[[drugs[1]]])(dataMean[, drugs[1]])
    dataMean[, paste0("pred", drugs[2])] <- mod2ll4(respInd[[drugs[2]]])(dataMean[, drugs[2]])
    meanData(object) <- dataMean
    validObject(object)
    return(object)
})

## computeBliss ----
setMethod("computeBliss", signature(object = "DrugSyn"), function(object){
    data <- meanData(object)
    if (is.null(data)) {
        stop("Please compute meanData before computing Bliss")
    }
    drugs <- drugNames(object)
    pred1 <- data[, paste0("pred", drugs[1])]
    pred2 <- data[, paste0("pred", drugs[2])]
    data$Bliss <- (pred1 + pred2 - pred1 * pred2) / data[, "value"]
    meanData(object) <- data
    validObject(object)
    return(object)
})

## computeHSA ----
setMethod("computeHSA", signature(object = "DrugSyn"), function(object){
    data <- meanData(object)
    content <- content(object)
    pminmax <- ifelse(content == "Death", pmax, pmin)
    drugs <- drugNames(object)
    if (is.null(data)) {
        stop("Please compute meanData before computing Bliss")
    }
    data$HSA <- pminmax(data[, paste0("pred", drugs[1])], data[, paste0("pred", drugs[2])]) / data[, "value"]
    meanData(object) <- data
    validObject(object)
    return(object)
})

## computeLoewe ----
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
    data$LoeweExcess <- data[, "value"] - data[, "Loewe"]
    meanData(object) <- data
    validObject(object)
    return(object)
    })

## computeChou ----
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
    value <- data[, "value"]
    data$Chou <- dose1 / EDpred(value, mods[[1]]) + dose2 / EDpred(value, mods[[2]])
    meanData(object) <- data
    validObject(object)
    return(object)
})

## synergy indexes getters ----
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

## extractData ----
extractData <- function(object, what){
    data <- meanData(object)
    drugs <- drugNames(object)
    d <- data[, c(drugs, paste0("pred", drugs), "value", what)]
    return(d)
}

## content ----
#' @describeIn DrugSyn get content type
#' @export
setMethod("content", signature(object = "DrugSyn"), function(object){
    object@content
    })
#' @describeIn DrugSyn set content type
#' @export
setReplaceMethod("content", signature(object = "DrugSyn", value = "character"), function(object, value = c("Death", "Survival")){
    value <- match.arg(value, c("Death", "Survival"))
    object@content <- value
    validObject(object)
    return(object)
})

## plot -----
#' @param x DrugSyn object
#' @param type type of plot
#' @param ... arguments given to plots functions
#'
#' @return
#' @aliases plot,DrugSyn-method
#' @describeIn DrugSyn plot DrugSyn object
#' @export
#'
#' @examples
setMethod("plot", signature(x = "DrugSyn"), function(x, type = c("heatmap", "parallel", "ind", "surface"), ...){
    type <- match.arg(type, c("heatmap", "parallel", "ind", "surface"))
    if (type == "surface") {
        plotSurface(x, ...)
        return(invisible())
    }
    if (type == "heatmap") g <- plotHeatmap(x, ...)
    if (type == "parallel") g <- parPlot(x, ...)
    if (type == "ind") g <- respPlot(x, ...)
    return(g)
    })
