setClass("DrugSyn", representation(data = "data.frame", doses = "list", respInd = "list"))


setValidity("DrugSyn", function(object){

    return(TRUE)
})


setMethod("show", signature(object = "DrugSyn"),
          function(object){
            cat(dim(getData(object)))
          })

setMethod("doses", signature(object = "DrugSyn"), function(object){
    object@doses
})
setReplaceMethod("doses", signature(object = "DrugSyn", value = "list"), function(object, value){
    object@doses <- value
    validObject(object)
    return(object)
})


setMethod("getData", signature(object = "DrugSyn"), function(object){
    object@data
})
setMethod("setData", signature(object = "DrugSyn", value = "data.frame"), function(object, value){
    object@data <- data
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


setMethod("makeRespInd", signature(object = "DrugSyn"), function(object){
    data <- getData(object)
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

makeDrugSyn <- function(data, doses){
    object <- new("DrugSyn", data = data, doses = doses)
    validObject(object)
    return(object)
}
