
## get/set doses slot
setGeneric("doses", function(object){
    standardGeneric("doses")
})
setGeneric("doses<-", function(object, value){
    standardGeneric("doses<-")
})

## get/set data slot
setGeneric("getData", function(object){
    standardGeneric("getData")
})
setGeneric("setData", function(object, value){
    standardGeneric("setData")
})

## get/set individual response slot
setGeneric("respInd", function(object){
    standardGeneric("respInd")
})
setGeneric("respInd<-", function(object, value){
    standardGeneric("respInd<-")
})

setGeneric("makeRespInd", function(object){
  standardGeneric("makeRespInd")
})
