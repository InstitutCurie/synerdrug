
## get/set doses slot
setGeneric("doses", function(object){
    standardGeneric("doses")
})
setGeneric("doses<-", function(object, value){
    standardGeneric("doses<-")
})

## get/set data slot
setGeneric("expData", function(object){
    standardGeneric("expData")
})
setGeneric("expData<-", function(object, value){
    standardGeneric("expData<-")
})

## get/set data mean slot
setGeneric("meanData", function(object){
    standardGeneric("meanData")
})
setGeneric("meanData<-", function(object, value){
    standardGeneric("meanData<-")
})


setGeneric("computeDataMean", function(object){
    standardGeneric("computeDataMean")
})

setGeneric("computeBliss", function(object){
    standardGeneric("computeBliss")
})

setGeneric("computeHSA", function(object){
    standardGeneric("computeHSA")
})

setGeneric("computeLoewe", function(object){
    standardGeneric("computeLoewe")
    })


setGeneric("computeChou", function(object){
    standardGeneric("computeChou")
})

setGeneric("HSA", function(object){
    standardGeneric("HSA")
})
setGeneric("Bliss", function(object){
    standardGeneric("Bliss")
})

setGeneric("Loewe", function(object){
    standardGeneric("Loewe")
})

setGeneric("Chou", function(object){
    standardGeneric("Chou")
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


setGeneric("drugNames", function(object){
     standardGeneric("drugNames")
    })
setGeneric("drugNames<-", function(object, value){
    standardGeneric("drugNames<-")
})

setGeneric("content", function(object){
    standardGeneric("content")
    })
setGeneric("content<-", function(object, value){
    standardGeneric("content<-")
})

setGeneric("typeHill", function(object){
    standardGeneric("typeHill")
})
setGeneric("typeHill<-", function(object, value){
    standardGeneric("typeHill<-")
})

