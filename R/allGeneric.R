## get/set doses slot ----

#' doses
#'
#' Extract doses values
#' @param object a DrugSyn object
#' @export
setGeneric("doses",
           function(object)
               standardGeneric("doses")
)

#' Set doses values
#' @param object a DrugSyn object
#' @param value list
#'
#' @export
setGeneric("doses<-", function(object, value)
    standardGeneric("doses<-")
)

## get/set data slot ----
#' Get experimental data
#' @param object a \code{DrugSyn} object
setGeneric("expData", function(object)
    standardGeneric("expData")
)

#' Set experimental data
#' @param object a \code{DrugSyn} object
#' @param value a \code{data.frame}
setGeneric("expData<-", function(object, value)
    standardGeneric("expData<-")
)

## get/set data mean slot
setGeneric("meanData", function(object)
    standardGeneric("meanData")
)
setGeneric("meanData<-", function(object, value)
    standardGeneric("meanData<-")
)


setGeneric("computeDataMean", function(object)
    standardGeneric("computeDataMean")
)

setGeneric("computeBliss", function(object)
    standardGeneric("computeBliss")
)

setGeneric("computeHSA", function(object)
    standardGeneric("computeHSA")
)

setGeneric("computeLoewe", function(object)
    standardGeneric("computeLoewe")
    )


setGeneric("computeChou", function(object)
    standardGeneric("computeChou")
)

setGeneric("HSA", function(object)
    standardGeneric("HSA")
)
setGeneric("Bliss", function(object)
    standardGeneric("Bliss")
)

setGeneric("Loewe", function(object)
    standardGeneric("Loewe")
)

setGeneric("Chou", function(object)
    standardGeneric("Chou")
    )

## individual response slot ----
#' Get individual responses parameters
#' @param object a \code{DrugSyn} object
setGeneric("respInd", function(object)
    standardGeneric("respInd")
)
#' Set individual responses
#' @param object a \code{DrugSyn} object
#' @param value list
setGeneric("respInd<-", function(object, value)
    standardGeneric("respInd<-")
)

setGeneric("makeRespInd", function(object)
  standardGeneric("makeRespInd")
)

## drug names ----
#' Get drug names
#' @param object a \code{DrugSyn} object
setGeneric("drugNames", function(object)
     standardGeneric("drugNames")
    )
#' Set drug names
#' @param object a \code{DrugSyn} object
#' @param value character vector
setGeneric("drugNames<-", function(object, value)
    standardGeneric("drugNames<-")
)

#' Get content type
#' @param object a \code{DrugSyn} object
setGeneric("content", function(object)
    standardGeneric("content")
)
#' Set content type
#' @param object a \code{DrugSyn} object
#' @param value "death" or "survival"
setGeneric("content<-", function(object, value)
    standardGeneric("content<-")
)


setGeneric("typeHill", function(object)
    standardGeneric("typeHill")
)
setGeneric("typeHill<-", function(object, value)
    standardGeneric("typeHill<-")
)

