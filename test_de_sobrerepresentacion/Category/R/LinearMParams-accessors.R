
setMethod("geneIds", "LinearMParams", function(object, ...) object@geneIds)
setReplaceMethod("geneIds", "LinearMParams", function(object, value) {
    object@geneIds <- value
    object
})


setMethod("testDirection", "LinearMParams", function(r) r@testDirection)
setReplaceMethod("testDirection", "LinearMParams", function(r, value) {
    r@testDirection <- value
    r
})

setMethod("universeGeneIds", "LinearMParams", function(r) r@universeGeneIds)
setReplaceMethod("universeGeneIds", "LinearMParams", function(r, value) {
    r@universeGeneIds <- value
    r
})

setMethod("categoryName", "LinearMParams", function(r) r@categoryName)
setReplaceMethod("categoryName", "LinearMParams", function(r, value) {
    r@categoryName <- value
    r
})

setMethod("annotation", "LinearMParams", function(object) object@annotation)
setReplaceMethod("annotation", c("LinearMParams", "character"),
                 function(object, value) {
                   object@annotation <- value
                   object@datPkg <- DatPkgFactory(value)
                   object
                 })

setMethod("pvalueCutoff", "LinearMParams", function(r) r@pvalueCutoff)
setReplaceMethod("pvalueCutoff", "LinearMParams", function(r, value) {
    r@pvalueCutoff <- value
    r
})


setMethod("conditional", "LinearMParams", function(r) r@conditional)

setReplaceMethod("conditional", c("LinearMParams", "logical"),
                 function(r, value) {
                     if (is.na(value))
                         stop("value must be TRUE or FALSE")
                     r@conditional <- value
                     r
                 })

