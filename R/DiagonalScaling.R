
require(Matrix)
# Two draft S4 classes for Max
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Diagonal scaling is a standard term for the operation method "upscale"
# performs; I'm adopting it for the associated data object.
# "downscale" is "upscaling" by the inverse scaling, without having to
# explicitly represent the inverse as an object; "upscale" is to "downscale"
# as multiplication is to division.
#
.DiagonalScaling <-
setClass ("DiagonalScaling",
          slots = c(colFactor="numeric",
                    rowFactor="numeric")
)

#-------------------------------------------------------------------------------
# Slot accessor methods
#
setGeneric("getColFactor", function(obj) standardGeneric("getColFactor"))
setMethod("getColFactor", "DiagonalScaling", function(obj){ obj@colFactor })

setGeneric("getRowFactor", function(obj) standardGeneric("getRowFactor"))
setMethod("getRowFactor", "DiagonalScaling", function(obj){ obj@rowFactor })

#-------------------------------------------------------------------------------
# show()
#
setGeneric("show")
setMethod("show", "DiagonalScaling",
    function (object){
        cat(sprintf("Diagonal Scaling object of dimension (%d x %d)\n",
                    length(object@rowFactor), length(object@colFactor)))
    }
)

#-------------------------------------------------------------------------------
# summary()
#
setGeneric("summary")
setMethod("summary", "DiagonalScaling",
    function(object,
             digits = max(3L, getOption("digits") - 3L),
             display = 5){
        m <- length(object@rowFactor)
        n <- length(object@colFactor)
        cat("(",paste(m,"x",n),")\ncols: ")
        if (m <= display) {
            cat(paste(signif(getRowFactor(object),digits)))
        } else {
            cat(paste(signif(getRowFactor(object)[1:display],digits)))
            cat(" ...")
        }
        cat("\ncols: ")
        if (n <= display) {
            cat(paste(signif(getColFactor(object),digits)))
        } else {
            cat(paste(signif(getColFactor(object)[1:display],digits)))
            cat("... ")
        }
        cat("\n")
    }
)

#-------------------------------------------------------------------------------
# t(): the transpose of an operator with rows and columns
#
setGeneric("t")
setMethod("t", "DiagonalScaling",
function(x) {
  .DiagonalScaling(colFactor = x@rowFactor, rowFactor = x@colFactor)
})

#-------------------------------------------------------------------------------
# inverse(): the inverse of a matrix operator
#
setGeneric("inverse",signature="obj", function(obj) standardGeneric("inverse"))

setMethod("inverse", "DiagonalScaling",
function(obj) {
  .DiagonalScaling(colFactor = 1.0 / obj@colFactor,
                   rowFactor = 1.0 / obj@rowFactor)
})


#-------------------------------------------------------------------------------
# upscale(), downscale():
#
# Methods defining how a DiagonalScaling operates on a matrix.
# These return the matrix class which they are given, whether or not
# it is an S4 class.
#
# upscale():  "Multiply" by the diagonal scaling.
#
setGeneric("upscale",      signature=c("obj","matrix"),
function(obj,matrix) standardGeneric("upscale"))

setMethod("upscale", c("DiagonalScaling","matrix"),
function(obj, matrix){
    matrix * obj@rowFactor[row(matrix)] * obj@colFactor[col(matrix)]
})

setMethod("upscale", c("DiagonalScaling","dMatrix"),
function(obj, matrix){
    matrix@x * obj@rowFactor[row(matrix)] * obj@colFactor[col(matrix)]
})

#
# downscale(): "Divide" by the diagonal scaling.
#
setGeneric("downscale",    signature=c("obj","matrix"),
           function(obj,matrix) {stop("not meaningful for arbitrary objects")})

setMethod("downscale", c("DiagonalScaling","matrix"),
function(obj,matrix) {
    matrix / (obj@rowFactor[row(matrix)] * obj@colFactor[col(matrix)])
})

setMethod("downscale", c("DiagonalScaling","dMatrix"),
function(obj,matrix) {
    matrix / (obj@rowFactor[row(matrix)] * obj@colFactor[col(matrix)])
})

# end of class DiagonalScaling
