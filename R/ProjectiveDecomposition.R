
require(Matrix)
# Two draft S4 classes for Max
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
# Diagonal scaling is a standard term for what method "upscale" does; I'm adopting it for the associated data object.
# "downscale" is upscaling by the inverse, without having to explicitly represent the inverse; if "upscale" is a
# multiplication, "downscale" is the corresponding division.
#
.DiagonalScaling <-
setClass ("DiagonalScaling",
          slots = c(colFactor="numeric",
                    rowFactor="numeric")
)

#------------------------------------------------------------------------------------------------------------------------
# Slot accessor methods
#
setGeneric("getColFactor", function(obj) standardGeneric("getColFactor"))
setMethod("getColFactor", "DiagonalScaling", function(obj){ obj@colFactor })

setGeneric("getRowFactor", function(obj) standardGeneric("getRowFactor"))
setMethod("getRowFactor", "DiagonalScaling", function(obj){ obj@rowFactor })

#------------------------------------------------------------------------------------------------------------------------
# show()
#
setGeneric("show")
setMethod("show", "DiagonalScaling",
    function (object){
        cat(sprintf("Diagonal Scaling object of dimension (%d x %d)\n",
                    length(object@rowFactor), length(object@colFactor)))
    }
)

#------------------------------------------------------------------------------------------------------------------------
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

#------------------------------------------------------------------------------------------------------------------------
# t(): the transpose of an operator with rows and columns
#
setGeneric("t")
setMethod("t", "DiagonalScaling",
function(x) {
  .DiagonalScaling(colFactor = x@rowFactor, rowFactor = x@colFactor)
})

#------------------------------------------------------------------------------------------------------------------------
# inverse(): the inverse of a matrix operator
#
setGeneric("inverse",signature="obj",
function(obj) {
    stop("Objects are not invertible in general")
})

setMethod("inverse", "DiagonalScaling",
function(obj) {
  .DiagonalScaling(colFactor = 1.0 / obj@colFactor, rowFactor = 1.0 / obj@rowFactor)
})


#------------------------------------------------------------------------------------------------------------------------
# upscale(), downscale():
#
# Methods defining how a DiagonalScaling operates on a matrix. These are written to return the
# matrix class which they are given, whether it is an S4 class or not.
#
# upscale():  "Multiply" by the diagonal scaling.
#
setGeneric("upscale",      signature=c("obj","matrix"),
           function(obj,matrix) {stop("not meaningful for arbitrary objects")})

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

#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
# ScaleDecomposition:
#
# An application area for diagonal scalings. It is _not_ standard to separate the magnitude out from the
# rest of a diagonal scaling, however for purposes of decomposition of a matrix into different scaling factors,
# and especially for data normalization, it is a useful thing to do. Note that many, if not most, data normalization
# methods do not only involve scalings, so it wouldn't be appropriate to call this class (something)Normalization.
#
# Likewise, most matrix decompositions are not limited to scaling.
#------------------------------------------------------------------------------------------------------------------------
.ScaleDecomposition <-
setClass ("ScaleDecomposition",
          slots = c(canonicalForm = "dMatrix",
                    magnitude     = "numeric",
                    scaling       = "DiagonalScaling")
)

#------------------------------------------------------------------------------------------------------------------------
# Slot accessor methods
#
setGeneric("getCanonicalForm", signature="obj", function(obj) standardGeneric("getCanonicalForm"))
setMethod("getCanonicalForm", "ScaleDecomposition", function(obj){ obj@canonicalForm })

setGeneric("getMagnitude", signature="obj", function(obj) standardGeneric("getMagnitude"))
setMethod("getMagnitude", "ScaleDecomposition", function(obj){ obj@magnitude })

setGeneric("getScaling", signature="obj", function(obj) standardGeneric("getScaling"))
setMethod("getScaling", "ScaleDecomposition", function(obj){ obj@scaling })

#------------------------------------------------------------------------------------------------------------------------
# t(): transposition of a scale decomposition
#
setGeneric("t")
setMethod("t", "ScaleDecomposition",
function(x){
    .ScalingDecomposition(
        canonicalForm = t(x@canonicalForm),
        magnitude     = x@magnitude,
        scaling       = t(x@scaling))
})

#------------------------------------------------------------------------------------------------------------------------
# inverse():  A scale decomposition has a well-defined inverse, which relates a matrix in a canonically scaled
#             form to a particular non-canonically scaled equivalent matrix instead of vice-versa.
#
setGeneric("inverse")
setMethod("inverse", "ScaleDecomposition",
function(obj){
    .ScalingDecomposition(
        canonicalForm = upscale(obj@scale * obj@canonicalForm, obj@scaling),
        magnitude     = 1/obj@magnitude,
        scaling       = inverse(obj@scaling))
})

#------------------------------------------------------------------------------------------------------------------------
# rms(obj): The root-mean-square of all elements in a numeric object.
#
setGeneric("rms", signature="obj", function(obj) standardGeneric("rms"))
setMethod("rms", "numeric", function(obj) { sqrt(mean(as.vector(obj) * as.vector(obj))) })

#------------------------------------------------------------------------------------------------------------------------
# Projective Decomposition (an L2-based scaling)
# This gets the magnitude part right, but doesn't do the row and column scalings
# Satisfies the criterion: rms(row) = rms(col) = 1.0.
#
setGeneric("projectiveDecomposition",signature="obj",
function(obj, ...) standardGeneric("projectiveDecomposition"))

setMethod("projectiveDecomposition", "dgTMatrix",
function(obj, ...) {
    mag <- rms(obj)
    .ScalingDecomposition(canonicalForm = obj / mag,
                          magnitude     = mag,
                          scaling       =
        .DiagonalScaling(colFactor = rep(1.0,dim(obj)[2]),
                         rowFactor = rep(1.0,dim(obj)[1])))
}) # projectiveDecomposition

#------------------------------------------------------------------------------------------------------------------------
# Doubly-Stochastic Decomposition (L1-based scaling, defined for some non-negative matrices)
# A version of Sinkhorn's algorithm computes the following decomposition.
# This gets the magnitude part right, but doesn't do the row and column scalings
# Unscaled matrix satisfies the criteria, for all rows and columns:
#    sum(row) = length(row)
#    sum(col) = length(col)
#
setGeneric("doublyStochasticDecomposition",signature="obj",
function(obj, ...) standardGeneric("doublyStochasticDecomposition"))

setMethod("doublyStochasticDecomposition", "dgTMatrix",
function(obj, ...) { # Must be a real-valued matrix, broadest virtual class is dMatrix
    mag <- sum(obj)/prod(dim(obj))
    .ScalingDecomposition(canonicalForm = obj / mag,
                          magnitude     = mag,
                          scaling       = 
        .DiagonalScaling(colFactor = rep(1.0,dim(obj)[2]),
                         rowFactor = rep(1.0,dim(obj)[1])))
}) # doublyStochasticDecomposition

# end of ProjectiveDecomposition.R
