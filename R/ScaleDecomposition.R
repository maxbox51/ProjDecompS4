
require(Matrix)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# ScaleDecomposition:
#
# An application area for diagonal scalings. It is _not_ standard to separate
# the magnitude out from the rest of a diagonal scaling, however for purposes
# of decomposition of a matrix into different scaling factors, and especially
# for data normalization, it is a useful thing to do. Note that many, if not
# most, data normalization methods do not only involve scalings; it wouldn't
# be appropriate to call this class (something)Normalization.
#
# Likewise, most matrix decompositions are not limited to scaling.
#-------------------------------------------------------------------------------
.ScaleDecomposition <-
setClass ("ScaleDecomposition",
          slots = c(canonicalForm = "dMatrix",
                    magnitude     = "numeric",
                    scaling       = "DiagonalScaling")
)

#-------------------------------------------------------------------------------
# Slot accessor methods
#
setGeneric("getCanonicalForm", signature="obj",
           function(obj) standardGeneric("getCanonicalForm"))
setMethod("getCanonicalForm", "ScaleDecomposition",
          function(obj){ obj@canonicalForm })

setGeneric("getMagnitude", signature="obj",
           function(obj) standardGeneric("getMagnitude"))
setMethod("getMagnitude", "ScaleDecomposition",
          function(obj){ obj@magnitude })

setGeneric("getScaling", signature="obj",
           function(obj) standardGeneric("getScaling"))
setMethod("getScaling", "ScaleDecomposition",
          function(obj){ obj@scaling })

#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------
# inverse():  A scale decomposition has a well-defined inverse, which relates
#             a matrix in a canonically scaled form to a particular
#             non-canonically scaled equivalent matrix (instead of vice-versa).
#
setGeneric("inverse")
setMethod("inverse", "ScaleDecomposition",
function(obj){
    .ScalingDecomposition(
        canonicalForm = upscale(obj@scale * obj@canonicalForm, obj@scaling),
        magnitude     = 1/obj@magnitude,
        scaling       = inverse(obj@scaling))
})

#-------------------------------------------------------------------------------
# rms(obj): The root-mean-square of all elements in a numeric object.
#
setGeneric("rms", signature="obj",
           function(obj) standardGeneric("rms"))
setMethod("rms", "numeric",
function(obj) {
    sqrt(mean(as.vector(obj) * as.vector(obj)))
})

# end of ScaleDecomposition.R
