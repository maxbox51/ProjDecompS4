
require(Matrix)

#-------------------------------------------------------------------------------
# Projective Decomposition (an L2-based scaling; data normalization)
#     Satisfies the criterion: rms(row) = rms(col) = 1.0.
#
# This stub implementation gets the magnitude part right, but doesn't do
# the row and column scalings.
#
setGeneric("projectiveDecomposition",signature="obj",
function(obj, ...) standardGeneric("projectiveDecomposition"))

setMethod("projectiveDecomposition", "matrix",
function(obj, ...) {
    if (is.numeric(obj)) {
        projectiveDecomposition(as(obj,"dgTMatrix"), ...)
    } else {
        stop("only applies to numeric matrices")
    }
})

setMethod("projectiveDecomposition", "dMatrix",
function(obj, ...) {
    projectiveDecomposition(as(obj,"dgTMatrix"), ...)
})

setMethod("projectiveDecomposition", "dgTMatrix",
function(obj, ...) {
    mag <- rms(obj)
    .ScalingDecomposition(canonicalForm = obj / mag,
                          magnitude     = mag,
                          scaling       =
        .DiagonalScaling(colFactor = rep(1.0,dim(obj)[2]),
                         rowFactor = rep(1.0,dim(obj)[1])))
})

# end of projectiveDecomposition.R
