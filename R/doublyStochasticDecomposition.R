
require(Matrix)

#-------------------------------------------------------------------------------
# Doubly-Stochastic Decomposition (L1-based scaling)
#    The matrix in canonical form satisfies the constraints:
#        sum(row) = length(row)
#        sum(col) = length(col)
#    for all rows and columns. Note that this is a different overall
#    scaling for square matrices than is usually meant by doubly-stochastic
#    matrices, however it is (1) compatible with a rational definition of
#    canonical scale for non-square matrices, (2) unbiased toward either the
#    row-stochastic or column-stochastic version of doubly-stochastic, and
#    easy to convert to either one (just multiply the magnitude
#    by 1/length(row) for row-stochastic or 1/length(col) for column-
#    stochastic; for square matrices, the result is the same).
#
# This stub implementation gets the magnitude part right, but doesn't do
# the relative row- and column-scalings.
#
setGeneric("doublyStochasticDecomposition",signature="obj",
function(obj, ...) standardGeneric("doublyStochasticDecomposition"))

setMethod("doublyStochasticDecomposition", "matrix",
function(obj, ...) {
    if (is.numeric(obj)) {
        doublyStochasticDecomposition(as(obj,"dgTMatrix"))
    } else {
        stop("")
    }
})

setMethod("doublyStochasticDecomposition", "dMatrix",
function(obj, ...) {
    doublyStochasticDecomposition(as(obj,"dgTMatrix"))
})

setMethod("doublyStochasticDecomposition", "dgTMatrix",
function(obj, ...) {
    mag <- sum(obj)/prod(dim(obj))
    .ScalingDecomposition(canonicalForm = obj / mag,
                          magnitude     = mag,
                          scaling       = 
        .DiagonalScaling(colFactor = rep(1.0,dim(obj)[2]),
                         rowFactor = rep(1.0,dim(obj)[1])))
})

# end of ProjectiveDecomposition.R
