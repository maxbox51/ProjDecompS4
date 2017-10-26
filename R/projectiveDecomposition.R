
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

scaleHadamardSquare <- function(obj, ...) {
    H   <- obj
    H@x <- H@x*H@x # Hadamard square
    getScaling(doublyStochasticDecomposition(H, ...))
}

setMethod("projectiveDecomposition", "dgTMatrix",
function(obj, ...) {
    mag <- rms(obj)
    hsc <- scaleHadamardSquare(obj/mag, ...)
    sc  <- .DiagonalScaling(colFactor = sqrt(getColFactor(hsc)),
                            rowFactor = sqrt(getRowFactor(hsc)))
    .ScaleDecomposition(canonicalForm = downscale(sc, obj/mag),
                        magnitude     = mag,
                        scaling       = sc)
})

# end of projectiveDecomposition.R
