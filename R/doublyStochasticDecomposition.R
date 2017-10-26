
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
#-------------------------------------------------------------------------------
# Interfaces for different calling signatures.
#
setGeneric("doublyStochasticDecomposition",signature="obj",
function(obj, ...) standardGeneric("doublyStochasticDecomposition"))

setMethod("doublyStochasticDecomposition", "matrix",
function(obj, ...) {
    if (is.numeric(obj)) {
        doublyStochasticDecomposition(as(obj,"dgTMatrix"), ...)
    } else {
        stop("")
    }
})

setMethod("doublyStochasticDecomposition", "dMatrix",
function(obj, ...) {
    doublyStochasticDecomposition(as(obj,"dgTMatrix"), ...)
})

#-------------------------------------------------------------------------------
# Internal support functions
#
###
# Reuse the parts of an existing scaling that fit the size of
# the matrix, as a starting point to reduce the number of iterations
# required to converge.
###
initial.scaling <-
function(rows, cols, start, verbose) {
    if (is.null(start)) {
       rf <- rep(1.0, rows)
       cf <- rep(1.0, rows)
    } else {
       if (rows == length(start@rowFactor)) {
           rf <- start@rowFactor
           if (cols == length(start@colFactor)) {
               cf <- start@colFactor
               if (verbose) {
                   cat("Reusing row and column factors\n")
               }
           } else {
               cf <- rep(1.0, rows)
               if (verbose) {
                   cat("Reusing row factors\n")
               }
           }
       } else {
           rf <- rep(1.0, rows)
           if (cols == length(start@colFactor)) {
               cf <- start@colFactor
               cat("Reusing column factors\n")
           } else {
               cf <- rep(1.0, rows)
           }
       }
    }
    if (min(rf) <= 0.0) {
        stop("Row factors must be strictly positive")
    }
    if (min(cf) <= 0.0) {
        stop("Column factors must be strictly positive")
    }
    .DiagonalScaling(rowFactor = rf, colFactor = cf)
}

# Allows e.g. cv(x, na.rm=TRUE)
cv <- function(x,...) {
    sqrt(var(x,...))/mean(x,...)
}


#-------------------------------------------------------------------------------
# Full implementation using the dgTMatrix class.
#
setMethod("doublyStochasticDecomposition", "dgTMatrix",
function(obj,
         start     = NULL,
         precision = 1e-6,
         max.iter  = 1000,
         verbose   = FALSE,
         digits    = 3) {
    if (min(obj@x) >= 0.0) {
        m   <- dim(obj)[1]
        n   <- dim(obj)[2]
        sc  <- initial.scaling(m,n,start,verbose)
        rf  <- getRowFactor(sc)
        cf  <- getColFactor(sc)
        mag <- sum(obj@x) / (m*n)
        a   <- obj@x / mag
        can <- obj
        can@x <- a/(rf[1+can@i]*cf[1+can@j])
        for (iter in c(1:max.iter)) {
            rfa <- rowSums(can) / n
            cfa <- colSums(can) / m
            coeff.v <- max(cv(rfa[rfa > 0.0]), cv(cfa[cfa > 0.0]))
            if (verbose) {
                cat(sprintf("%d %s %s\n",
                            iter,signif(coeff.v,digits), date()))
            }
            if (coeff.v < precision) { break }
            rf    <- rf * rfa
            can@x <- a/(rf[1+can@i]*cf[1+can@j])
            cf    <- cf * colSums(can) / m
            can@x <- a/(rf[1+can@i]*cf[1+can@j])
        }
        .ScaleDecomposition(canonicalForm = can,
                            magnitude = mag,
                            scaling = 
            .DiagonalScaling(rowFactor=rf,colFactor=cf))
    } else {
        stop("requires a non-negative matrix")
    }
})

# end of ProjectiveDecomposition.R
