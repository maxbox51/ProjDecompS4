# a draft S4 class for Max
#------------------------------------------------------------------------------------------------------------------------
.ProjectiveDecomposition <- setClass ("ProjectiveDecomposition",
                                      representation = representation(
                                         matrix="Matrix",
                                         horizontalVector="numeric",
                                         verticalVector="numeric",
                                         scale="numeric"
                                         ),
                            )

setGeneric("show", signature="obj", function(obj) standardGeneric("show"))
setGeneric("getHorizontalVector", signature="obj", function(obj) standardGeneric("getHorizontalVector"))
setGeneric("getVerticalVector", signature="obj", function(obj) standardGeneric("getVerticalVector"))
#------------------------------------------------------------------------------------------------------------------------
ProjectiveDecomposition <- function(matrix=Matrix(),
                                    horizontalVector=NA_numeric_,
                                    verticalVector=NA_numeric_,
                                    scale=NA_numeric_)
{
   horizontalVector <- sample(as.numeric(matrix), size=ncol(matrix))        # dumb example
   verticalVector <- exp(sample(as.numeric(matrix), size=nrow(matrix)))     # another

   scale <- sum(exp(matrix))  # and another

   obj <- .ProjectiveDecomposition(matrix=matrix, horizontalVector=horizontalVector,
                                   verticalVector=verticalVector,
                                   scale=scale)

} # ProjectiveDecomposition
#------------------------------------------------------------------------------------------------------------------------
setMethod("show", "ProjectiveDecomposition",

      function(obj){
          msg = sprintf("Projective Decomposition object, base matrix of dimension %d %d",
                        nrow(obj@matrix), ncol(obj@matrix))
          cat(msg, "\n", sep="")
          cat(sprintf("horiontalVector of length %d", length(obj@horizontalVector)), "\n", sep="")
          })

#------------------------------------------------------------------------------------------------------------------------
setMethod("getHorizontalVector", "ProjectiveDecomposition",

      function(obj){
          obj@horizontalVector
          })

#------------------------------------------------------------------------------------------------------------------------
setMethod("getVerticalVector", "ProjectiveDecomposition",

      function(obj){
          obj@verticalVector
          })

#------------------------------------------------------------------------------------------------------------------------
