library(RUnit)
library(ProjectiveDecomposition)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test.defaultConstructor()
   test.constructorWithSmallData()

} # runTests
#----------------------------------------------------------------------------------------------------
test.defaultConstructor <- function()
{
   printf("--- test.defaultConstructor")
   pd <- ProjectiveDecomposition()
   checkTrue("ProjectiveDecomposition" %in% is(pd))

} # test.defaultConstructor
#----------------------------------------------------------------------------------------------------
test.constructorWithSmallData <- function()
{
   printf("--- test.constructorWithSmallData")
   set.seed(17)

   mtx <- Matrix(data=runif(100, min=-10, max=10), nrow=10, byrow=TRUE)
   pd <- ProjectiveDecomposition(mtx)

   hVec <- getHorizontalVector(pd)
   expected <- c(1.490591, -9.837341, 6.294315, 8.884366, 7.757625, 3.681117, 2.045356, 6.062666, -6.786408, 4.796362)
   checkEqualsNumeric(hVec, expected, tol=1e-3)

   vVec <- getVerticalVector(pd)
   expected <-  c(1.143827e-03, 1.316429e-01, 7.449221e+00, 1.452078e+01, 4.772257e-03, 1.747349e+01, 4.737491e+00,
                   2.339341e+03, 2.195992e-03, 2.044964e+03)
   checkEqualsNumeric(vVec, expected, tol=1e-3)

} # test.constructor
#----------------------------------------------------------------------------------------------------

