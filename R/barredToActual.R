barredToActual <-
function(vectorL, boundaryVec, barred)
{
  size <- nrow(barred)
  dim <- ncol(barred)
  
  actual <- matrix(nrow = size, ncol = dim)
  
  for(i in 1:size)
  {
    for(j in 1:dim)
    {
      actual[i,j] <- boundaryVec[2*j - 1] +
        (boundaryVec[2*j] - boundaryVec[2*j - 1])/(vectorL[j] - 1)*(barred[i,j]-1)
    }
  }
  return(actual)
}
