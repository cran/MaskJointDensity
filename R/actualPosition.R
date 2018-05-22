actualPosition <-
function(vectorL, prob, boundaryVec, size = 1)
{
  len <- length(vectorL)
  n <- c()
  
  for(i in 1:len)
  {
    n[i] <- vectorL[i]
  }
  
  maxSize = 1
  for(i in 1:len)
  {
    maxSize <- maxSize*n[i]
  }
  w <- c(0:(maxSize-1))
  k <- sample(w, size = size, replace = TRUE, prob = prob)
  
  barredPoints <- matrix(nrow = size, ncol = len)
  
  for(i in 1:size)
  {
    maxSize <- 1
    for(l in 1:len)
    {
      maxSize <- maxSize*n[l]
    }
    for(j in 1:(len-1))
    {
      maxSize <- maxSize/n[len+1-j]
      if(k[i] > maxSize)
      {
        barredPoints[i,(len+1-j)] <- floor(k[i]/maxSize) + 1
        k[i] <- k[i]%%maxSize +1  # This line is modified.
      }
      else
      {
        barredPoints[i,(len+1-j)] <- 1
      }
    }
    barredPoints[i,1] <- k[i]+1  # This line is modified.
  }
  
  return(barredToActual(vectorL, boundaryVec, barredPoints))
  
}
