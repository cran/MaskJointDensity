mask <-
function (vectorToBeMasked, noisefile, noise = createNoise(length(vectorToBeMasked)), lowerBoundAsGivenByProvider=min(vectorToBeMasked), upperBoundAsGivenByProvider=max(vectorToBeMasked),
                  maxorder = 100, EPS = 1e-06) 
  
{
  
  n <- length(vectorToBeMasked)
  
  if(length(noise) != n) stop("'vectorToBeMasked' and 'noise' lengths differ")
  
  if(sum(noise < 0) > 0) stop("Negative noise not allowed")
  
  lvls <- NULL
  
  FACTOR <- is.factor(vectorToBeMasked)
  
  if (FACTOR) {
    
    k <- length(levels(vectorToBeMasked))
    
    lvls <- levels(vectorToBeMasked)
    
    y1 <- rep(NA, length(vectorToBeMasked))
    
    for (i in 1:k) y1[vectorToBeMasked == lvls[i]] <- i
    
    vectorToBeMasked <- y1
    
  }
  
  if (sum(noise < 0) > 0) 
    
    stop("Negative noise not allowed")
  
  ndens <- density(noise)
  
  prob <- ndens$y/sum(ndens$y)
  
  ystar = vectorToBeMasked * noise
  
  C2 <- sample(noise, size = 10*n, replace = TRUE)
  
  if (FACTOR) {
    
    a <- 0
    
    b <- k + 1
    
  }
  
  else {
    
    #   a <- min(vectorToBeMasked) - runif(1, 0, sd(vectorToBeMasked))
    
    #  b <- max(vectorToBeMasked) + runif(1, 0, sd(vectorToBeMasked))
    
    a<- lowerBoundAsGivenByProvider
    b<- upperBoundAsGivenByProvider
    
    k <- 0
    
  }
  
  encriptNoise(C2, a, b, maxorder, lvls, EPS, noisefile)
  
  return(list(ystar = ystar, noisefile = noisefile))
  
}
