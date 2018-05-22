createNoise <-
function (n, mean = rep(NA, 5), sd = rep(1, length(mean)),
                         
                         prob = rep(1, length(mean))/length(mean)) 
  
{
  
  k <- length(mean)
  
  if (is.na(sum(mean))) {
    
    best <- rnorm(k, 500, 200)
    
    MAX <- min(dist(best))
    
    for (i in 1:1000) {
      
      mean <- rnorm(k, 500, 200)
      
      if ((min(dist(mean)) > MAX) && (min(mean) > 50)) {
        
        MAX <- min(dist(mean))
        
        best <- mean
        
      }
      
    }
    
    mean <- best
    
  }
  
  C <- rnorm(n)
  
  grp <- sample(1:length(mean), n, prob = prob, replace = TRUE)
  
  for (i in 1:length(mean)) C[grp == i] = mean[i] + sd[i] * 
    
    C[grp == i]
  
  C[C < 0] <- -C[C < 0]
  
  return(C)
  
}
