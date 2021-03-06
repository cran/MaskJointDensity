qkdeSorted <- function (p, fhat) 
{
  if (any(p > 1) | any(p < 0)) 
    stop("p must be <= 1 and >= 0")
  cumul.prob <- ks::pkde(q = fhat$eval.points, fhat = fhat) # we make sure it is sorted
  if (is.unsorted(cumul.prob)) {  
    if(max(abs(sort(cumul.prob)-cumul.prob)) > 10^(-9)) {
      print(which(abs(sort(cumul.prob)-cumul.prob) > 10^(-9) ))
      stop("not sorted and sorting would make a large difference")
    } else {
      warning("cumul.prob was unsorted, but the difference is not large so it was sorted")
      cumul.prob <- sort(cumul.prob)
    } 
  }
  ind <- findInterval(x = p, vec = cumul.prob)
  quant <- rep(0, length(ind))
  for (j in 1:length(ind)) {
    i <- ind[j]
    if (i == 0) 
      quant[j] <- fhat$eval.points[1]
    else if (i >= length(fhat$eval.points)) 
      quant[j] <- fhat$eval.points[length(fhat$eval.points)]
    else {
      quant1 <- fhat$eval.points[i]
      quant2 <- fhat$eval.points[i + 1]
      prob1 <- cumul.prob[i]
      prob2 <- cumul.prob[i + 1]
      alpha <- (p[j] - prob2)/(prob1 - prob2)
      quant[j] <- quant1 * alpha + quant2 * (1 - alpha)
    }
  }
  return(quant)
}