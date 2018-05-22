density_Rmask <-
function (moments, a, b, n = 512)
{
  span <- 25
  fYseq <- yseq<- seq(a, b, length.out = n)
  for(i in 1:n) {
    jseq <- i + (-span):span
    if(i <= span) jseq <- 1:i
    if(i > n-span) jseq <- i:n
    fYseq[i] <- length(jseq)/(2*span+1) * mean(fY_Rmask(yseq[jseq], moments, a, b))
  }
  scale=sum(fYseq)
  negative <- (1:n)[fYseq <= 0]
  fYseq[negative] <- 0
  output <- list(x = yseq, y = fYseq, call=sys.call())
  class(output) <- "density"
  return(output)
}
