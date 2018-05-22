rmulti <-
function(n, mean, sd, p)
{
  x <- rnorm(n)
  u <- sample(1:length(mean), n, prob=p, replace=T)
  for(i in 1:length(mean)) x[u==i]=mean[i]+sd[i]*x[u==i]
  return(x)
}
