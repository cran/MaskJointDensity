P_Rmask <-
function (x, k) 
{
  SUM <- 0 * x
  for (i in 0:floor(k/2)) SUM <- SUM + (-1)^i/2^k *
    exp(lgamma(2 * k - 2 * i + 1) - lgamma(i + 1) -
          lgamma(k - i + 1) - lgamma(k - 2 * i + 1)) * x^(k - 2 * i)
  return(SUM)
}
