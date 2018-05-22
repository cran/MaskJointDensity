lambda_Rmask <-
function (muX, k) 
{
  i <- 0:floor(k/2)
  val <- (2 * k + 1)/2 * sum((-1)^i/2^k * exp(lgamma(2 *
                                                       k - 2 * i + 1) - lgamma(i + 1) - lgamma(k - i + 1) -
                                                lgamma(k - 2 * i + 1)) * muX[k - 2 * i + 1])
  return(val)
}
