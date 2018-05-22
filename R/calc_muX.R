calc_muX <-
function (muY, a, b) 
{
  muX <- 0 * muY
  muX[1] <- 1
  for (j in 1:(length(muY) - 1)) {
    k <- 0:j
    muX[j + 1] <- sum(choose(j, k) * 2^k * muY[k + 1] *
                        (-1)^(j - k) * (a + b)^(j - k))/(b - a)^j
  }
  return(muX)
}
