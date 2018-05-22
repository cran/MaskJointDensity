fY_Rmask <-
function(y, muY, a, b) 
{
  x <- (2 * y - (a + b))/(b - a)
  muX <- calc_muX(muY, a, b)
  SUM <- 0 * y
  for (k in 0:(length(muY) - 1)) SUM <- SUM + lambda_Rmask(muX, k) * P_Rmask(x, k)
  SUM[(y < a) || (y > b)] <- 0
  val <- SUM * 2/(b - a)
  val[is.na(val)] <- 0
  return(val)
}
