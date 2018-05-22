moments <-
function (y, order = 8) 
{
  val <- 1
  if (order > 0) for (k in 1:order)
    val <- c(val, mean(y^k))
  return(val)
}
