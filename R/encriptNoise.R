encriptNoise <-
function(noise, a, b, maxorder, levels, EPS, noisefile)
{
  lvls <- levels
  save(a, b, maxorder, lvls, noise, EPS, file=noisefile)
}
