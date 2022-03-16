### Simple log Odds Function Profile

logOddPro <- function(s   = 1,
                      Sxx = 1,
                      Sxy = 1,
                      t   = 1,
                      ver = TRUE){
  if(ver){
    bh <- Sxy / (Sxx + 1 / s)
    lo <- - log(s) / 2 - log(Sxx + 1 / s) / 2 + bh^2 * (Sxx + 1 / s) / (2 * t)
  } else {
    lo <- - log(Sxx * s + 1) / 2 + Sxy^2 / (2 * t * (Sxx + 1 / s))
  }
  ### Returns the function
  return(lo)
}