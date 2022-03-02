group_iterator <- function(Say,
                           Sgy,
                           Syy,
                           LT,
                           LB,
                           s2,
                           g,
                           p){
  # Checks if there are elements in the Symmetric Network
  if(!prod(!(g == 1)[-p])){
    # Obtains the Relevant Values
    Sxy <- c(Say[p, -p][g[-p] == 1], Sgy[, p])
    L   <- c(LT[p, -p][g[-p] == 1], LB[, p])
    nT  <- sum((g == 1)[-p])
  } else {
    Sxy <- Sgy[, p]
    L   <- LB[, p]
    nT  <- 0
  }
  
  # Computes the Log-Odds
  bh <- Sxy / (Syy + 1 / L)
  lo <- -(1 / 2) * (log(L) + log(Syy + 1 / L)) + bh^2 * (Syy + 1 / L) / (2 * s2)
  
  # Computes the Probability
  od <- exp(sum(lo))
  if(is.infinite(od)){
    pr <- 1
  } else {
    pr <- od / (1 + od)
  }
  
  # Updates g
  g[p] <- rbinom(n = 1, size = 1, prob = pr)
  
  # Updates Theta and B
  if(g[p] == 1){
    b      <- rnorm(n    = length(L),
                    mean = 0,
                    sd   = sqrt(s2 / (Syy + 1 / L))) + bh
  } else {
    b      <- rep(0, length(L))
  }
  
  return(list(bh = bh,
              lo = lo,
              pr = pr,
              g  = g,
              b  = b))
}