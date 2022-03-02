iboom_iterator <- function(Say,
                           Sgy,
                           Syy,
                           G,
                           A,
                           y,
                           Theta,
                           B,
                           t2T,
                           l2T,
                           vT,
                           xiT,
                           t2B,
                           l2B,
                           vB,
                           xiB,
                           s2,
                           g){
  # Problem dimensions
  n <- length(y)  # Number of Observations
  V <- dim(B)[1]  # Voxel Size
  P <- dim(B)[2]  # Number of ROI's
  Q <- sum(g)     # Number of Active Regions
  gp <- numeric(length = P)
  
  # Samples s2
  Q   <- sum(g == 1)
  bs2 <-       sum((array(data = Theta %x% y, dim = dim(A)) - A)^2) / 2
  bs2 <- bs2 + sum((array(data = B %x% y, dim = dim(G)) - G)^2)
  bs2 <- bs2 + sum(Theta[g == 1, g == 1]^2 / l2T[g == 1, g == 1]) / 2
  bs2 <- bs2 + sum(B[, g == 1]^2 / l2B[, g == 1])
  bs2 <- bs2 / 2
  as2 <- n * P * (P - 1) / 2
  as2 <- as2 + n * V * P
  as2 <- as2 + Q * (Q - 1) / 2
  as2 <- as2 + Q * V
  as2 <- as2 / 2
  s2  <- 1 / rgamma(n = 1, shape = as2, rate = bs2)
    
  # # Samples r
  # if(!is.numeric(r)){
  #   r <- rbeta(n      = 1,
  #              shape1 = sum(g) + ar,
  #              shape2 = P - sum(g) + ar)
  # }
  
  # Horseshoe Structure for B
  # Samples l2B
  for(i in 1:P){
    if(g[i] == 1){
      l2B[, i] <- 1 / rgamma(n     = V,
                             shape = 1,
                             rate  = 1 / vB[, i] + B[, i]^2 / (2 * s2 * t2B))
    } else {
      l2B[, i] <- 1 / rgamma(n     = V,
                             shape = 1 / 2,
                             rate  = 1 / vB[, i])
    }
  }

  # # Samples r2B
  # for(i in 1:P){
  #   if(g[i] > 0){
  #     r2B[i] <- 1 / rgamma(n     = 1,
  #                          shape = (V + 1) / 2,
  #                          rate  = 1 / wB[i] +
  #                            sum(B[, i]^2 / (2 * s2 * l2B[, i] * t2B)))
  #   } else {
  #     r2B[i] <- 1 / rgamma(n     = 1,
  #                          shape = 1 / 2,
  #                          rate  = 1 / wB[i])
  #   }
  #   r2B[i] <- 1
  # }
  # 
  # # Samples t2B
  # if(sum(g) > 0){
  #   t2B <- 1 / rgamma(n     = 1,
  #                     shape = (sum(g) * V + 1) / 2,
  #                     rate  = 1 / xiB +
  #                       sum(B^2 / (2 * s2 * t(t(l2B) * r2B))))
  # } else {
  #   t2B <- 1 / rgamma(n     = 1,
  #                     shape = 1 / 2,
  #                     rate  = 1 / xiB)
  # }
  
  # Samples t2B
  if(sum(g) > 0){
    t2B <- 1 / rgamma(n     = 1,
                      shape = (sum(g) * V + 1) / 2,
                      rate  = 1 / xiB +
                        sum(B^2 / (2 * s2 * l2B)))
  } else {
    t2B <- 1 / rgamma(n     = 1,
                      shape = 1 / 2,
                      rate  = 1 / xiB)
  }
  

  # Samples vB
  vB  <- 1 / rgamma(n     = P * V,
                    shape = 1,
                    rate  = 1 + 1 / c(l2B))
  vB <- matrix(data = vB, nrow = V, ncol = P)
  # 
  # # Samples w2B
  # wB  <- 1 / rgamma(n     = P,
  #                   shape = 1,
  #                   rate  = 1 + 1 / r2B)

  # Samples xiB
  xiB <- 1 / rgamma(n     = 1,
                    shape = 1,
                    rate  = 1 + 1 / t2B)

  # Horseshoe Structure for Theta
  # Samples l2T
  for(i in 2:P){
    for(j in 1:(i - 1)){
      if(g[i] * g[j] == 1){
        l2T[i, j] <- 1 / rgamma(n     = 1,
                                shape = 1,
                                rate  = 1 / vT[i, j] + Theta[i, j]^2 / (2 * s2 * t2T))
        l2T[j, i] <- l2T[i, j]
      } else {
        l2T[i, j] <- 1 / rgamma(n     = 1,
                                shape = 1 / 2,
                                rate  = 1 / vT[i, j])
        l2T[j, i] <- l2T[i, j]
      }
    }
  }

  # Samples t2T
  if(sum(g) > 0){
    t2T <- 1 / rgamma(n     = 1,
                      shape = (sum(g) * (sum(g) - 1) / 2 + 1) / 2,
                      rate  = 1 / xiT +
                        sum(Theta^2 / (2 * s2 * l2T)) / 2)
  } else {
    t2T <- 1 / rgamma(n     = 1,
                      shape = 1 / 2,
                      rate  = 1 / xiT)
  }

  # Samples vT
  out <- 1 / rgamma(n     = P * (P - 1) / 2,
                    shape = 1,
                    rate  = 1 + 1 / l2T[lower.tri(l2T)])
  vT[upper.tri(vT, diag = TRUE)] <- 0
  vT[lower.tri(vT)]              <- out
  vT                             <- vT + t(vT)

  # Samples xiT
  xiT <- 1 / rgamma(n     = 1,
                    shape = 1,
                    rate  = 1 + 1 / t2T)
  
  # Samples g, Theta and B
  for(p in 1:P){
    res <- group_iterator(Say = Say, 
                          Sgy = Sgy,
                          Syy = Syy,
                          LT  = t2T * l2T,
                          LB  = t2B * l2B,
                          s2  = s2,
                          g   = g,
                          p   = p)
    gp[p] <- res$pr
    # Updates g
    g <- res$g
    # Updates Theta
    Q <- sum(g[-p])
    if(Q > 0){
      Theta[-p, p][g[-p] == 1] <- res$b[1:Q]
      Theta[p, -p][g[-p] == 1] <- Theta[-p, p][g[-p] == 1]
    }
    # Updates B
    B[, p] <- res$b[(Q + 1):(Q + V)]
  }
  
  # Returns Values
  return(list(Theta = Theta,
              B     = B,
              g     = g,
              s2    = s2,
              as2   = as2,
              bs2   = bs2,
              l2T   = l2T,
              t2T   = t2T,
              vT    = vT,
              xiT   = xiT,
              l2B   = l2B,
              t2B   = t2B,
              vB    = vB,
              xiB   = xiB,
              gp    = gp))
}