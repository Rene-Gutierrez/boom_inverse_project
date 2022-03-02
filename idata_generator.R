idata_generator <- function(PP = 20,
                            VV = 10,
                            pB = 0.5,
                            pT = 0.5,
                            cB = 1,
                            cT = 1,
                            s2 = 0.5,
                            n  = 25){
  # Computes Active Regions
  gT              <- rep(0, PP)
  gT[1:(PP * pT)] <- 1
  # Computes Active Voxels
  gB              <- rep(0, VV)
  gB[1:(VV * pB)] <- 1
  
  # Samples Theta
  Theta                   <- (gT %*% t(gT)) * cT
  Theta[upper.tri(Theta)] <- 0
  diag(Theta)             <- 0
  Theta                   <- Theta + t(Theta)
  
  # Sample B
  B <- gB %*% t(gT) * cB
  
  # Samples y
  y <- rnorm(n    = n,
             mean = 0,
             sd   = 1)
  # Standarizes
  y <- (y - mean(y)) / sd(y)
  
  # Samples A
  A <- array(data = Theta %x% y,
             dim  = c(n, PP, PP)) + array(data = rnorm(n    = PP * PP,
                                                       mean = 0,
                                                       sd   = sqrt(s2)),
                                          dim  = c(n, PP, PP))
  for(i in 1:n){
    # A[i,,] <- Theta * y[i]
    A[i,,][lower.tri(A[i,,])] <- 0
    A[i,,] <- A[i,,] + t(A[i,,]) 
  }
  
  # Samples B
  G <- array(data = B %x% y,
             dim  = c(n, VV, PP)) + array(data = rnorm(n    = VV * PP,
                                                       mean = 0,
                                                       sd   = sqrt(s2)),
                                          dim  = c(n, VV, PP))
  
  # Centers
  mA <- apply(X = A, MARGIN = c(2, 3), FUN = mean)
  mG <- apply(X = G, MARGIN = c(2, 3), FUN = mean)
  A  <- A - array(data = mA %x% rep(1, n), dim = c(n, PP, PP))
  G  <- G - array(data = mG %x% rep(1, n), dim = c(n, VV, PP))
  
  # Returns Values
  return(list(Theta = Theta,
              B     = B,
              gT    = gT,
              gB    = gB,
              A     = A,
              G     = G,
              y     = y,
              P     = PP,
              V     = VV,
              n     = n,
              s2    = s2))
  
}
