iboom_sampler <- function(y,
                          G,
                          A,
                          nmcmc  = 1000,
                          burnin = 0){
  # Computes Sufficient Statistics and Problem Dimensions
  n   <- dim(G)[1]
  V   <- dim(G)[2]
  P   <- dim(G)[3]
  Say <- apply(X = (A * y), MARGIN = c(2, 3), sum)
  Sgy <- apply(X = (G * y), MARGIN = c(2, 3), sum)
  Syy <- sum(y^2)
  
  # Sample Holders
  # Non-Burn-in
  sTheta <- array(data = NA, dim = c(burnin + nmcmc, P, P))
  sB     <- array(data = NA, dim = c(burnin + nmcmc, V, P))
  sg     <- matrix(data = NA, nrow = burnin + nmcmc, ncol = P)
  sgp    <- matrix(data = NA, nrow = burnin + nmcmc, ncol = P)
  ss2    <- numeric(length = burnin + nmcmc)
  sl2T   <- array(data = NA, dim = c(burnin + nmcmc, P, P))
  st2T   <- numeric(length = burnin + nmcmc)
  svT    <- array(data = NA, dim = c(burnin + nmcmc, P, P))
  sxiT   <- numeric(length = burnin + nmcmc)
  sl2B   <- array(data = NA, dim = c(burnin + nmcmc, V, P))
  st2B   <- numeric(length = burnin + nmcmc)
  svB    <- array(data = NA, dim = c(burnin + nmcmc, V, P))
  sxiB   <- numeric(length = burnin + nmcmc)
  
  # Initialization
  Theta <- matrix(data = 0, nrow = P, ncol = P)
  B     <- matrix(data = 0, nrow = V, ncol = P)
  t2T   <- 1
  l2T   <- matrix(data = 1, nrow = P, ncol = P)
  xiT   <- 1
  vT    <- matrix(data = 1, nrow = P, ncol = P)
  t2B   <- 1
  l2B   <- matrix(data = 1, nrow = V, ncol = P)
  xiB   <- 1
  vB    <- matrix(data = 1, nrow = V, ncol = P)
  s2    <- 1
  g     <- rep(1, P)
  
  # Progress Bar
  pb <- txtProgressBar(min     = 0,
                       max     = 1,
                       initial = 0,
                       style   = 3,
                       width   = 72)
  # First Run
  out <- iboom_iterator(Say   = Say,
                        Sgy   = Sgy,
                        Syy   = Syy,
                        A     = A,
                        G     = G,
                        y     = y,
                        Theta = Theta,
                        B     = B,
                        t2T   = t2T,
                        l2T   = l2T,
                        xiT   = xiT,
                        vT    = vT,
                        t2B   = t2B,
                        l2B   = l2B,
                        xiB   = xiB,
                        vB    = vB,
                        s2    = s2,
                        g     = g)
  # Sampling
  for(s in 1:(burnin + nmcmc)){
    # Runs Iterator
    out <- iboom_iterator(Say   = Say,
                          Sgy   = Sgy,
                          Syy   = Syy,
                          A     = A,
                          G     = G,
                          y     = y,
                          Theta = out$Theta,
                          B     = out$B,
                          t2T   = out$t2T,
                          l2T   = out$l2T,
                          xiT   = out$xiT,
                          vT    = out$vT,
                          t2B   = out$t2B,
                          l2B   = out$l2B,
                          xiB   = out$xiB,
                          vB    = out$vB,
                          s2    = out$s2,
                          g     = out$g)
    # Saves Samples
    sTheta[s,,] <- out$Theta
    sB[s,,]     <- out$B
    sg[s,]      <- out$g
    sgp[s,]     <- out$gp
    ss2[s]      <- out$s2
    sl2T[s,,]   <- out$l2T
    st2T[s]     <- out$t2T
    svT[s,,]    <- out$vT
    sxiT[s]     <- out$xiT
    sl2B[s,,]   <- out$l2B
    st2B[s]     <- out$t2B
    svB[s,,]    <- out$vB
    sxiB[s]     <- out$xiB
    
    # Progress Bar Update
    setTxtProgressBar(pb    = pb,
                      value = s / (burnin + nmcmc))
  }
  
  # Splits Burn-In and After-Burin Samples
  sam <- list(Theta = sTheta[(burnin + 1):(burnin + nmcmc),,],
              B     = sB[(burnin + 1):(burnin + nmcmc),,],
              t2T   = st2T[(burnin + 1):(burnin + nmcmc)],
              l2T   = sl2T[(burnin + 1):(burnin + nmcmc),,],
              xiT   = sxiT[(burnin + 1):(burnin + nmcmc)],
              vT    = svT[(burnin + 1):(burnin + nmcmc),,],
              t2B   = st2B[(burnin + 1):(burnin + nmcmc)],
              l2B   = sl2B[(burnin + 1):(burnin + nmcmc),,],
              xiB   = sxiB[(burnin + 1):(burnin + nmcmc)],
              vB    = svB[(burnin + 1):(burnin + nmcmc),,],
              s2    = ss2[(burnin + 1):(burnin + nmcmc)],
              g     = sg[(burnin + 1):(burnin + nmcmc),],
              gp    = sgp[(burnin + 1):(burnin + nmcmc),])
  if(burnin > 0){
    bin <- list(Theta = sTheta[1:burnin,,],
                B     = sB[1:burnin,,],
                t2T   = st2T[1:burnin],
                l2T   = sl2T[1:burnin,,],
                xiT   = sxiT[1:burnin],
                vT    = svT[1:burnin,,],
                t2B   = st2B[1:burnin],
                l2B   = sl2B[1:burnin,,],
                xiB   = sxiB[1:burnin],
                vB    = svB[1:burnin,,],
                s2    = ss2[1:burnin],
                g     = sg[1:burnin,],
                gp    = sgp[1:burnin,])
  } else {
    bin <- list()
  }
  
  # Returns Values
  return(list(sam = sam,
              bin = bin,
              nmcmc = nmcmc))
}