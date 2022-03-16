mreg <- function(A,
                 G,
                 y,
                 B,
                 Theta){
  ### Sufficient Statistics
  Say <- apply(X = (A * y), MARGIN = c(2, 3), sum)
  Sgy <- apply(X = (G * y), MARGIN = c(2, 3), sum)
  Syy <- sum(y^2)
  P   <- dim(G)[3]
  V   <- dim(G)[2]
  n   <- dim(G)[1]
  
  ### Theta Estimates
  # Estimates
  that  <- Say / Syy
  # Auxiliary Variable
  mthat <- aperm(a    = array(data = that,
                              dim  = c(P, P, n)),
                 perm = c(3, 1, 2))
  # Response Estimate
  ahat  <- mthat * y
  # Errors
  terr  <- A - ahat
  # Response Variance Estimate
  tshat <- apply(X      = (terr)^2,
                 MARGIN = c(2, 3),
                 FUN    = sum) / (n - 1) 
  # Estimate SE
  tsehat <- sqrt( tshat / var(y) / (n - 1))
  # t-Statistic
  tstat <- that / tsehat
  # p-values
  tpval <- (1 - pt(q = tstat, df = n - 1))
  
  ### B Estimates
  # Estimates
  bhat  <- Sgy / Syy
  # Auxiliary Variable
  mbhat <- aperm(a    = array(data = bhat,
                              dim  = c(V, P, n)),
                 perm = c(3, 1, 2))
  # Response Estimate
  ghat  <- mbhat * y
  # Errors
  berr  <- G - ghat
  # Response Variance Estimate
  bshat <- apply(X      = (berr)^2,
                 MARGIN = c(2, 3),
                 FUN    = sum) / (n - 1)
  # Estimate SE
  bsehat <- sqrt( bshat / var(y) / (n - 1))
  # t-Statistic
  tstat <- bhat / bsehat
  # p-values
  bpval <- (1 - pt(q = tstat, df = n - 1))
  
  # Joint P-Values
  tpval[lower.tri(x = tpval, diag = TRUE)] <- NA
  pval <- rbind(tpval, bpval)
  
  # Joint Estimates
  that[lower.tri(x = that, diag = TRUE)] <- NA
  hat  <- rbind(that, bhat)
  
  # Errors
  tmse <- (Theta - that)^2
  bmse <- (B     - bhat)^2
  tmse[lower.tri(x = tmse, diag = TRUE)] <- NA
  mse  <- rbind(tmse, bmse)
  
  # CI
  tlci <- that - qnorm(p = 0.975) * tsehat
  tlci[lower.tri(x = tlci, diag = TRUE)] <- NA
  tuci <- that + qnorm(p = 0.975) * tsehat
  tuci[lower.tri(x = tuci, diag = TRUE)] <- NA
  blci <- bhat - qnorm(p = 0.975) * bsehat
  buci <- bhat + qnorm(p = 0.975) * bsehat
  lci  <- rbind(tlci, blci)
  uci  <- rbind(tuci, buci)
  
  # Coverage
  tcov <- (tlci < Theta) & (Theta < tuci)
  bcov <- (blci < B    ) & (B     < buci)
  cov  <- rbind(tcov, bcov)
  
  ### Returns Values
  return(list(bpval = bpval,
              tpval = tpval,
              pval  = pval,
              that  = that,
              bhat  = bhat,
              hat   = hat,
              tmse  = tmse,
              bmse  = bmse,
              mse   = mse,
              tlci  = tlci,
              tuci  = tuci,
              blci  = blci,
              buci  = buci,
              lci   = lci,
              uci   = uci,
              tcov  = tcov,
              bcov  = bcov,
              cov   = cov))
}