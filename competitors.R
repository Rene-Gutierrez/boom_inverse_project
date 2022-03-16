source("./idata_generator.R")
source("./iboom_iterator.R")
source("./iboom_sampler.R")
source("./group_iterator.R")
source("./sas_iterator.R")
source("./sas_sampler.R")
source("./mreg.R")
source("./ihorseshoe_iterator.R")
source("./ihorseshoe_sampler.R")

iteNum <- 1000

# TPR and FDR Storing
aTPR <-      array(data = NA,
                   dim  = c(3, iteNum, 10))
aFDR <-      array(data = NA,
                   dim  = c(3, iteNum, 10))
aTNR <-      array(data = NA,
                   dim  = c(3, iteNum, 10))
aT_TP_MSE <- matrix(data = NA,
                    nrow = iteNum,
                    ncol = 10)
aT_FP_MSE <- matrix(data = NA,
                    nrow = iteNum,
                    ncol = 10)
aT_TN_MSE <- matrix(data = NA,
                    nrow = iteNum,
                    ncol = 10)
aT_FN_MSE <- matrix(data = NA,
                    nrow = iteNum,
                    ncol = 10)
aT_MSE    <- matrix(data = NA,
                    nrow = iteNum,
                    ncol = 10)
aB_TP_MSE <- matrix(data = NA,
                    nrow = iteNum,
                    ncol = 10)
aB_FP_MSE <- matrix(data = NA,
                    nrow = iteNum,
                    ncol = 10)
aB_TN_MSE <- matrix(data = NA,
                    nrow = iteNum,
                    ncol = 10)
aB_FN_MSE <- matrix(data = NA,
                    nrow = iteNum,
                    ncol = 10)
aB_MSE    <- matrix(data = NA,
                    nrow = iteNum,
                    ncol = 10)
aMSE      <- array(data = NA,
                   dim  = c(3, iteNum, 10))
aLab <- c("No Correction",
          "Global Bonferroni",
          "Regional Bonferroni 1",
          "Regional Bonferroni 2",
          "Global GLM",
          "Regional GLM 1",
          "Regional GLM 2",
          "Spike and Slab",
          "Horseshoe",
          "iBoom")
# General Simulation Settings
PP <- 20
VV <- 100
cB <- 1
cT <- 1
n  <- 16
pB <- 0.2
pT <- 0.2
cB1 <- 2^(-2)
cT1 <- 2^(-2)
cB2 <- 2^(0)
cT2 <- 2^(0)
s2  <- 1
q0  <- 0.05

# Iterations
for (ite in 1:iteNum){
  print("")
  print(ite)
  ### Data Generation
  dat <- idata_generator(PP     = PP,
                         VV     = VV,
                         pB     = pB,
                         pT     = pT,
                         cB     = c(cB1, cB2),
                         cT     = c(cT1, cT2),
                         s2     = s2,
                         n      = n,
                         method = 'random')
  
  # Computes the Multiple Regression
  mr <- mreg(A     = dat$A,
             G     = dat$G,
             y     = dat$y,
             Theta = dat$Theta,
             B     = dat$B)
  # Obtains the Results
  pval  <- mr$pval
  tpval <- mr$tpval
  bpval <- mr$bpval
  that  <- mr$that
  bhat  <- mr$bhat
  
  #### No Corrections
  ##  Both Objects
  dp  <- colSums(pval < q0, na.rm = TRUE) > 0
  #   Computes the Rates and different
  out <- rates(sel = dp, pos = dat$gT)
  #   Saves Values
  aTNR[1, ite, 1] <- out$tnr
  aTPR[1, ite, 1] <- out$tpr
  aFDR[1, ite, 1] <- out$fdr
  ##  Region Level Network Matrix
  dp  <- colSums(tpval < q0, na.rm = TRUE) > 0
  #   Computes the Rates and different
  out <- rates(sel = dp, pos = dat$gT)
  #   Saves Values
  aTNR[2, ite, 1] <- out$tnr
  aTPR[2, ite, 1] <- out$tpr
  aFDR[2, ite, 1] <- out$fdr
  ##  Region Level Grey Matter
  dp  <- colSums(bpval < q0, na.rm = TRUE) > 0
  #   Computes the Rates and different
  out <- rates(sel = dp, pos = dat$gT)
  #   Saves Values
  aTNR[3, ite, 1] <- out$tnr
  aTPR[3, ite, 1] <- out$tpr
  aFDR[3, ite, 1] <- out$fdr
  ### MSE Analysis
  # MSE Theta
  sel <- dp
  pos <- dat$gT == 1
  t_tp_se <- (that[upper.tri(that)][sel * pos == 1] - 
                dat$Theta[upper.tri(dat$Theta)][sel * pos == 1])^2
  t_fp_se <- (that[upper.tri(that)][sel * neg == 1] - 
                dat$Theta[upper.tri(dat$Theta)][sel * neg == 1])^2
  t_tn_se <- (that[upper.tri(that)][uns * neg == 1] - 
                dat$Theta[upper.tri(dat$Theta)][uns * neg == 1])^2
  t_fn_se <- (that[upper.tri(that)][uns * pos == 1] - 
                dat$Theta[upper.tri(dat$Theta)][uns * pos == 1])^2
  aT_TP_MSE[ite, 1] <- mean(t_tp_se)
  aT_FP_MSE[ite, 1] <- mean(t_fp_se)
  aT_TN_MSE[ite, 1] <- mean(t_tn_se)
  aT_FN_MSE[ite, 1] <- mean(t_fn_se)
  aT_MSE[ite, 1]    <- mean(c(t_tp_se, t_fp_se, t_tn_se, t_fn_se))
  # MSE B
  sel <- bpval <  q0
  uns <- bpval >= q0
  pos <- dat$B != 0
  neg <- dat$B == 0
  b_tp_se <- (bhat[sel * pos == 1] - dat$B[sel * pos == 1])^2
  b_fp_se <- (bhat[sel * neg == 1] - dat$B[sel * neg == 1])^2
  b_tn_se <- (bhat[uns * neg == 1] - dat$B[uns * neg == 1])^2
  b_fn_se <- (bhat[uns * pos == 1] - dat$B[uns * pos == 1])^2
  aB_TP_MSE[ite, 1] <- mean(b_tp_se)
  aB_FP_MSE[ite, 1] <- mean(b_fp_se)
  aB_TN_MSE[ite, 1] <- mean(b_tn_se)
  aB_FN_MSE[ite, 1] <- mean(b_fn_se)
  aB_MSE[ite, 1]    <- mean(c(b_tp_se, b_fp_se, b_tn_se, b_fn_se))
  # Global
  aMSE[1, ite, 1] <- mean(c(t_tp_se, t_fp_se, t_tn_se, t_fn_se,
                            b_tp_se, b_fp_se, b_tn_se, b_fn_se))
  aMSE[2, ite, 1] <- mean(c(t_tp_se, t_fp_se, t_tn_se, t_fn_se))
  aMSE[3, ite, 1] <- mean(c(b_tp_se, b_fp_se, b_tn_se, b_fn_se))
  
  ### Global Bonferroni Both Objects
  dp  <- colSums(rbind(tpval < q0 / (PP * (PP - 1) / 2 + PP * VV),
                       bpval < q0 / (PP * (PP - 1) / 2 + PP * VV)), na.rm = TRUE) > 0
  rp  <- dat$gT
  tp  <- dp * rp
  tpr <- sum(tp) / sum(rp)
  fp  <- dp * (1 - rp)
  if(sum(dp) > 0){
    fdr <- sum(fp) / sum(dp)
  } else {
    fdr <- 0
  }
  dn  <- 1 - dp
  rn  <- 1 - tp
  tn  <- dn * rn
  tnr <- sum(tn) / sum(rn)
  aTNR[1, ite, 2] <- tnr
  aTPR[1, ite, 2] <- tpr
  aFDR[1, ite, 2] <- fdr
  ### Global Bonferroni Network Matrix
  dp  <- colSums(tpval < q0 / (PP * (PP - 1) / 2 + PP * VV), na.rm = TRUE) > 0
  rp  <- dat$gT
  tp  <- dp * rp
  tpr <- sum(tp) / sum(rp)
  fp  <- dp * (1 - rp)
  if(sum(dp) > 0){
    fdr <- sum(fp) / sum(dp)
  } else {
    fdr <- 0
  }
  dn  <- 1 - dp
  rn  <- 1 - tp
  tn  <- dn * rn
  tnr <- sum(tn) / sum(rn)
  aTNR[2, ite, 2] <- tnr
  aTPR[2, ite, 2] <- tpr
  aFDR[2, ite, 2] <- fdr
  ### Global Bonferroni Grey Matter
  dp  <- colSums(bpval < q0 / (PP * (PP - 1) / 2 + PP * VV), na.rm = TRUE) > 0
  rp  <- dat$gT
  tp  <- dp * rp
  tpr <- sum(tp) / sum(rp)
  fp  <- dp * (1 - rp)
  if(sum(dp) > 0){
    fdr <- sum(fp) / sum(dp)
  } else {
    fdr <- 0
  }
  dn  <- 1 - dp
  rn  <- 1 - tp
  tn  <- dn * rn
  tnr <- sum(tn) / sum(rn)
  aTNR[3, ite, 2] <- tnr
  aTPR[3, ite, 2] <- tpr
  aFDR[3, ite, 2] <- fdr
  ### MSE Analysis
  # MSE Theta
  sel <- tpval[upper.tri(tpval)] <  q0 / (PP * (PP - 1) / 2 + PP * VV)
  uns <- tpval[upper.tri(tpval)] >= q0 / (PP * (PP - 1) / 2 + PP * VV)
  pos <- dat$Theta[upper.tri(dat$Theta)] != 0
  neg <- dat$Theta[upper.tri(dat$Theta)] == 0
  t_tp_se <- (that[upper.tri(that)][sel * pos == 1] - 
                dat$Theta[upper.tri(dat$Theta)][sel * pos == 1])^2
  t_fp_se <- (that[upper.tri(that)][sel * neg == 1] - 
                dat$Theta[upper.tri(dat$Theta)][sel * neg == 1])^2
  t_tn_se <- (that[upper.tri(that)][uns * neg == 1] - 
                dat$Theta[upper.tri(dat$Theta)][uns * neg == 1])^2
  t_fn_se <- (that[upper.tri(that)][uns * pos == 1] - 
                dat$Theta[upper.tri(dat$Theta)][uns * pos == 1])^2
  aT_TP_MSE[ite, 2] <- mean(t_tp_se)
  aT_FP_MSE[ite, 2] <- mean(t_fp_se)
  aT_TN_MSE[ite, 2] <- mean(t_tn_se)
  aT_FN_MSE[ite, 2] <- mean(t_fn_se)
  aT_MSE[ite, 2]    <- mean(c(t_tp_se, t_fp_se, t_tn_se, t_fn_se))
  # MSE B
  sel <- bpval <  q0 / (PP * (PP - 1) / 2 + PP * VV)
  uns <- bpval >= q0 / (PP * (PP - 1) / 2 + PP * VV)
  pos <- dat$B != 0
  neg <- dat$B == 0
  b_tp_se <- (bhat[sel * pos == 1] - dat$B[sel * pos == 1])^2
  b_fp_se <- (bhat[sel * neg == 1] - dat$B[sel * neg == 1])^2
  b_tn_se <- (bhat[uns * neg == 1] - dat$B[uns * neg == 1])^2
  b_fn_se <- (bhat[uns * pos == 1] - dat$B[uns * pos == 1])^2
  aB_TP_MSE[ite, 2] <- mean(b_tp_se)
  aB_FP_MSE[ite, 2] <- mean(b_fp_se)
  aB_TN_MSE[ite, 2] <- mean(b_tn_se)
  aB_FN_MSE[ite, 2] <- mean(b_fn_se)
  aB_MSE[ite, 2]    <- mean(c(b_tp_se, b_fp_se, b_tn_se, b_fn_se))
  # Global
  aMSE[1, ite, 2] <- mean(c(t_tp_se, t_fp_se, t_tn_se, t_fn_se,
                            b_tp_se, b_fp_se, b_tn_se, b_fn_se))
  aMSE[2, ite, 2] <- mean(c(t_tp_se, t_fp_se, t_tn_se, t_fn_se))
  aMSE[3, ite, 2] <- mean(c(b_tp_se, b_fp_se, b_tn_se, b_fn_se))
  
  ### Regional Bonferroni Version 1 Both Regions
  dp  <- colSums(rbind(tpval < q0 / (PP + VV),
                       bpval < q0 / (PP + VV)), na.rm = TRUE) > 0
  rp  <- dat$gT
  tp  <- dp * rp
  tpr <- sum(tp) / sum(rp)
  fp  <- dp * (1 - rp) 
  if(sum(dp) > 0){
    fdr <- sum(fp) / sum(dp)
  } else {
    fdr <- 0
  }
  dn  <- 1 - dp
  rn  <- 1 - tp
  tn  <- dn * rn
  tnr <- sum(tn) / sum(rn)
  aTNR[1, ite, 3] <- tnr
  aTPR[1, ite, 3] <- tpr
  aFDR[1, ite, 3] <- fdr
  ### Regional Bonferroni Version 1 Network Matrix
  dp  <- colSums(tpval < q0 / (PP + VV), na.rm = TRUE) > 0
  rp  <- dat$gT
  tp  <- dp * rp
  tpr <- sum(tp) / sum(rp)
  fp  <- dp * (1 - rp) 
  if(sum(dp) > 0){
    fdr <- sum(fp) / sum(dp)
  } else {
    fdr <- 0
  }
  dn  <- 1 - dp
  rn  <- 1 - tp
  tn  <- dn * rn
  tnr <- sum(tn) / sum(rn)
  aTNR[2, ite, 3] <- tnr
  aTPR[2, ite, 3] <- tpr
  aFDR[2, ite, 3] <- fdr
  ### Regional Bonferroni Version 1 Grey Matter
  dp  <- colSums(bpval < q0 / (PP + VV), na.rm = TRUE) > 0
  rp  <- dat$gT
  tp  <- dp * rp
  tpr <- sum(tp) / sum(rp)
  fp  <- dp * (1 - rp) 
  if(sum(dp) > 0){
    fdr <- sum(fp) / sum(dp)
  } else {
    fdr <- 0
  }
  dn  <- 1 - dp
  rn  <- 1 - tp
  tn  <- dn * rn
  tnr <- sum(tn) / sum(rn)
  aTNR[3, ite, 3] <- tnr
  aTPR[3, ite, 3] <- tpr
  aFDR[3, ite, 3] <- fdr
  ### MSE Analysis
  # MSE Theta
  sel <- tpval[upper.tri(tpval)] <  q0 / (PP + VV)
  uns <- tpval[upper.tri(tpval)] >= q0 / (PP + VV)
  pos <- dat$Theta[upper.tri(dat$Theta)] != 0
  neg <- dat$Theta[upper.tri(dat$Theta)] == 0
  t_tp_se <- (that[upper.tri(that)][sel * pos == 1] - 
                dat$Theta[upper.tri(dat$Theta)][sel * pos == 1])^2
  t_fp_se <- (that[upper.tri(that)][sel * neg == 1] - 
                dat$Theta[upper.tri(dat$Theta)][sel * neg == 1])^2
  t_tn_se <- (that[upper.tri(that)][uns * neg == 1] - 
                dat$Theta[upper.tri(dat$Theta)][uns * neg == 1])^2
  t_fn_se <- (that[upper.tri(that)][uns * pos == 1] - 
                dat$Theta[upper.tri(dat$Theta)][uns * pos == 1])^2
  aT_TP_MSE[ite, 3] <- mean(t_tp_se)
  aT_FP_MSE[ite, 3] <- mean(t_fp_se)
  aT_TN_MSE[ite, 3] <- mean(t_tn_se)
  aT_FN_MSE[ite, 3] <- mean(t_fn_se)
  aT_MSE[ite, 3]    <- mean(c(t_tp_se, t_fp_se, t_tn_se, t_fn_se))
  # MSE B
  sel <- bpval <  q0 / (PP + VV)
  uns <- bpval >= q0 / (PP + VV)
  pos <- dat$B != 0
  neg <- dat$B == 0
  b_tp_se <- (bhat[sel * pos == 1] - dat$B[sel * pos == 1])^2
  b_fp_se <- (bhat[sel * neg == 1] - dat$B[sel * neg == 1])^2
  b_tn_se <- (bhat[uns * neg == 1] - dat$B[uns * neg == 1])^2
  b_fn_se <- (bhat[uns * pos == 1] - dat$B[uns * pos == 1])^2
  aB_TP_MSE[ite, 3] <- mean(b_tp_se)
  aB_FP_MSE[ite, 3] <- mean(b_fp_se)
  aB_TN_MSE[ite, 3] <- mean(b_tn_se)
  aB_FN_MSE[ite, 3] <- mean(b_fn_se)
  aB_MSE[ite, 3]    <- mean(c(b_tp_se, b_fp_se, b_tn_se, b_fn_se))
  # Global
  aMSE[1, ite, 3] <- mean(c(t_tp_se, t_fp_se, t_tn_se, t_fn_se,
                            b_tp_se, b_fp_se, b_tn_se, b_fn_se))
  aMSE[2, ite, 3] <- mean(c(t_tp_se, t_fp_se, t_tn_se, t_fn_se))
  aMSE[3, ite, 3] <- mean(c(b_tp_se, b_fp_se, b_tn_se, b_fn_se))
  
  ### Regional Bonferroni Version 2 Both Regions
  dp  <- colSums(rbind(tpval < q0 / (PP / 2 + VV),
                       bpval < q0 / (PP / 2 + VV)), na.rm = TRUE) > 0
  rp  <- dat$gT
  tp  <- dp * rp
  tpr <- sum(tp) / sum(rp)
  fp  <- dp * (1 - rp) 
  if(sum(dp) > 0){
    fdr <- sum(fp) / sum(dp)
  } else {
    fdr <- 0
  }
  dn  <- 1 - dp
  rn  <- 1 - tp
  tn  <- dn * rn
  tnr <- sum(tn) / sum(rn)
  aTNR[1, ite, 4] <- tnr
  aTPR[1, ite, 4] <- tpr
  aFDR[1, ite, 4] <- fdr
  ### Regional Bonferroni Version 2 Network Matrix
  dp  <- colSums(tpval < q0 / (PP / 2 + VV), na.rm = TRUE) > 0
  rp  <- dat$gT
  tp  <- dp * rp
  tpr <- sum(tp) / sum(rp)
  fp  <- dp * (1 - rp) 
  if(sum(dp) > 0){
    fdr <- sum(fp) / sum(dp)
  } else {
    fdr <- 0
  }
  dn  <- 1 - dp
  rn  <- 1 - tp
  tn  <- dn * rn
  tnr <- sum(tn) / sum(rn)
  aTNR[2, ite, 4] <- tnr
  aTPR[2, ite, 4] <- tpr
  aFDR[2, ite, 4] <- fdr
  ### Regional Bonferroni Version 2 Grey Matter
  dp  <- colSums(bpval < q0 / (PP / 2 + VV), na.rm = TRUE) > 0
  rp  <- dat$gT
  tp  <- dp * rp
  tpr <- sum(tp) / sum(rp)
  fp  <- dp * (1 - rp) 
  if(sum(dp) > 0){
    fdr <- sum(fp) / sum(dp)
  } else {
    fdr <- 0
  }
  dn  <- 1 - dp
  rn  <- 1 - tp
  tn  <- dn * rn
  tnr <- sum(tn) / sum(rn)
  aTNR[3, ite, 4] <- tnr
  aTPR[3, ite, 4] <- tpr
  aFDR[3, ite, 4] <- fdr
  ### MSE Analysis
  # MSE Theta
  sel <- tpval[upper.tri(tpval)] <  q0 / (PP / 2 + VV)
  uns <- tpval[upper.tri(tpval)] >= q0 / (PP / 2 + VV)
  pos <- dat$Theta[upper.tri(dat$Theta)] != 0
  neg <- dat$Theta[upper.tri(dat$Theta)] == 0
  t_tp_se <- (that[upper.tri(that)][sel * pos == 1] - 
                dat$Theta[upper.tri(dat$Theta)][sel * pos == 1])^2
  t_fp_se <- (that[upper.tri(that)][sel * neg == 1] - 
                dat$Theta[upper.tri(dat$Theta)][sel * neg == 1])^2
  t_tn_se <- (that[upper.tri(that)][uns * neg == 1] - 
                dat$Theta[upper.tri(dat$Theta)][uns * neg == 1])^2
  t_fn_se <- (that[upper.tri(that)][uns * pos == 1] - 
                dat$Theta[upper.tri(dat$Theta)][uns * pos == 1])^2
  aT_TP_MSE[ite, 4] <- mean(t_tp_se)
  aT_FP_MSE[ite, 4] <- mean(t_fp_se)
  aT_TN_MSE[ite, 4] <- mean(t_tn_se)
  aT_FN_MSE[ite, 4] <- mean(t_fn_se)
  aT_MSE[ite, 4]    <- mean(c(t_tp_se, t_fp_se, t_tn_se, t_fn_se))
  # MSE B
  sel <- bpval <  q0 / (PP / 2 + VV)
  uns <- bpval >= q0 / (PP / 2 + VV)
  pos <- dat$B != 0
  neg <- dat$B == 0
  b_tp_se <- (bhat[sel * pos == 1] - dat$B[sel * pos == 1])^2
  b_fp_se <- (bhat[sel * neg == 1] - dat$B[sel * neg == 1])^2
  b_tn_se <- (bhat[uns * neg == 1] - dat$B[uns * neg == 1])^2
  b_fn_se <- (bhat[uns * pos == 1] - dat$B[uns * pos == 1])^2
  aB_TP_MSE[ite, 4] <- mean(b_tp_se)
  aB_FP_MSE[ite, 4] <- mean(b_fp_se)
  aB_TN_MSE[ite, 4] <- mean(b_tn_se)
  aB_FN_MSE[ite, 4] <- mean(b_fn_se)
  aB_MSE[ite, 4]    <- mean(c(b_tp_se, b_fp_se, b_tn_se, b_fn_se))
  # Global
  aMSE[1, ite, 4] <- mean(c(t_tp_se, t_fp_se, t_tn_se, t_fn_se,
                            b_tp_se, b_fp_se, b_tn_se, b_fn_se))
  aMSE[2, ite, 4] <- mean(c(t_tp_se, t_fp_se, t_tn_se, t_fn_se))
  aMSE[3, ite, 4] <- mean(c(b_tp_se, b_fp_se, b_tn_se, b_fn_se))
  
  ### Global GLM Both Regions
  Opval <- c(tpval[up], bpval)[order(c(tpval[up], bpval))]
  npval <- length(Opval)
  compv <- (1:npval) * q0 / npval
  cutof <- Opval < compv
  cutof <- (1:npval)[cutof]
  cutof <- max(Opval[cutof], warning = FALSE)
  dp    <- colSums(rbind(tpval <= cutof,
                         bpval <= cutof), na.rm = TRUE) > 0
  rp    <- dat$gT
  tp    <- dp * rp
  tpr   <- sum(tp) / sum(rp)
  fp    <- dp * (1 - rp) 
  if(sum(dp) > 0){
    fdr <- sum(fp) / sum(dp)
  } else {
    fdr <- 0
  }
  dn  <- 1 - dp
  rn  <- 1 - tp
  tn  <- dn * rn
  tnr <- sum(tn) / sum(rn)
  aTNR[1, ite, 5] <- tnr
  aTPR[1, ite, 5] <- tpr
  aFDR[1, ite, 5] <- fdr
  ### Global GLM Network Matrix
  Opval <- tpval[up][order(tpval[up])]
  npval <- length(Opval)
  compv <- (1:npval) * q0 / npval
  cutof <- Opval < compv
  cutof <- (1:npval)[cutof]
  cutof <- max(Opval[cutof], warning = FALSE)
  dp    <- colSums(tpval <= cutof, na.rm = TRUE) > 0
  rp    <- dat$gT
  tp    <- dp * rp
  tpr   <- sum(tp) / sum(rp)
  fp    <- dp * (1 - rp) 
  if(sum(dp) > 0){
    fdr <- sum(fp) / sum(dp)
  } else {
    fdr <- 0
  }
  dn  <- 1 - dp
  rn  <- 1 - tp
  tn  <- dn * rn
  tnr <- sum(tn) / sum(rn)
  aTNR[2, ite, 5] <- tnr
  aTPR[2, ite, 5] <- tpr
  aFDR[2, ite, 5] <- fdr
  ### Global GLM Grey Matter
  Opval <- bpval[order(bpval)]
  npval <- length(Opval)
  compv <- (1:npval) * q0 / npval
  cutof <- Opval < compv
  cutof <- (1:npval)[cutof]
  cutof <- max(Opval[cutof], warning = FALSE)
  dp    <- colSums(bpval <= cutof, na.rm = TRUE) > 0
  rp    <- dat$gT
  tp    <- dp * rp
  tpr   <- sum(tp) / sum(rp)
  fp    <- dp * (1 - rp) 
  if(sum(dp) > 0){
    fdr <- sum(fp) / sum(dp)
  } else {
    fdr <- 0
  }
  dn  <- 1 - dp
  rn  <- 1 - tp
  tn  <- dn * rn
  tnr <- sum(tn) / sum(rn)
  aTNR[3, ite, 5] <- tnr
  aTPR[3, ite, 5] <- tpr
  aFDR[3, ite, 5] <- fdr
  ### MSE Analysis
  # MSE Theta
  sel <- tpval[upper.tri(tpval)] <  cutof
  uns <- tpval[upper.tri(tpval)] >= cutof
  pos <- dat$Theta[upper.tri(dat$Theta)] != 0
  neg <- dat$Theta[upper.tri(dat$Theta)] == 0
  t_tp_se <- (that[upper.tri(that)][sel * pos == 1] - 
                dat$Theta[upper.tri(dat$Theta)][sel * pos == 1])^2
  t_fp_se <- (that[upper.tri(that)][sel * neg == 1] - 
                dat$Theta[upper.tri(dat$Theta)][sel * neg == 1])^2
  t_tn_se <- (that[upper.tri(that)][uns * neg == 1] - 
                dat$Theta[upper.tri(dat$Theta)][uns * neg == 1])^2
  t_fn_se <- (that[upper.tri(that)][uns * pos == 1] - 
                dat$Theta[upper.tri(dat$Theta)][uns * pos == 1])^2
  aT_TP_MSE[ite, 5] <- mean(t_tp_se)
  aT_FP_MSE[ite, 5] <- mean(t_fp_se)
  aT_TN_MSE[ite, 5] <- mean(t_tn_se)
  aT_FN_MSE[ite, 5] <- mean(t_fn_se)
  aT_MSE[ite, 5]    <- mean(c(t_tp_se, t_fp_se, t_tn_se, t_fn_se))
  # MSE B
  sel <- bpval <  cutof
  uns <- bpval >= cutof
  pos <- dat$B != 0
  neg <- dat$B == 0
  b_tp_se <- (bhat[sel * pos == 1] - dat$B[sel * pos == 1])^2
  b_fp_se <- (bhat[sel * neg == 1] - dat$B[sel * neg == 1])^2
  b_tn_se <- (bhat[uns * neg == 1] - dat$B[uns * neg == 1])^2
  b_fn_se <- (bhat[uns * pos == 1] - dat$B[uns * pos == 1])^2
  aB_TP_MSE[ite, 5] <- mean(b_tp_se)
  aB_FP_MSE[ite, 5] <- mean(b_fp_se)
  aB_TN_MSE[ite, 5] <- mean(b_tn_se)
  aB_FN_MSE[ite, 5] <- mean(b_fn_se)
  aB_MSE[ite, 5]    <- mean(c(b_tp_se, b_fp_se, b_tn_se, b_fn_se))
  # Global
  aMSE[1, ite, 5] <- mean(c(t_tp_se, t_fp_se, t_tn_se, t_fn_se,
                            b_tp_se, b_fp_se, b_tn_se, b_fn_se))
  aMSE[2, ite, 5] <- mean(c(t_tp_se, t_fp_se, t_tn_se, t_fn_se))
  aMSE[3, ite, 5] <- mean(c(b_tp_se, b_fp_se, b_tn_se, b_fn_se))
  
  ### Regional GLM 1 Both Regions
  pval  <- rbind(tpval, bpval)
  Opval <- apply(X = pval, MARGIN = 2, FUN = order)
  for(i in 1:PP){
    Opval[, i] <- pval[, i][Opval[, i]]
  }
  npval <- dim(pval)[1]
  compv <- (1:npval) * q0 / npval
  cutof <- Opval < compv
  pvcut <- numeric(length = PP)
  for(i in 1:PP){
    pvcut[i] <- max(Opval[, i][cutof[, i]], warning = FALSE, na.rm = TRUE)
  }
  cSel  <- t(t(pval) <= pvcut)
  tSel  <- cSel[1:PP, 1:PP]
  tSel  <- tSel & t(tSel)
  cSel[1:PP, 1:PP] <- tSel
  dp    <- colSums(cSel, na.rm = TRUE) > 0
  rp    <- dat$gT
  tp    <- dp * rp
  tpr   <- sum(tp) / sum(rp)
  fp    <- dp * (1 - rp) 
  if(sum(dp) > 0){
    fdr <- sum(fp) / sum(dp)
  } else {
    fdr <- 0
  }
  dn  <- 1 - dp
  rn  <- 1 - tp
  tn  <- dn * rn
  tnr <- sum(tn) / sum(rn)
  aTNR[1, ite, 6] <- tnr
  aTPR[1, ite, 6] <- tpr
  aFDR[1, ite, 6] <- fdr
  ### Regional GLM 1 Network Matrix
  pval  <- tpval
  Opval <- apply(X = pval, MARGIN = 2, FUN = order)
  for(i in 1:PP){
    Opval[, i] <- pval[, i][Opval[, i]]
  }
  npval <- dim(pval)[1]
  compv <- (1:npval) * q0 / npval
  cutof <- Opval < compv
  pvcut <- numeric(length = PP)
  for(i in 1:PP){
    pvcut[i] <- max(Opval[, i][cutof[, i]], warning = FALSE, na.rm = TRUE)
  }
  cSel  <- t(t(pval) <= pvcut)
  dp    <- colSums(cSel, na.rm = TRUE) > 0
  rp    <- dat$gT
  tp    <- dp * rp
  tpr   <- sum(tp) / sum(rp)
  fp    <- dp * (1 - rp) 
  if(sum(dp) > 0){
    fdr <- sum(fp) / sum(dp)
  } else {
    fdr <- 0
  }
  dn  <- 1 - dp
  rn  <- 1 - tp
  tn  <- dn * rn
  tnr <- sum(tn) / sum(rn)
  aTNR[2, ite, 6] <- tnr
  aTPR[2, ite, 6] <- tpr
  aFDR[2, ite, 6] <- fdr
  ### Regional GLM 1 Grey Matter
  pval  <- bpval
  Opval <- apply(X = pval, MARGIN = 2, FUN = order)
  for(i in 1:PP){
    Opval[, i] <- pval[, i][Opval[, i]]
  }
  npval <- dim(pval)[1]
  compv <- (1:npval) * q0 / npval
  cutof <- Opval < compv
  pvcut <- numeric(length = PP)
  for(i in 1:PP){
    pvcut[i] <- max(Opval[, i][cutof[, i]], warning = FALSE, na.rm = TRUE)
  }
  cSel  <- t(t(pval) <= pvcut)
  dp    <- colSums(cSel, na.rm = TRUE) > 0
  rp    <- dat$gT
  tp    <- dp * rp
  tpr   <- sum(tp) / sum(rp)
  fp    <- dp * (1 - rp) 
  if(sum(dp) > 0){
    fdr <- sum(fp) / sum(dp)
  } else {
    fdr <- 0
  }
  dn  <- 1 - dp
  rn  <- 1 - tp
  tn  <- dn * rn
  tnr <- sum(tn) / sum(rn)
  aTNR[3, ite, 6] <- tnr
  aTPR[3, ite, 6] <- tpr
  aFDR[3, ite, 6] <- fdr
  ### MSE Analysis
  # MSE Theta
  sel <- tSel[upper.tri(tSel)]
  uns <- (!tSel)[upper.tri(tpval)]
  pos <- dat$Theta[upper.tri(dat$Theta)] != 0
  neg <- dat$Theta[upper.tri(dat$Theta)] == 0
  t_tp_se <- (that[upper.tri(that)][sel * pos == 1] - 
                dat$Theta[upper.tri(dat$Theta)][sel * pos == 1])^2
  t_fp_se <- (that[upper.tri(that)][sel * neg == 1] - 
                dat$Theta[upper.tri(dat$Theta)][sel * neg == 1])^2
  t_tn_se <- (that[upper.tri(that)][uns * neg == 1] - 
                dat$Theta[upper.tri(dat$Theta)][uns * neg == 1])^2
  t_fn_se <- (that[upper.tri(that)][uns * pos == 1] - 
                dat$Theta[upper.tri(dat$Theta)][uns * pos == 1])^2
  aT_TP_MSE[ite, 6] <- mean(t_tp_se)
  aT_FP_MSE[ite, 6] <- mean(t_fp_se)
  aT_TN_MSE[ite, 6] <- mean(t_tn_se)
  aT_FN_MSE[ite, 6] <- mean(t_fn_se)
  aT_MSE[ite, 6]    <- mean(c(t_tp_se, t_fp_se, t_tn_se, t_fn_se))
  # MSE B
  sel <-  cSel
  uns <- !cSel
  pos <- dat$B != 0
  neg <- dat$B == 0
  b_tp_se <- (bhat[sel * pos == 1] - dat$B[sel * pos == 1])^2
  b_fp_se <- (bhat[sel * neg == 1] - dat$B[sel * neg == 1])^2
  b_tn_se <- (bhat[uns * neg == 1] - dat$B[uns * neg == 1])^2
  b_fn_se <- (bhat[uns * pos == 1] - dat$B[uns * pos == 1])^2
  aB_TP_MSE[ite, 6] <- mean(b_tp_se)
  aB_FP_MSE[ite, 6] <- mean(b_fp_se)
  aB_TN_MSE[ite, 6] <- mean(b_tn_se)
  aB_FN_MSE[ite, 6] <- mean(b_fn_se)
  aB_MSE[ite, 6]    <- mean(c(b_tp_se, b_fp_se, b_tn_se, b_fn_se))
  # Global
  aMSE[1, ite, 6] <- mean(c(t_tp_se, t_fp_se, t_tn_se, t_fn_se,
                            b_tp_se, b_fp_se, b_tn_se, b_fn_se))
  aMSE[2, ite, 6] <- mean(c(t_tp_se, t_fp_se, t_tn_se, t_fn_se))
  aMSE[3, ite, 6] <- mean(c(b_tp_se, b_fp_se, b_tn_se, b_fn_se))
  
  ### Regional GLM 2 Both Regions
  pval  <- rbind(tpval, bpval)
  Opval <- apply(X = pval, MARGIN = 2, FUN = order)
  for(i in 1:PP){
    Opval[, i] <- pval[, i][Opval[, i]]
  }
  npval <- dim(pval)[1]
  compv <- (1:npval) * q0 / npval
  cutof <- Opval < compv
  pvcut <- numeric(length = PP)
  for(i in 1:PP){
    pvcut[i] <- max(Opval[, i][cutof[, i]], warning = FALSE, na.rm = TRUE)
  }
  cSel  <- t(t(pval) <= pvcut)
  tSel  <- cSel[1:PP, 1:PP]
  tSel  <- tSel | t(tSel)
  cSel[1:PP, 1:PP] <- tSel
  dp    <- colSums(cSel, na.rm = TRUE) > 0
  rp    <- dat$gT
  tp    <- dp * rp
  tpr   <- sum(tp) / sum(rp)
  fp    <- dp * (1 - rp) 
  if(sum(dp) > 0){
    fdr <- sum(fp) / sum(dp)
  } else {
    fdr <- 0
  }
  dn  <- 1 - dp
  rn  <- 1 - tp
  tn  <- dn * rn
  tnr <- sum(tn) / sum(rn)
  aTNR[1, ite, 7] <- tnr
  aTPR[1, ite, 7] <- tpr
  aFDR[1, ite, 7] <- fdr
  ### Regional GLM 2 Network Matrix
  pval  <- tpval
  Opval <- apply(X = pval, MARGIN = 2, FUN = order)
  for(i in 1:PP){
    Opval[, i] <- pval[, i][Opval[, i]]
  }
  npval <- dim(pval)[1]
  compv <- (1:npval) * q0 / npval
  cutof <- Opval < compv
  pvcut <- numeric(length = PP)
  for(i in 1:PP){
    pvcut[i] <- max(Opval[, i][cutof[, i]], warning = FALSE, na.rm = TRUE)
  }
  cSel  <- t(t(pval) <= pvcut)
  dp    <- colSums(cSel, na.rm = TRUE) > 0
  rp    <- dat$gT
  tp    <- dp * rp
  tpr   <- sum(tp) / sum(rp)
  fp    <- dp * (1 - rp) 
  if(sum(dp) > 0){
    fdr <- sum(fp) / sum(dp)
  } else {
    fdr <- 0
  }
  dn  <- 1 - dp
  rn  <- 1 - tp
  tn  <- dn * rn
  tnr <- sum(tn) / sum(rn)
  aTNR[2, ite, 7] <- tnr
  aTPR[2, ite, 7] <- tpr
  aFDR[2, ite, 7] <- fdr
  ### Regional GLM 2 Grey Matter
  pval  <- bpval
  Opval <- apply(X = pval, MARGIN = 2, FUN = order)
  for(i in 1:PP){
    Opval[, i] <- pval[, i][Opval[, i]]
  }
  npval <- dim(pval)[1]
  compv <- (1:npval) * q0 / npval
  cutof <- Opval < compv
  pvcut <- numeric(length = PP)
  for(i in 1:PP){
    pvcut[i] <- max(Opval[, i][cutof[, i]], warning = FALSE, na.rm = TRUE)
  }
  cSel  <- t(t(pval) <= pvcut)
  dp    <- colSums(cSel, na.rm = TRUE) > 0
  rp    <- dat$gT
  tp    <- dp * rp
  tpr   <- sum(tp) / sum(rp)
  fp    <- dp * (1 - rp) 
  if(sum(dp) > 0){
    fdr <- sum(fp) / sum(dp)
  } else {
    fdr <- 0
  }
  dn  <- 1 - dp
  rn  <- 1 - tp
  tn  <- dn * rn
  tnr <- sum(tn) / sum(rn)
  aTNR[3, ite, 7] <- tnr
  aTPR[3, ite, 7] <- tpr
  aFDR[3, ite, 7] <- fdr
  ### MSE Analysis
  # MSE Theta
  sel <- tSel[upper.tri(tSel)]
  uns <- (!tSel)[upper.tri(tpval)]
  pos <- dat$Theta[upper.tri(dat$Theta)] != 0
  neg <- dat$Theta[upper.tri(dat$Theta)] == 0
  t_tp_se <- (that[upper.tri(that)][sel * pos == 1] - 
                dat$Theta[upper.tri(dat$Theta)][sel * pos == 1])^2
  t_fp_se <- (that[upper.tri(that)][sel * neg == 1] - 
                dat$Theta[upper.tri(dat$Theta)][sel * neg == 1])^2
  t_tn_se <- (that[upper.tri(that)][uns * neg == 1] - 
                dat$Theta[upper.tri(dat$Theta)][uns * neg == 1])^2
  t_fn_se <- (that[upper.tri(that)][uns * pos == 1] - 
                dat$Theta[upper.tri(dat$Theta)][uns * pos == 1])^2
  aT_TP_MSE[ite, 7] <- mean(t_tp_se)
  aT_FP_MSE[ite, 7] <- mean(t_fp_se)
  aT_TN_MSE[ite, 7] <- mean(t_tn_se)
  aT_FN_MSE[ite, 7] <- mean(t_fn_se)
  aT_MSE[ite, 7]    <- mean(c(t_tp_se, t_fp_se, t_tn_se, t_fn_se))
  # MSE B
  sel <-  cSel
  uns <- !cSel
  pos <- dat$B != 0
  neg <- dat$B == 0
  b_tp_se <- (bhat[sel * pos == 1] - dat$B[sel * pos == 1])^2
  b_fp_se <- (bhat[sel * neg == 1] - dat$B[sel * neg == 1])^2
  b_tn_se <- (bhat[uns * neg == 1] - dat$B[uns * neg == 1])^2
  b_fn_se <- (bhat[uns * pos == 1] - dat$B[uns * pos == 1])^2
  aB_TP_MSE[ite, 7] <- mean(b_tp_se)
  aB_FP_MSE[ite, 7] <- mean(b_fp_se)
  aB_TN_MSE[ite, 7] <- mean(b_tn_se)
  aB_FN_MSE[ite, 7] <- mean(b_fn_se)
  aB_MSE[ite, 7]    <- mean(c(b_tp_se, b_fp_se, b_tn_se, b_fn_se))
  # Global
  aMSE[1, ite, 7] <- mean(c(t_tp_se, t_fp_se, t_tn_se, t_fn_se,
                            b_tp_se, b_fp_se, b_tn_se, b_fn_se))
  aMSE[2, ite, 7] <- mean(c(t_tp_se, t_fp_se, t_tn_se, t_fn_se))
  aMSE[3, ite, 7] <- mean(c(b_tp_se, b_fp_se, b_tn_se, b_fn_se))
  
  # Spike and Slab
  # Data Manipulation
  Y <- rbind(matrix(data  = c(dat$A),
                    nrow  = PP * PP,
                    ncol  = n,
                    byrow = TRUE)[c(upper.tri(dat$A[1,,])), ],
             matrix(data  = c(dat$G),
                    nrow  = VV * PP,
                    ncol  = n,
                    byrow = TRUE))
  X <- matrix(data  = dat$y,
              nrow  = (PP * (PP - 1) / 2 + VV * PP),
              ncol  = n,
              byrow = TRUE)
  
  # Sas Sampler
  out <- sas_sampler(Y      = Y,
                     X      = X,
                     burnin = 100,
                     nmcmc  = 1000)
  
  # Recover Estimates
  eb  <- apply(X = out$sam$b, MARGIN = 2, FUN = median)
  exi <- colMeans(x = out$sam$xi)
  # Theta Estimate
  eTheta <- matrix(data = 0, nrow = PP, ncol = PP)
  eTheta[upper.tri(eTheta)] <- eb[1:(PP * (PP - 1) / 2)]
  eTheta                    <- eTheta + t(eTheta) 
  # Indicators Theta
  egT <- matrix(data = 0, nrow = PP, ncol = PP)
  egT[upper.tri(egT)] <- exi[1:(PP * (PP - 1) / 2)]
  egT <- egT + t(egT)
  # B Estimate
  eB <- matrix(data = eb[(PP * (PP - 1) / 2 + 1):(PP * (PP - 1) / 2 + PP * VV)],
               nrow = VV,
               ncol = PP)
  # Indicators B
  egB <- matrix(data = exi[(PP * (PP - 1) / 2 + 1):(PP * (PP - 1) / 2 + PP * VV)],
                nrow = VV,
                ncol = PP)
  
  ### Region Selection Both Regions
  dp  <- colSums(x = rbind(egT, egB) > 0.5) > 0
  rp  <- dat$gT
  tp  <- dp * rp
  tpr <- sum(tp) / sum(rp)
  fp  <- dp * (1 - rp)
  if(sum(dp) > 0){
    fdr <- sum(fp) / sum(dp)
  } else {
    fdr <- 0
  }
  dn  <- 1 - dp
  rn  <- 1 - tp
  tn  <- dn * rn
  tnr <- sum(tn) / sum(rn)
  aTNR[1, ite, 8] <- tnr
  aTPR[1, ite, 8] <- tpr
  aFDR[1, ite, 8] <- fdr
  ### Region Selection Network Matrix
  dp  <- colSums(x = egT > 0.5) > 0
  rp  <- dat$gT
  tp  <- dp * rp
  tpr <- sum(tp) / sum(rp)
  fp  <- dp * (1 - rp)
  if(sum(dp) > 0){
    fdr <- sum(fp) / sum(dp)
  } else {
    fdr <- 0
  }
  dn  <- 1 - dp
  rn  <- 1 - tp
  tn  <- dn * rn
  tnr <- sum(tn) / sum(rn)
  aTNR[2, ite, 8] <- tnr
  aTPR[2, ite, 8] <- tpr
  aFDR[2, ite, 8] <- fdr
  ### Region Selection Grey Matrix
  dp  <- colSums(x = egB > 0.5) > 0
  rp  <- dat$gT
  tp  <- dp * rp
  tpr <- sum(tp) / sum(rp)
  fp  <- dp * (1 - rp)
  if(sum(dp) > 0){
    fdr <- sum(fp) / sum(dp)
  } else {
    fdr <- 0
  }
  dn  <- 1 - dp
  rn  <- 1 - tp
  tn  <- dn * rn
  tnr <- sum(tn) / sum(rn)
  aTNR[3, ite, 8] <- tnr
  aTPR[3, ite, 8] <- tpr
  aFDR[3, ite, 8] <- fdr
  ### MSE Analysis
  # MSE Theta
  sel <- egT[upper.tri(that)] > 0.5
  uns <- egT[upper.tri(that)] <= 0.5
  pos <- dat$Theta[upper.tri(dat$Theta)] != 0
  neg <- dat$Theta[upper.tri(dat$Theta)] == 0
  t_tp_se <- (that[upper.tri(that)][sel * pos == 1] - 
                dat$Theta[upper.tri(dat$Theta)][sel * pos == 1])^2
  t_fp_se <- (that[upper.tri(that)][sel * neg == 1] - 
                dat$Theta[upper.tri(dat$Theta)][sel * neg == 1])^2
  t_tn_se <- (that[upper.tri(that)][uns * neg == 1] - 
                dat$Theta[upper.tri(dat$Theta)][uns * neg == 1])^2
  t_fn_se <- (that[upper.tri(that)][uns * pos == 1] - 
                dat$Theta[upper.tri(dat$Theta)][uns * pos == 1])^2
  aT_TP_MSE[ite, 8] <- mean(t_tp_se)
  aT_FP_MSE[ite, 8] <- mean(t_fp_se)
  aT_TN_MSE[ite, 8] <- mean(t_tn_se)
  aT_FN_MSE[ite, 8] <- mean(t_fn_se)
  aT_MSE[ite, 8]    <- mean(c(t_tp_se, t_fp_se, t_tn_se, t_fn_se))
  # MSE B
  sel <- egB > 0.5
  uns <- egB <= 0.5
  pos <- dat$B != 0
  neg <- dat$B == 0
  b_tp_se <- (bhat[sel * pos == 1] - dat$B[sel * pos == 1])^2
  b_fp_se <- (bhat[sel * neg == 1] - dat$B[sel * neg == 1])^2
  b_tn_se <- (bhat[uns * neg == 1] - dat$B[uns * neg == 1])^2
  b_fn_se <- (bhat[uns * pos == 1] - dat$B[uns * pos == 1])^2
  aB_TP_MSE[ite, 8] <- mean(b_tp_se)
  aB_FP_MSE[ite, 8] <- mean(b_fp_se)
  aB_TN_MSE[ite, 8] <- mean(b_tn_se)
  aB_FN_MSE[ite, 8] <- mean(b_fn_se)
  aB_MSE[ite, 8]    <- mean(c(b_tp_se, b_fp_se, b_tn_se, b_fn_se))
  # Global
  aMSE[1, ite, 8] <- mean(c(t_tp_se, t_fp_se, t_tn_se, t_fn_se,
                            b_tp_se, b_fp_se, b_tn_se, b_fn_se))
  aMSE[2, ite, 8] <- mean(c(t_tp_se, t_fp_se, t_tn_se, t_fn_se))
  aMSE[3, ite, 8] <- mean(c(b_tp_se, b_fp_se, b_tn_se, b_fn_se))
  
  # Global Horseshoe
  # Horseshoe Sampler
  out <- ihorseshoe_sampler(Y      = Y,
                            X      = X,
                            burnin = 0,
                            nmcmc  = 1000)
  # Generates cxi (Via Clustering)
  clu <- kmeans(x = abs(apply(X = out$sam$b, MARGIN = 2, FUN = median)), centers = 2)
  clu1Mea <- mean(abs(out$sam$b)[,clu$cluster == 1])
  clu2Mea <- mean(abs(out$sam$b)[,clu$cluster == 2])
  # Cluster Classification
  xi <- rep(0, out$p)
  if(clu1Mea > clu2Mea){
    xi[clu$cluster == 1] <- 1
  } else {
    xi[clu$cluster == 2] <- 1
  }
  
  # Recover Estimates
  eb  <- apply(X = out$sam$b, MARGIN = 2, FUN = median)
  exi <- xi
  # Theta Estimate
  eTheta <- matrix(data = 0, nrow = PP, ncol = PP)
  eTheta[upper.tri(eTheta)] <- eb[1:(PP * (PP - 1) / 2)]
  eTheta                    <- eTheta + t(eTheta) 
  # Indicators Theta
  egT <- matrix(data = 0, nrow = PP, ncol = PP)
  egT[upper.tri(egT)] <- exi[1:(PP * (PP - 1) / 2)]
  egT <- egT + t(egT)
  # B Estimate
  eB <- matrix(data = eb[(PP * (PP - 1) / 2 + 1):(PP * (PP - 1) / 2 + PP * VV)],
               nrow = VV,
               ncol = PP)
  # Indicators B
  egB <- matrix(data = exi[(PP * (PP - 1) / 2 + 1):(PP * (PP - 1) / 2 + PP * VV)],
                nrow = VV,
                ncol = PP)
  # Global Horseshoe
  ### Region Selection Both Objects
  dp  <- colSums(x = rbind(egT, egB) > 0.5) > 0
  rp  <- dat$gT
  tp  <- dp * rp
  tpr <- sum(tp) / sum(rp)
  fp  <- dp * (1 - rp)
  if(sum(dp) > 0){
    fdr <- sum(fp) / sum(dp)
  } else {
    fdr <- 0
  }
  dn  <- 1 - dp
  rn  <- 1 - tp
  tn  <- dn * rn
  tnr <- sum(tn) / sum(rn)
  aTNR[1, ite, 9] <- tnr
  aTPR[1, ite, 9] <- tpr
  aFDR[1, ite, 9] <- fdr
  ### Region Selection Network Matrix
  dp  <- colSums(x = egT > 0.5) > 0
  rp  <- dat$gT
  tp  <- dp * rp
  tpr <- sum(tp) / sum(rp)
  fp  <- dp * (1 - rp)
  if(sum(dp) > 0){
    fdr <- sum(fp) / sum(dp)
  } else {
    fdr <- 0
  }
  dn  <- 1 - dp
  rn  <- 1 - tp
  tn  <- dn * rn
  tnr <- sum(tn) / sum(rn)
  aTNR[2, ite, 9] <- tnr
  aTPR[2, ite, 9] <- tpr
  aFDR[2, ite, 9] <- fdr
  ### Region Selection Both Objects
  dp  <- colSums(x = egB > 0.5) > 0
  rp  <- dat$gT
  tp  <- dp * rp
  tpr <- sum(tp) / sum(rp)
  fp  <- dp * (1 - rp)
  if(sum(dp) > 0){
    fdr <- sum(fp) / sum(dp)
  } else {
    fdr <- 0
  }
  dn  <- 1 - dp
  rn  <- 1 - tp
  tn  <- dn * rn
  tnr <- sum(tn) / sum(rn)
  aTNR[3, ite, 9] <- tnr
  aTPR[3, ite, 9] <- tpr
  aFDR[3, ite, 9] <- fdr
  ### MSE Analysis
  # MSE Theta
  sel <- egT[upper.tri(that)] > 0.5
  uns <- egT[upper.tri(that)] <= 0.5
  pos <- dat$Theta[upper.tri(dat$Theta)] != 0
  neg <- dat$Theta[upper.tri(dat$Theta)] == 0
  t_tp_se <- (that[upper.tri(that)][sel * pos == 1] - 
                dat$Theta[upper.tri(dat$Theta)][sel * pos == 1])^2
  t_fp_se <- (that[upper.tri(that)][sel * neg == 1] - 
                dat$Theta[upper.tri(dat$Theta)][sel * neg == 1])^2
  t_tn_se <- (that[upper.tri(that)][uns * neg == 1] - 
                dat$Theta[upper.tri(dat$Theta)][uns * neg == 1])^2
  t_fn_se <- (that[upper.tri(that)][uns * pos == 1] - 
                dat$Theta[upper.tri(dat$Theta)][uns * pos == 1])^2
  aT_TP_MSE[ite, 9] <- mean(t_tp_se)
  aT_FP_MSE[ite, 9] <- mean(t_fp_se)
  aT_TN_MSE[ite, 9] <- mean(t_tn_se)
  aT_FN_MSE[ite, 9] <- mean(t_fn_se)
  aT_MSE[ite, 9]    <- mean(c(t_tp_se, t_fp_se, t_tn_se, t_fn_se))
  # MSE B
  sel <- egB > 0.5
  uns <- egB <= 0.5
  pos <- dat$B != 0
  neg <- dat$B == 0
  b_tp_se <- (bhat[sel * pos == 1] - dat$B[sel * pos == 1])^2
  b_fp_se <- (bhat[sel * neg == 1] - dat$B[sel * neg == 1])^2
  b_tn_se <- (bhat[uns * neg == 1] - dat$B[uns * neg == 1])^2
  b_fn_se <- (bhat[uns * pos == 1] - dat$B[uns * pos == 1])^2
  aB_TP_MSE[ite, 9] <- mean(b_tp_se)
  aB_FP_MSE[ite, 9] <- mean(b_fp_se)
  aB_TN_MSE[ite, 9] <- mean(b_tn_se)
  aB_FN_MSE[ite, 9] <- mean(b_fn_se)
  aB_MSE[ite, 9]    <- mean(c(b_tp_se, b_fp_se, b_tn_se, b_fn_se))
  # Global
  aMSE[1, ite, 9] <- mean(c(t_tp_se, t_fp_se, t_tn_se, t_fn_se,
                            b_tp_se, b_fp_se, b_tn_se, b_fn_se))
  aMSE[2, ite, 9] <- mean(c(t_tp_se, t_fp_se, t_tn_se, t_fn_se))
  aMSE[3, ite, 9] <- mean(c(b_tp_se, b_fp_se, b_tn_se, b_fn_se))
  
  #iBOOM
  out <- iboom_sampler(y      = dat$y,
                       G      = dat$G,
                       A      = dat$A,
                       nmcmc  = 100,
                       burnin = 100)
  dp    <- colMeans(out$sam$g) > 0.5
  rp    <- dat$gT
  tp    <- dp * rp
  tpr   <- sum(tp) / sum(rp)
  fp    <- dp * (1 - rp) 
  if(sum(dp) > 0){
    fdr <- sum(fp) / sum(dp)
  } else {
    fdr <- 0
  }
  dn  <- 1 - dp
  rn  <- 1 - tp
  tn  <- dn * rn
  tnr <- sum(tn) / sum(rn)
  aTNR[1, ite, 10] <- tnr
  aTPR[1, ite, 10] <- tpr
  aFDR[1, ite, 10] <- fdr
  ### MSE Analysis
  # MSE Theta
  reg <- colMeans(out$sam$g)
  sel <- (reg %*% t(reg))[upper.tri(that)]
  uns <- (!(reg %*% t(reg)))[upper.tri(that)]
  pos <- dat$Theta[upper.tri(dat$Theta)] != 0
  neg <- dat$Theta[upper.tri(dat$Theta)] == 0
  t_tp_se <- (that[upper.tri(that)][sel * pos == 1] - 
                dat$Theta[upper.tri(dat$Theta)][sel * pos == 1])^2
  t_fp_se <- (that[upper.tri(that)][sel * neg == 1] - 
                dat$Theta[upper.tri(dat$Theta)][sel * neg == 1])^2
  t_tn_se <- (that[upper.tri(that)][uns * neg == 1] - 
                dat$Theta[upper.tri(dat$Theta)][uns * neg == 1])^2
  t_fn_se <- (that[upper.tri(that)][uns * pos == 1] - 
                dat$Theta[upper.tri(dat$Theta)][uns * pos == 1])^2
  aT_TP_MSE[ite, 10] <- mean(t_tp_se)
  aT_FP_MSE[ite, 10] <- mean(t_fp_se)
  aT_TN_MSE[ite, 10] <- mean(t_tn_se)
  aT_FN_MSE[ite, 10] <- mean(t_fn_se)
  aT_MSE[ite, 10]    <- mean(c(t_tp_se, t_fp_se, t_tn_se, t_fn_se))
  # MSE B
  sel <-  rep(1, 100) %*% t(reg)
  uns <- !(rep(1, 100) %*% t(reg))
  pos <- dat$B != 0
  neg <- dat$B == 0
  b_tp_se <- (bhat[sel * pos == 1] - dat$B[sel * pos == 1])^2
  b_fp_se <- (bhat[sel * neg == 1] - dat$B[sel * neg == 1])^2
  b_tn_se <- (bhat[uns * neg == 1] - dat$B[uns * neg == 1])^2
  b_fn_se <- (bhat[uns * pos == 1] - dat$B[uns * pos == 1])^2
  aB_TP_MSE[ite, 10] <- mean(b_tp_se)
  aB_FP_MSE[ite, 10] <- mean(b_fp_se)
  aB_TN_MSE[ite, 10] <- mean(b_tn_se)
  aB_FN_MSE[ite, 10] <- mean(b_fn_se)
  aB_MSE[ite, 10]    <- mean(c(b_tp_se, b_fp_se, b_tn_se, b_fn_se))
  # Global
  aMSE[1, ite, 10] <- mean(c(t_tp_se, t_fp_se, t_tn_se, t_fn_se,
                         b_tp_se, b_fp_se, b_tn_se, b_fn_se))
}

# Saves Data Sets
#saveRDS(object = aTPR, file = paste0("tpr_comp_", n, "_", PP, "_", VV, "_", pT, "_", pB, "_", ".rds"))
#saveRDS(object = aTNR, file = paste0("tnr_comp_", n, "_", PP, "_", VV, "_", pT, "_", pB, "_", ".rds"))
saveRDS(object = aTPR, file = paste0("tpr_comp.rds"))
saveRDS(object = aTNR, file = paste0("tnr_comp.rds"))
saveRDS(object = aLab, file = paste0("lab_comp.rds"))
saveRDS(object = aT_TP_MSE, file = paste0("t_tp_mse_comp.rds"))
saveRDS(object = aT_FP_MSE, file = paste0("t_fp_mse_comp.rds"))
saveRDS(object = aT_TN_MSE, file = paste0("t_tn_mse_comp.rds"))
saveRDS(object = aT_FN_MSE, file = paste0("t_fn_mse_comp.rds"))
saveRDS(object = aB_TP_MSE, file = paste0("b_tp_mse_comp.rds"))
saveRDS(object = aB_FP_MSE, file = paste0("b_fp_mse_comp.rds"))
saveRDS(object = aB_TN_MSE, file = paste0("b_tn_mse_comp.rds"))
saveRDS(object = aB_FN_MSE, file = paste0("b_fn_mse_comp.rds"))
saveRDS(object = aT_MSE, file = paste0("t_mse_comp.rds"))
saveRDS(object = aB_MSE, file = paste0("b_mse_comp.rds"))
saveRDS(object = aMSE, file = paste0("mse_comp.rds"))
