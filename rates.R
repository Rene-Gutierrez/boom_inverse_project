rates <- function(sel, pos){
  # Segments
  tp <- (sel   * pos)     == 1
  fp <- (sel   * (!pos))  == 1
  tn <- ((!sel) * (!pos)) == 1
  fn <- ((!sel) * pos)    == 1
  
  # Rates
  # TPR
  if(sum(pos) == 0){
    tpr <- 0
  } else {
    tpr <- sum(tp) / sum(pos)
  }
  # TNR
  if(sum(!pos) == 0){
    tnr <- 0
  } else {
    tnr <- sum(tn) / sum(!pos)
  }
  # FDR
  if(sum(sel) == 0){
    fdr <- 0
  } else {
    fdr <- sum(fp) / sum(sel)
  }
  # Precision
  if(sum(sel) == 0){
    pre <- 0
  } else {
    pre <- sum(tp) / sum(sel)
  }
  # Recall
  if(sum(pos) == 0){
    rec <- 0
  } else {
    rec <- sum(tp) / sum(pos)
  }
  # Accuracy
  acc <- (sum(tp) + sum(tn)) / sum(pos + neg)
  
  
  # Returns
  return(list(tp  = tp,
              fp  = fp,
              tn  = tn,
              fn  = fn,
              tpr = tpr,
              tnr = tnr,
              fdr = fdr,
              pre = pre,
              rec = rec,
              acc = acc))
}