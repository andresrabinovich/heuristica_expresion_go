cbh <- function(pvalues, alfa){
  pvalues<-sort(pvalues)
  m <- length(pvalues)
  cbh<-rep(FALSE, m)
  for(i in 1:m){
    cbh[i] <- (pvalues[i] <= i*alfa/m)
  }
  return(cbh)
}
