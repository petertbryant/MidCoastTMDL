lrtSSN <- function(object1, object2) {
  df1.n <-  object1$sampinf$obs.sample.size
  df1 <- df1.n - object1$sampinf$rankX
  
  df2.n <-  object2$sampinf$obs.sample.size
  df2 <- df2.n - object2$sampinf$rankX
  
  lrt <- object2$estimates$m2LL - object1$estimates$m2LL
  diff.df <- df2 - df1
  if (lrt < 0) {
    lrt <- -lrt
    diff.df <- -diff.df
  }
  if (lrt * diff.df < 0) {
    stop("Likelihood gets worse with more variables. Test not executed")
  }

output <- list(Chisquared = lrt, df = diff.df, 
               p.value = pchisq(lrt, diff.df, lower.tail = FALSE))
return(output)
}