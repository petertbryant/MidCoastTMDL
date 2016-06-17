library(SSN)

kfold_CrossValidationSSN <- function(fit) {
  z <- fit$sampinfo$z
  X <- fit$sampinfo$X
  V <- fit$estimates$V
  Vi <- fit$estimates$Vi
  n <- fit$sampinfo$obs.sample.size
  #These are in here to pull out the locations
  odf <- fit$ssn.object@obspoints@SSNPoints[[1]]@point.data
  stns <- levels(odf$locID)
  nstn <- length(levels(odf$locID))
  #But if you wanted to specify the size of the k fold you could do it here
  #k <- 50
  #stns <- split(sample(1:n), ceiling(seq_along(sample(1:n))/k))
  #nstn <- length(stns)
  cdd.out <- matrix(-999.9, nrow = n, ncol = 3)
  cdd.out[,1] <- attributes(fit$sampinfo$z)$pid
  
  for(i in 1:nstn) {
    #Define j using the observations at each location
    j <- which(rownames(X) %in% odf[odf$locID == stns[i],'pid'])
    #Using the above code to set the k-fold you can grab those indices here
    #j <- stns[[i]]
    
    if (length(j) == 1) {
      Vi.i <- Vi[(1:n) != j,(1:n) != j] -
        matrix(Vi[(1:n) != j,j],ncol = 1) %*%
        matrix(Vi[j,(1:n) != j],nrow = 1)/Vi[j,j]
      c.i <- matrix(V[(1:n) != j,j],ncol = 1)
      xi <- matrix(X[j,], ncol = 1)
      X.i <- X[(1:n) != j,]
      z.i <- matrix(z[(1:n) != j], ncol = 1)
      xxi <- xi - t(X.i) %*% Vi.i %*% c.i
      covb.i <- solve(t(X.i) %*% Vi.i %*% X.i)
      si <- V[i,i]  - t(c.i) %*% Vi.i %*% c.i
      lam <- t(c.i + X.i %*% covb.i %*% xxi) %*% Vi.i
    } else {
      Vi.i <- Vi[!(1:n) %in% j,!(1:n) %in% j] -
        matrix(Vi[!(1:n) %in% j,j],ncol = length(j)) %*%
        solve(Vi[j,j]) %*%
        matrix(Vi[j,!(1:n) %in% j],nrow = length(j))
      c.i <- V[!(1:n) %in% j,(1:n) %in% j, drop=F]
      xi <- t(X[(1:n) %in% j,])
      X.i <- X[!(1:n) %in% j,]
      z.i <- matrix(z[!(1:n) %in% j], ncol = 1)
      xxi <- xi - t(X.i) %*% Vi.i %*% c.i
      covb.i <- solve(t(X.i) %*% Vi.i %*% X.i)
      si <- V[i,i]  - t(c.i) %*% Vi.i %*% c.i
      lam <- t(c.i + X.i %*% covb.i %*% xxi) %*% Vi.i
    }
    
    cdd.out[j,2] <- lam %*% z.i
    cdd.out[j,3] <- sqrt(diag(si + t(xxi) %*% covb.i %*% xxi))
    
  }
  
  cdd.out <- as.data.frame(cdd.out)
  names(cdd.out) <- c("pid","cv.pred","cv.se")
  cdd.out
  cv.out <- cdd.out
  cv.out
}

# 
#   zt <- fit$sampinf$z
#   n <-  fit$sampinf$obs.sample.size
#   dfobs <- n - fit$sampinf$rankX
#   cv.stats <- data.frame(bias = sum(cv.out[,"cv.pred"] - zt)/n,
#                          std.bias = sum((cv.out[,"cv.pred"] - zt)/sqrt(cv.out[,"cv.se"]))/n,
#                          RMSPE = sqrt(sum((zt - cv.out[,"cv.pred"])^2)/n),
#                          RAV = sqrt(sum(cv.out[,"cv.se"]^2)/n),
#                          std.MSPE = sum((cv.out[,"cv.pred"] - zt)^2/cv.out[,"cv.se"]^2)/n,
#                          cov.80 = sum(abs((zt - cv.out[,"cv.pred"])/cv.out[,"cv.se"]) < qt(.90, dfobs))/n,
#                          cov.90 = sum(abs((zt - cv.out[,"cv.pred"])/cv.out[,"cv.se"]) < qt(.95, dfobs))/n,
#                          cov.95 = sum(abs((zt - cv.out[,"cv.pred"])/cv.out[,"cv.se"]) < qt(.975, dfobs))/n)
#   cv.stats
#   
# }
