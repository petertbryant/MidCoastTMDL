InfoCritCompare2 <- function (model.list) 
{
  IC <- NULL
  for (i in 1:length(model.list)) {
    if (class(model.list[[i]]) != "glmssn") {
      stop("All models must be of type glmssn")
    }
    model.name <- NULL
    ind <- !duplicated(attributes(model.list[[i]]$estimates$theta)$terms)
    terms <- attributes(model.list[[i]]$estimates$theta)$terms[ind]
    model.name <- paste(terms, collapse = " + ")
    if (model.list[[i]]$args$family != "gaussian") {
      model.AIC <- NA
      model.neg2LogL <- NA
    }
    if (model.list[[i]]$args$family == "gaussian") {
      model.AIC <- AIC(model.list[[i]])
      model.neg2LogL <- model.list[[i]]$estimates$m2LL
    }
    IC <- rbind(IC, data.frame(formula = deparse(model.list[[i]]$args$formula, 
                                                 width.cutoff = 500), 
                               EstMethod = model.list[[i]]$args$EstMeth, 
                               Variance_Components = model.name, 
                               neg2LogL = model.neg2LogL, 
                               AIC = model.AIC, 
                               CrossValidationStatsSSN_ptb(model.list[[i]])))
  }
  IC
}

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

CrossValidationStatsSSN_ptb <- function (object) {
  cv.out <- kfold_CrossValidationSSN(object)
  zt <- object$sampinf$z
  n <- object$sampinf$obs.sample.size
  dfobs <- n - object$sampinf$rankX
  cv.stats <- data.frame(bias = sum(cv.out[, "cv.pred"] - zt)/n, 
                         std.bias = sum((cv.out[, "cv.pred"] - zt)/sqrt(cv.out[, 
                                                                               "cv.se"]))/n, RMSPE = sqrt(sum((zt - cv.out[, "cv.pred"])^2)/n), 
                         RAV = sqrt(sum(cv.out[, "cv.se"]^2)/n), std.MSPE = sum((cv.out[, 
                                                                                        "cv.pred"] - zt)^2/cv.out[, "cv.se"]^2)/n, cov.80 = sum(abs((zt - 
                                                                                                                                                       cv.out[, "cv.pred"])/cv.out[, "cv.se"]) < qt(0.9, 
                                                                                                                                                                                                    dfobs))/n, cov.90 = sum(abs((zt - cv.out[, "cv.pred"])/cv.out[, 
                                                                                                                                                                                                                                                                  "cv.se"]) < qt(0.95, dfobs))/n, cov.95 = sum(abs((zt - 
                                                                                                                                                                                                                                                                                                                      cv.out[, "cv.pred"])/cv.out[, "cv.se"]) < qt(0.975, 
                                                                                                                                                                                                                                                                                                                                                                   dfobs))/n)
  cv.stats
}


format.perc <- function(probs, digits)
  ## Not yet exported, maybe useful in other contexts:
  ## quantile.default() sometimes uses a version of it
  paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits),
        "%")

confint.glmssn <- function (object, parm, level = 0.95, ...) 
{
  cf <- t(object$estimates$betahat)
  names(cf) <- colnames(cf)
  pnames <- names(cf)
  if (missing(parm)) 
    parm <- pnames
  else if (is.numeric(parm)) 
    parm <- pnames[parm]
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  pct <- format.perc(a, 3)
  fac <- qnorm(a)
  ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, 
                                                             pct))
  ses <- sqrt(diag(object$estimates$covb))[parm]
  ci[] <- cf[parm] + ses %o% fac
  ci
}

predict.vary <- function(betahat, ss, r_vec) {
  options(warn = -1)
  betahat <- plyr::rename(betahat, c('HDWTR1' = 'HDWTR'))
  bsubm <- melt(betahat)
  r_vec <- plyr::rename(r_vec, c('HDWTR1' = 'HDWTR'))
  r_vec <- melt(r_vec)
  r_vec$variable <- rownames(r_vec)
  r_vec <- plyr::rename(r_vec, c("value" = "value.betahat.vary"))
  for (i in 1:nrow(ss)) {
    vals <- ss[i, names(betahat)[-1]]
    vals <- melt(vals, measure.vars = 1:length(names(betahat[-1])))
    vals$value <- as.numeric(vals$value)
    inter <- merge(bsubm, vals, by = 'variable', suffixes = c('.betahat',''))
    inter$inter <- inter$value.betahat * inter$value
    BSTI <- ss[i, 'log10_BSTI']
    Z <- BSTI - betahat$`(Intercept)`[1] - sum(inter$inter)
    
    inter <- merge(inter, r_vec, by = 'variable')
    inter$inter3 <- inter$value * inter$value.betahat.vary
    BSTI_prd <- Z + betahat$`(Intercept)` + sum(inter$inter3)
    
    if (i == 1) {
      BSTI_prd_df <- data.frame("pid" = ss[i, 'pid'], "BSTI_prd" = BSTI_prd)
    } else {
      BSTI_prd_df <- rbind(BSTI_prd_df, data.frame("pid" = ss[i, 'pid'], "BSTI_prd" = BSTI_prd))
    }
  }
  
  return(BSTI_prd_df)
  options(warn = 0)
}

simplify_target_equation <- function(betahat, ss, station) {
  options(warn = -1)
  betahat <- plyr::rename(betahat, c('HDWTR1' = 'HDWTR'))
  bsubm <- melt(betahat)
  vals <- ss[ss$STATION_KEY == station, names(betahat)[-1]]
  vals <- melt(vals, measure.vars = 1:length(names(betahat[-1])))
  vals$value <- as.numeric(vals$value)
  inter <- merge(bsubm, vals, by = 'variable', suffixes = c('.betahat',''))
  inter$inter <- inter$value.betahat * inter$value
  BSTI <- ss[ss$STATION_KEY == station, 'log10_pred']
  Z <- BSTI - betahat$`(Intercept)`[1] - sum(inter$inter)
  
  inter2 <- inter[inter$variable != 'sum_1095_days',]
  b <- betahat$`(Intercept)`[1] + sum(inter2$inter) + Z
  m <- betahat$sum_1095_days
  
  
  x <- list("b" = b, "m" = m, "Z" = Z)
  attr(x, 'inter') <- inter
  attr(x, 'inter2') <- inter2
  
  return(x)
  options(warn = 0)
}

simplify_target_equation_all <- function(betahat, ss, var_means, var_sd) {
  options(warn = -1)
  betahat <- plyr::rename(betahat, c('HDWTR1.21463119808789' = 'HDWTR'))
  bsubm <- melt(betahat)
  vals <- ss[, names(betahat)[-1]]
  vals <- melt(vals, measure.vars = 1:length(names(betahat[-1])))
  vals$value <- as.numeric(vals$value)
  inter <- merge(bsubm, vals, by = 'variable', suffixes = c('.betahat',''))
  inter$inter <- inter$value.betahat * inter$value
  BSTI <- ss[, 'log10_BSTI']
  Z <- BSTI - betahat$`(Intercept)`[1] - sum(inter$inter)
  
  inter2 <- inter[inter$variable != 'sum_1095_days',]
  b <- betahat$`(Intercept)`[1] + sum(inter2$inter) + Z
  b_u <- (b*(2*var_sd["BSTI"]))+var_means["BSTI"]
  m_u <- (betahat$sum_1095_days*(2*var_sd["BSTI"]))+var_means["BSTI"]
  
  
  x <- list("b_u" = b_u, "m_u" = m_u, "Z" = Z, "SVN" = ss[,'SVN'])
  attr(x, 'inter') <- inter
  attr(x, 'inter2') <- inter2
  
  return(x)
  options(warn = 0)
}

vif_func<-function(in_frame,ob.ssn){
  
  vif_init<-NULL
  for(val in names(in_frame)){
    form_in<-formula(paste(val,' ~ ',paste(setdiff(names(in_frame),val),collapse = ' + ')))
    tmp.ssn <- glmssn(form_in, 
                      ob.ssn,
                      EstMeth = "ML",
                      CorModels = c("locID", "Exponential.taildown", "Exponential.Euclid"),
                      addfunccol = "afvArea",
                      family = "Gaussian")
    ssn.summary <- summary(tmp.ssn)
    VIF.val <- 1/(1 - ssn.summary$Rsquared)
    new.row <- c(val, VIF.val)
    vif_init <- rbind(vif_init,new.row)
  }
  return(vif_init)
}

#### Correlation Matrix Functions ####
# source
# http://stackoverflow.com/questions/15271103/how-to-modify-this-correlation-matrix-plot
#
panel.cor <- function(x, y, digits=2, cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  test <- cor.test(x,y)
  Signif <- ifelse(round(test$p.value,3)<0.001,"p<0.001",paste("p=",round(test$p.value,3)))  
  text(0.5, 0.25, paste("r=",txt))
  text(.5, .75, Signif)
}

panel.smooth<-function (x, y, col = "black", bg = NA, pch = 18, 
                        cex = 0.8, col.smooth = "red", span = 2/3, iter = 3, ...) 
{
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
          col = col.smooth, ...)
}

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="gray", ...)
}

AICc.glmssn <- function (object, ..., k = 2) 
{
  if (class(object) != "glmssn") 
    return("Not a glmssn object")
  if (!missing(k)) 
    cat("This argument has no effect\n")
  if (object$args$family == "poisson") 
    return(NA)
  if (object$args$family == "binomial") 
    return(NA)
  cp <- covparms(object)[, 3]
  nparmstheta <- length(cp)
  rankX <- object$sampinfo$rankX
  if (object$args$EstMeth == "REML") 
    nparms <- nparmstheta
  if (object$args$EstMeth == "ML") 
    nparms <- rankX + nparmstheta
  object$estimates$m2LL + (k*nparms*(nparms+1))/(object$sampinfo$sample.size-nparms-1)
}

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
