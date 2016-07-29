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