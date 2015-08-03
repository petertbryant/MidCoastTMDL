options(na.action = NULL)
#global.model <- glm(as.formula(obs.vars[,c('log10_FSS_26Aug14',names(obs.fss2))]),data = obs.vars,na.action = NULL)
global.model <- glmssn(as.formula(obs.vars[,c('log10_FSS_26Aug14',"sum_1095_days","DIS_1YR_PARSA","OWN_FED_PRCA","POP_DARCA")]),ssn.object = ssn1,CorModels = c('Exponential.Euclid','Exponential.taildown'))
#global.model <- test

hasSubset <- 1
ct.args = NULL
trace = FALSE
evaluate = TRUE
betaMode = NULL
#gmCall <- substitute(global.model)
#gmCall <- get_call(global.model)
gmCall <- global.model$args$call
gmEnv <- parent.frame()
gmNobs <- global.model$sampinfo$sample.size
m.min <- 0
m.max <- NA
gmCoefNames <- global.model$sampinfo$effnames
gmFormulaEnv <- environment(as.formula(global.model$args$formula, 
                                       env = gmEnv))
allTerms <- allTerms0 <- getAllTerms(global.model, intercept = TRUE, 
                                     data = eval(gmCall$data, envir = gmEnv))
deps <- attr(allTerms0, "deps")
interceptLabel <- attr(allTerms,"interceptLabel")
nIntercepts <- sum(attr(allTerms, "intercept"))
fixed <- NULL
fixed <- union(fixed, rownames(deps)[rowSums(deps, na.rm = TRUE) == 
                                       ncol(deps)])
fixed <- c(fixed, allTerms[allTerms %in% interceptLabel])
termsOrder <- order(allTerms %in% fixed)
allTerms <- allTerms[termsOrder]
nFixed <- length(fixed)
nVars <- length(allTerms)
nVariants <- 1L
nVarying <- 0L
nov <- as.integer(nVars - nFixed)
#ncomb <- (2L^nov) * nVariants
ncomb <- 3
nmax <- ncomb * nVariants
rvNcol <- nVars + nVarying + 3L
rvChunk <- 25L
rval <- matrix(NA_real_, ncol = rvNcol, nrow = rvChunk)
coefTables <- vector(rvChunk, mode = "list")
comb.sfx <- rep(TRUE, nFixed)
comb.seq <- if (nov != 0L) {
  seq_len(nov)
} else {0L}
k <- 0L
extraResult1 <- integer(0L)
calls <- vector(mode = "list", length = rvChunk)
ord <- integer(rvChunk)
lik <- list(logLik = logLik, name = "logLik")
logLik <- lik$logLik
argsOptions <- list(response = attr(allTerms0, "response"), 
                    intercept = nIntercepts, interceptLabel = interceptLabel, 
                    random = attr(allTerms0, "random"), gmCall = gmCall, 
                    gmEnv = gmEnv, allTerms = allTerms0, gmCoefNames = gmCoefNames, 
                    gmDataHead = if (!is.null(gmCall$data)) {
                      if (eval(call("is.data.frame", gmCall$data), gmEnv)) eval(call("head", 
                                                                                     gmCall$data, 1L), gmEnv) else gmCall$data
                    } else NULL, gmFormulaEnv = gmFormulaEnv)
matchCoefCall <- as.call(c(alist(matchCoef, fit1, all.terms = allTerms, 
                                 allCoef = TRUE), ct.args))
retColIdx <- if (nVarying) {-nVars - seq_len(nVarying)
  } else TRUE
iComb <- -1L
while ((iComb <- iComb + 1L) < ncomb) {
  varComb <- iComb%%nVariants
  jComb <- (iComb - varComb)/nVariants
  if (varComb == 0L) {
    isok <- TRUE
    comb <- c(as.logical(intToBits(jComb)[comb.seq]), 
              comb.sfx)
    nvar <- sum(comb) - nIntercepts
    if (nvar > Inf || nvar < m.min || !formula_margin_check(comb, 
                                                              deps) || switch(hasSubset, FALSE, !all(subset[comb, 
                                                                                                            comb], na.rm = TRUE), !evalExprInEnv(subsetExpr, 
                                                                                                                                                 env = ssEnv, enclos = parent.frame(), comb = comb, 
                                                                                                                                                 `*nvar*` = nvar), FALSE)) {
      isok <- FALSE
      next
    }
    newArgs <- makeArgs(global.model, allTerms[comb], 
                        comb, argsOptions)
    formulaList <- if (is.null(attr(newArgs, "formulaList"))) {newArgs
      } else {
        attr(newArgs, "formulaList")
      }
    if (!is.null(attr(newArgs, "problems"))) {
      print.warnings(structure(vector(mode = "list", 
                                      length = length(attr(newArgs, "problems"))), 
                               names = attr(newArgs, "problems")))
    }
    cl <- gmCall
    cl[names(newArgs)] <- newArgs
    if (!isok) 
      next
    clVariant <- cl
    if (evaluate) {
      start.time <- Sys.time()
      print(start.time)
      fit1 <- tryCatch(eval(clVariant, gmEnv), error = function(err) {
        err$message <- paste(conditionMessage(err), "(model", 
                             iComb, "skipped)", collapse = "")
        class(err) <- c("simpleError", "warning", "condition")
        warning(err)
        return(NULL)
      })
      end.time <- Sys.time()
      print(end.time)
      print(end.time - start.time)
      if (is.null(fit1)) 
        next
      ll1 <- fit1$estimates$m2LL
      nobs1 <- fit1$sampinfo$sample.size
      mcoef1 <- eval(matchCoefCall)
      if (nobs1 != gmNobs) 
        cry(, "number of observations in model #%d (%d) different from global model (%d)", 
            iComb, nobs1, gmNobs, warn = TRUE)
      row1 <- c(mcoef1[allTerms], extraResult1, df = nobs1 - fit1$sampinf$rankX, ll = ll1, ic = AICc.glmssn(fit1))
      k <- k + 1L
      rvlen <- nrow(rval)
      if (retNeedsExtending <- k > rvlen) {
        nadd <- min(rvChunk, nmax - rvlen)
        rval <- rbind(rval, matrix(NA_real_, ncol = rvNcol, 
                                   nrow = nadd), deparse.level = 0L)
        addi <- seq.int(rvlen + 1L, length.out = nadd)
        coefTables[addi] <- vector("list", nadd)
      }
      rval[k, retColIdx] <- row1
      coefTables[[k]] <- attr(mcoef1, "coefTable")
    }
    else {
      k <- k + 1L
      rvlen <- length(ord)
      if (retNeedsExtending <- k > rvlen) {
        nadd <- min(rvChunk, nmax - rvlen)
        addi <- seq.int(rvlen + 1L, length.out = nadd)
      }
    }
    if (retNeedsExtending) {
      calls[addi] <- vector("list", nadd)
      ord[addi] <- integer(nadd)
    }
    ord[k] <- iComb
    calls[[k]] <- clVariant
  }
}  
if (k == 0L) 
  stop("result is empty")
ord <- ord + 1L
names(calls) <- ord
if (!evaluate) 
  return(calls[seq_len(k)])
if (k < nrow(rval)) {
  i <- seq_len(k)
  rval <- rval[i, , drop = FALSE]
  ord <- ord[i]
  calls <- calls[i]
  coefTables <- coefTables[i]
}
if (nVarying) {
  varlev <- ord%%nVariants
  varlev[varlev == 0L] <- nVariants
  rval[, nVars + seq_len(nVarying)] <- variants[varlev, 
                                                ]
}
rval <- as.data.frame(rval)
row.names(rval) <- ord
tfac <- which(!(allTerms %in% gmCoefNames))
rval[tfac] <- lapply(rval[tfac], factor, levels = NaN, labels = "+")
rval[, seq_along(allTerms)] <- rval[, v <- order(termsOrder)]
allTerms <- allTerms[v]
colnames(rval) <- c(allTerms, "df", 
                    lik$name, 'AICc')
rval <- rval[o <- order(rval[, 'AICc'], decreasing = FALSE), 
             ]
coefTables <- coefTables[o]
rval$delta <- rval[, 'AICc'] - min(rval[, 'AICc'])
rval$weight <- exp(-rval$delta/2)/sum(exp(-rval$delta/2))
mode(rval$df) <- "integer"
rval.out <- structure(rval, model.calls = calls[o], global = global.model, 
          global.call = gmCall, terms = structure(allTerms, interceptLabel = interceptLabel), 
          rank = AICc, call = match.call(expand.dots = TRUE), 
          coefTables = coefTables, nobs = gmNobs,  
          column.types = {
            colTypes <- c(terms = length(allTerms),  
                          df = 1, loglik = 1, 
                          ic = 1, delta = 1, weight = 1)
            column.types <- rep(1L:length(colTypes), colTypes)
            names(column.types) <- colnames(rval)
            lv <- 1L:length(colTypes)
            factor(column.types, levels = lv, labels = names(colTypes)[lv])
          }, class = c("model.selection", "data.frame"))
