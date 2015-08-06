library(MuMIn)
library(SSN)
library(parallel)
library(snow)
options(na.action = NULL)
#global.model <- glm(as.formula(obs.vars[,c('log10_FSS_26Aug14',names(obs.fss2))]),data = obs.vars,na.action = NULL)
#global.model <- glmssn(as.formula(obs.vars[,c('log10_FSS_26Aug14',"sum_1095_days","DIS_1YR_PARSA","OWN_FED_PRCA","POP_DARCA")]),ssn.object = ssn1,CorModels = c('Exponential.Euclid','Exponential.taildown'))
#save(file='ssn1_testing.Rdata',global.model)
load('ssn1_testing.Rdata')
#global.model <- test

#Cluster testing
cluster <- makeCluster(7)

source('funMuMInhelpers.R')

nextra <- NULL
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
m.max <- Inf
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
qi <- 0L
qlen <- 25L
queued <- vector(qlen, mode = "list")
props <- list(gmEnv = gmEnv, IC = 'AICc.glmssn', matchCoefCall = as.call(c(list(as.name("matchCoef"), 
                                                                                      as.name("fit1"), all.terms = allTerms, beta = betaMode, 
                                                                                      allCoef = TRUE), ct.args)))
clusterVExport(cluster, pdredge_props = props, .pdredge_process_model = .pdredge_process_model)
clusterCall(cluster, eval, call("options", options("na.action")), 
            env = 0L)
clusterExport(cluster, lsf.str())
clusterCall(cluster, function() {library(SSN); library(MuMIn)})
clusterExport(cluster, 'ssn1')
clusterExport(cluster, 'allTerms')
matchCoefCall <- as.call(c(alist(matchCoef, fit1, all.terms = allTerms, 
                                 allCoef = TRUE), ct.args))
retColIdx <- if (nVarying) {-nVars - seq_len(nVarying)
  } else TRUE
warningList <- list()
iComb <- -1L
while ((iComb <- iComb + 1L) < ncomb) {
  varComb <- iComb%%nVariants
  jComb <- (iComb - varComb)/nVariants
  if (varComb == 0L) {
    isok <- TRUE
    comb <- c(as.logical(intToBits(jComb)[comb.seq]), 
              comb.sfx)
    nvar <- sum(comb) - nIntercepts
    if ((nvar >= m.min && nvar <= m.max) && formula_margin_check(comb, 
         deps) && switch(hasSubset, TRUE, all(subset[comb, 
         comb], na.rm = TRUE), evalExprInEnv(subsetExpr, 
         env = ssEnv, enclos = parent.frame(), comb = comb, 
         `*nvar*` = nvar), TRUE)) {
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
    } else {
      isok <- FALSE
    }
  }
  if (isok) {
    clVariant <- cl
    if (evaluate) {
      qi <- qi + 1L
      queued[[(qi)]] <- list(call = clVariant, id = iComb)
    }
    else {
      k <- k + 1L
      rvlen <- length(ord)
      if (k > rvlen) {
        nadd <- min(rvChunk, nmax - rvlen)
        addi <- seq.int(rvlen + 1L, length.out = nadd)
        calls[addi] <- vector("list", nadd)
        ord[addi] <- integer(nadd)
      }
      calls[[k]] <- clVariant
      ord[k] <- iComb
    }
  }
  if (evaluate && qi && (qi > qlen || (iComb + 1L) == ncomb)) {
    qseq <- seq_len(qi)
    qresult <- .getRow(queued[qseq])
    utils::flush.console()
    if (any(vapply(qresult, is.null, TRUE))) 
      stop("some results returned from cluster node(s) are NULL. \n", 
           "This should not happen and indicates problems with ", 
           "the cluster node", domain = "R-MuMIn")
    haveProblems <- logical(qi)
    nadd <- sum(sapply(qresult, function(x) inherits(x$value, 
                                                     "condition") + length(x$warnings)))
    wi <- length(warningList)
    if (nadd) 
      warningList <- c(warningList, vector(nadd, mode = "list"))
    for (i in qseq) for (cond in c(qresult[[i]]$warnings, 
                                   if (inherits(qresult[[i]]$value, "condition")) list(qresult[[i]]$value))) {
      wi <- wi + 1L
      warningList[[wi]] <- if (is.null(conditionCall(cond))) 
        queued[[i]]$call
      else conditionCall(cond)
      if (inherits(cond, "error")) {
        haveProblems[i] <- TRUE
        msgsfx <- "(model %d skipped)"
      }
      else msgsfx <- "(in model %d)"
      names(warningList)[wi] <- paste(conditionMessage(cond), 
                                      gettextf(msgsfx, queued[[i]]$id))
      attr(warningList[[wi]], "id") <- queued[[i]]$id
    }
    withoutProblems <- which(!haveProblems)
    qrows <- lapply(qresult[withoutProblems], "[[", "value")
    qresultLen <- length(qrows)
    rvlen <- nrow(rval)
    if (retNeedsExtending <- k + qresultLen > rvlen) {
      nadd <- min(max(rvChunk, qresultLen), nmax - 
                    rvlen)
      rval <- rbind(rval, matrix(NA_real_, ncol = rvNcol, 
                                 nrow = nadd), deparse.level = 0L)
      addi <- seq.int(rvlen + 1L, length.out = nadd)
      coefTables[addi] <- vector("list", nadd)
      calls[addi] <- vector("list", nadd)
      ord[addi] <- integer(nadd)
    }
    qseqOK <- seq_len(qresultLen)
    for (m in qseqOK) rval[k + m, retColIdx] <- qrows[[m]]
    ord[k + qseqOK] <- vapply(queued[withoutProblems], 
                              "[[", 1L, "id")
    calls[k + qseqOK] <- lapply(queued[withoutProblems], 
                                "[[", "call")
    coefTables[k + qseqOK] <- lapply(qresult[withoutProblems], 
                                     "[[", "coefTable")
    k <- k + qresultLen
    qi <- 0L
  }
  }

if (k == 0L) 
  stop("result is empty")
#ord <- ord + 1L
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
  rval[, nVars + seq_len(nVarying)] <- variants[varlev,]
}

rval <- as.data.frame(rval)
row.names(rval) <- ord
tfac <- which(!(allTerms %in% gmCoefNames))
rval[tfac] <- lapply(rval[tfac], factor, levels = NaN, labels = "+")
rval[, seq_along(allTerms)] <- rval[, v <- order(termsOrder)]
allTerms <- allTerms[v]
colnames(rval) <- c(allTerms, "df", lik$name, 'AICc')
rval <- rval[o <- order(rval[, 'AICc'], decreasing = FALSE),]
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

stopCluster(cluster)
