object<- rval.out[rval.out$delta < 2,]
fit = FALSE
revised.var = TRUE
# if (!missing(subset)) {
#   cl <- match.call()
#   cl[[1L]] <- as.name("subset")
#   names(cl)[2L] <- "x"
#   object <- eval.parent(cl[1L:3L])
# }
# if (fit || !missing(...)) {
#   cl <- match.call()
#   cl$fit <- NULL
#   arg1 <- names(cl)[-(1L:2L)] %in% names(formals("model.avg.default"))
#   cl1 <- cl[c(TRUE, TRUE, !arg1)]
#   cl1[[1L]] <- as.name("get.models")
#   if (is.null(cl1[["subset"]])) 
#     cl1[["subset"]] <- NA
#   cl2 <- cl[c(TRUE, TRUE, arg1)]
#   cl2[[2L]] <- cl1
#   cl2[[1L]] <- as.name("model.avg")
#   return(eval(cl2, parent.frame()))
# }
if (nrow(object) <= 1L) 
  stop("'object' consists of only one model")
ct <- attr(object, "coefTables")
cfarr <- coefArray(ct)
weight <- Weights(object)
cfmat <- cfarr[, 1L, ]
cfmat[is.na(cfmat)] <- 0
coefMat <- array(dim = c(2L, ncol(cfmat)), dimnames = list(c("full", 
                                                             "subset"), colnames(cfmat)))
coefMat[1L, ] <- drop(weight %*% cfmat)
coefMat[2L, ] <- coefMat[1L, ]/colSums(array(weight * as.numeric(!is.na(cfarr[, 
                                                                              1L, ])), dim = dim(cfmat)))
coefMat[is.nan(coefMat)] <- NA_real_
all.terms <- attr(object, "terms")
all.vterms <- all.terms[!(all.terms %in% attr(all.terms, 
                                              "interceptLabel") | apply(is.na(object[, all.terms]), 
                                                                        2L, all))]
allterms1 <- applyrns(!is.na(object[, all.vterms, drop = FALSE]), 
                      function(x) all.vterms[x])
allmodelnames <- .modelNames(allTerms = allterms1, uqTerms = all.vterms)
mstab <- itemByType(object, c("df", "loglik", "ic", "delta", 
                              "weight"))
rownames(mstab) <- allmodelnames
ret <- list(msTable = structure(as.data.frame(mstab), term.codes = attr(allmodelnames, 
                                                                        "variables")), coefficients = coefMat, coefArray = cfarr, 
            importance = importance(object), x = NULL, residuals = NULL, 
            formula = if (!is.null(attr(object, "global"))) formula(attr(object, 
                                                                         "global")) else NULL, call = match.call())
attr(ret, "rank") <- attr(object, "rank")
attr(ret, "modelList") <- attr(object, "modelList")
if (is.null(attr(ret, "modelList"))) 
  attr(ret, "model.calls") <- attr(object, "model.calls")
attr(ret, "beta") <- attr(object, "beta")
attr(ret, "nobs") <- attr(object, "nobs")
attr(ret, "revised.var") <- revised.var
class(ret) <- "averaging"
