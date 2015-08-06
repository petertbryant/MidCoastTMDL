library(MuMIn)
library(SSN)

`coefArray` <- function(object) {
  coefNames <- fixCoefNames(unique(unlist(lapply(object, rownames),
                                          use.names = FALSE)))
  nCoef <- length(coefNames)
  nModels <- length(object)
  rval <- array(NA_real_, dim = c(nModels, 3L, nCoef),
                dimnames = list(names(object), c("Estimate", "Std. Error", "df"), coefNames))
  for(i in seq_along(object)) {
    z <- object[[i]]
    rval[i, seq_len(ncol(z)), ] <- t(z[match(coefNames, fixCoefNames(rownames(z))), ])
  }
  rval
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

formula_margin_check <- function(j, m) {
  stopifnot(is.logical(j))
  !any(m[!j[-length(j)], j-length(j)], na.rm = TRUE)
}

makeArgs <- function(obj, termNames, comb, opt, ...) {
  reportProblems <- character(0L)
  termNames[termNames %in% opt$interceptLabel] <- "1"
  ## XXX: what if 'opt$intercept' is of length > 1 ???
  f <- reformulate(c(if(!opt$intercept) "0", termNames), response = opt$response)
  environment(f) <- opt$gmFormulaEnv
  ret <- list(formula = f)
  if(!is.null(opt$gmCall$start)) {
    coefNames <- fixCoefNames(.getCoefNames(f, opt$gmDataHead,
                                            opt$gmCall$contrasts, envir = opt$gmEnv))
    idx <- match(coefNames, opt$gmCoefNames)
    if(anyNA(idx)) reportProblems <-
      append(reportProblems, "cannot subset 'start' argument. Coefficients in model do not exist in global.model")
    else ret$start <- substitute(start[idx], list(start = opt$gmCall$start,
                                                  idx = idx))
  }
  attr(ret, "formulaList") <- list(f)
  attr(ret, "problems") <- reportProblems
  ret
}

matchCoef <- function(m1, m2,
                      all.terms = getAllTerms(m2, intercept = TRUE),
                      beta = 0L,
                      terms1 = getAllTerms(m1, intercept = TRUE),
                      coef1 = NULL,
                      allCoef = TRUE,
                      ...
) {
  if(is.null(coef1)) {
    ct <- if (beta != 0L) std.coef(m1, beta == 2L, ...) else coefTable(m1, ...)
    ct.df <- summary(m1)$fixed.effects.estimates[,1:3]
    ct <- as.matrix(ct.df[,2:3])
    rownames(ct) <- ct.df[,1]
    coef1 <- ct[, 1L]
    names(coef1) <- rownames(ct)
  } else if(allCoef) stop("'coef1' is given and 'allCoef' is not FALSE")
  
  if(any((terms1 %in% all.terms) == FALSE)) stop("'m1' is not nested within 'm2'")
  row <- structure(rep(NA_real_, length(all.terms)), names = all.terms)
  
  fxdCoefNames <- fixCoefNames(names(coef1))
  row[terms1] <- NaN
  pos <- match(terms1, fxdCoefNames, nomatch = 0L)
  row[fxdCoefNames[pos]] <- coef1[pos]
  if(allCoef) {
    i <- match(names(coef1), rownames(ct))
    j <- !is.na(i)
    rownames(ct)[i[j]] <- fxdCoefNames[j]
    attr(row, "coefTable") <- ct
  }
  row
}

fixCoefNames <-
  function(x, peel = TRUE) {
    if(!length(x)) return(x)
    ox <- x
    ia <- grep(":", x, fixed = TRUE)
    if(!length(ia)) return(structure(x, order = rep.int(1L, length(x))))
    x <- ret <- x[ia]
    if(peel) {
      # case of pscl::hurdle, cf are prefixed with count_|zero_
      if(all(substr(x, 1L, pos <- regexpr("_", x, fixed = TRUE)) %in%
             c("count_", "zero_"))) {
        ret <- substr(ret, pos + 1L, 256L)
        k <- TRUE
        suffix <- ""
      } else { # unmarkedFit with its phi(...), lambda(...) etc...
        k <- grepl("^\\w+\\(.+\\)$", x, perl = TRUE)
        fname <- substring(x[k], 1L, attr(regexpr("^\\w+(?=\\()", x[k],
                                                  perl = TRUE),"match.length"))
        
        # do not peel off if a function
        k[k] <- !vapply(fname, exists, FALSE, mode = "function", envir = .GlobalEnv)
        if(any(k)) {
          pos <- vapply(x[k], function(z) {
            parens <- lapply(lapply(c("(", ")"),
                                    function(s) gregexpr(s, z, fixed = TRUE)[[1L]]),
                             function(y) y[y > 0L])
            parseq <- unlist(parens, use.names = FALSE)
            p <- cumsum(rep(c(1L, -1L), sapply(parens, length))[order(parseq)])
            if(any(p[-length(p)] == 0L)) -1L else parseq[1L]
          }, 1L, USE.NAMES = FALSE)
          k[k] <- pos != -1L
          pos <- pos[pos != -1]
          if(any(k)) ret[k] <- substring(x[k], pos + 1L, nchar(x[k]) - 1L)
        }
        suffix <- ")"
      }
    } else	k <- FALSE
    
    ## prepare = replace multiple ':' to avoid splitting by '::' and ':::'
    spl <- expr.split(ret, ":", prepare = function(x) gsub("((?<=:):|:(?=:))", "_", x, perl = TRUE))
    ret <- vapply(lapply(spl, base::sort), paste0, "", collapse = ":")
    if(peel && any(k))
      ret[k] <- paste0(substring(x[k], 1L, pos), ret[k], suffix)
    ox[ia] <- ret
    ord <- rep.int(1, length(ox))
    ord[ia] <- sapply(spl, length)
    structure(ox, order = ord)
  }

formula.glmssn <- function(x, ...) {
  form <- x$args$formula
  environment(form) <- environment(form)
  form
}

# like apply(, 2) but returns a list (does not do any checking)
`applyrns` <- function (X, FUN, ...) {
  n <- nrow(X)
  ret <- vector(n, mode = "list")
  for(i in seq_len(n)) if(!is.null(z <- FUN(X[i, ], ...))) ret[[i]] <- z
  ret
}

type2col <-
  function (x, type) {
    if (inherits(x, "model.selection")) 
      x <- attr(x, "column.types")
    k <- match(x, type, nomatch = 0L)
    i <- k != 0
    which(i)[order(k[i])]
  }

`itemByType` <- function(x, type, i, ...) 
  `[.data.frame`(x, i, type2col(x, type), ...)

`itemByType<-` <- function(x, type, i, value)
  `[<-.data.frame`(x, i, type2col(x, type), value)

`.pdredge_process_model` <- function(modv, envir = get("pdredge_props", .GlobalEnv)) {
  ### modv == list(call = clVariant, id = modelId)
  result <- tryCatchWE(eval(modv$call, get("gmEnv", envir)))
  if (inherits(result$value, "condition")) return(result)
  
  fit1 <- result$value
#   if(get("nextra", envir) != 0L) {
#     extraResult1 <- get("applyExtras", envir)(fit1)
#     nextra <- get("nextra", envir)
#     if(length(extraResult1) < nextra) {
#       tmp <- rep(NA_real_, nextra)
#       tmp[match(names(extraResult1), get("extraResultNames", envir))] <-
#         extraResult1
#       extraResult1 <- tmp
#     }
#   } else 
    extraResult1 <- NULL
  #ll <- .getLik(fit1)$logLik(fit1)
  ll <- fit1$estimates$m2LL
  aic <- AICc.glmssn(fit1)
  mcoef <- matchCoef(fit1, all.terms = allTerms)
  # beta = get("beta", envir), allCoef = TRUE)
  #mcoef <- eval(get("matchCoefCall", envir))
  
  
  list(value = c(mcoef, extraResult1, df = fit1$sampinfo$sample.size - fit1$sampinfo$rankX, ll = ll,
                 ic = aic),
       nobs = fit1$sampinfo$sample.size,
       coefTable = attr(mcoef, "coefTable"),
       warnings = result$warnings)
  
}

`tryCatchWE` <- function (expr) {
  Warnings <- NULL
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                   warning = function(w) {
                                     Warnings <<- c(Warnings, list(w))
                                     invokeRestart("muffleWarning")
                                   }), warnings = Warnings)
}

`clusterVExport` <- local({
  
  `getv` <- function(obj, env = as.environment(1L))
    for (i in names(obj)) assign(i, obj[[i]], envir = env)
  function(cluster, ...) {
    Call <- match.call()
    Call$cluster <- NULL
    Call <- Call[-1L]
    vars <- list(...)
    vnames <- names(vars)
    #if(!all(sapply(Call, is.name))) warning("at least some elements do not have syntactic name")
    if(is.null(vnames)) {
      names(vars) <- vapply(Call, asChar, "")
    } else if (any(vnames == "")) {
      names(vars) <- ifelse(vnames == "", vapply(Call, asChar, ""), vnames)
    }
    get("clusterCall")(cluster, getv, vars)
    # clusterCall(cluster, getv, vars)
  }
})

.getRow <- function(X) clusterApply(cluster, X, fun = ".pdredge_process_model")
