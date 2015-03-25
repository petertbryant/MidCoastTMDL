predict.PTB <- function (Xs, pid, target, object, predpointsID, ...) 
{
  if (length(object$estimates$Warnlog) > 0 && length(grep("Algorithm diverging", 
                                                          object$estimates$Warnlog)) > 0) 
    stop("No predictions for diverging algorithm")
  datao <- object$ssn.object@obspoints@SSNPoints[[1]]@point.data
  ocoord <- object$ssn.object@obspoints@SSNPoints[[1]]@point.coords
  theta <- object$estimates$theta
  ind <- object$sampinfo$ind.obs
  nobs <- length(ind)
  CorModels <- object$args$CorModels
  useTailDownWeight <- object$args$useTailDownWeight
  REs <- object$sampinfo$REs
  distord <- order(as.integer(as.character(datao[, "netID"])), 
                   datao[, "pid"])
  a.mat <- NULL
  b.mat <- NULL
  net.zero <- NULL
  w.matrix <- NULL
  dist.hydro <- NULL
  xcoord <- NULL
  ycoord <- NULL
  xyobs <- NULL
  xypred <- NULL
  rnames <- NULL
  if (predpointsID == "_MissingObs_") {
    if (length(grep("tail", CorModels)) > 0 | useTailDownWeight == 
          TRUE) {
      if (length(grep("taildown", CorModels)) > 1) 
        stop("Cannot have more than 1 tailup model")
      n.all <- length(datao[, 1])
      dist.junc <- matrix(0, nrow = nobs, ncol = nobs)
      net.zero <- matrix(0, nrow = nobs, ncol = nobs)
      nsofar <- 0
      nIDs <- sort(as.integer(as.character(unique(datao[, 
                                                        "netID"]))))
      for (i in nIDs) {
        workspace.name <- paste("dist.net", i, ".RData", 
                                sep = "")
        path <- file.path(object$ssn.object@path, "distance", 
                          "obs", workspace.name)
        file_handle <- file(path, "rb")
        distmat <- unserialize(file_handle)
        close(file_handle)
        oi <- order(as.numeric(rownames(distmat)))
        distmat <- distmat[oi, oi, drop = F]
        rnames <- c(rnames, rownames(distmat))
        ni <- length(distmat[1, ])
        dist.junc[(nsofar + 1):(nsofar + ni), (nsofar + 
                                                 1):(nsofar + ni)] <- distmat
        net.zero[(nsofar + 1):(nsofar + ni), (nsofar + 
                                                1):(nsofar + ni)] <- 1
        nsofar <- nsofar + ni
      }
      rownames(dist.junc) <- rnames
      a.mat <- pmax(dist.junc, t(dist.junc))
      b.mat <- pmin(dist.junc, t(dist.junc))
      dist.hydro <- as.matrix(dist.junc + t(dist.junc)) * 
        net.zero
      flow.con.mat <- 1 - (b.mat > 0) * 1
      addfunccol <- object$args$addfunccol
      w.matrix <- sqrt(pmin(outer(datao[distord, addfunccol], 
                                  rep(1, times = nobs)), t(outer(datao[distord, 
                                                                       addfunccol], rep(1, times = nobs))))/pmax(outer(datao[distord, 
                                                                                                                             addfunccol], rep(1, times = nobs)), t(outer(datao[distord, 
                                                                                                                                                                               addfunccol], rep(1, times = nobs))))) * flow.con.mat * 
        net.zero
      if (length(grep("tailup", CorModels)) > 0) {
        if (length(grep("tailup", CorModels)) > 1) 
          stop("Cannot have more than 1 tailup model")
        flow.con.mat <- 1 - (b.mat > 0) * 1
      }
    }
    xcoord <- ocoord[distord, 1, drop = F]
    ycoord <- ocoord[distord, 2, drop = F]
    REind <- which(names(datao) %in% CorModels)
    if (length(REind)) {
      REs <- list()
      REnames <- sort(names(datao)[REind])
      for (ii in 1:length(REind)) REs[[REnames[ii]]] <- model.matrix(~datao[distord, 
                                                                            REnames[ii]] - 1)
      rownames(REs[[REnames[ii]]]) <- datao[distord, "pid"]
      for (ii in 1:length(REind)) REs[[ii]] <- REs[[ii]] %*% 
        t(REs[[ii]])
    }
    Vpred <- makeCovMat(theta = theta, dist.hydro = dist.hydro, 
                        a.mat = a.mat, b.mat = b.mat, w.matrix = w.matrix, 
                        net.zero = net.zero, x.row = xcoord, y.row = ycoord, 
                        x.col = xcoord, y.col = ycoord, CorModels, useTailDownWeight = useTailDownWeight, 
                        use.nugget = object$args$use.nugget, use.anisotropy = FALSE, 
                        REs)
    for (i in 1:length(object$ssn.object@predpoints@SSNPoints)) if (object$ssn.object@predpoints@ID[i] == 
                                                                      "_MissingObs_") {
      datap <- object$ssn.object@predpoints@SSNPoints[[i]]@point.data
      distordp <- order(as.integer(as.character(datap[, 
                                                      "netID"])), datap[, "pid"])
      datap <- datap[distordp, , drop = F]
      pcoord <- object$ssn.object@predpoints@SSNPoints[[i]]@point.coords
      pcoord <- pcoord[distordp, , drop = F]
      netpcoord <- object$ssn.object@predpoints@SSNPoints[[i]]@network.point.coords
      netpcoord <- netpcoord[distordp, , drop = F]
    }
    Vpred <- Vpred[!(datao[distord, "pid"] %in% datap[, "pid"]), 
                   datao[distord, "pid"] %in% datap[, "pid"], drop = F]
  }
  if (predpointsID != "_MissingObs_") {
    for (i in 1:length(object$ssn.object@predpoints@SSNPoints)) if (object$ssn.object@predpoints@ID[i] == 
                                                                      predpointsID) {
      datap <- object$ssn.object@predpoints@SSNPoints[[i]]@point.data
      pcoord <- object$ssn.object@predpoints@SSNPoints[[i]]@point.coords
      netpcoord <- object$ssn.object@predpoints@SSNPoints[[i]]@network.point.coords
    }
    npred <- length(datap[, 1])
    distordp <- order(as.integer(as.character(datap[, "netID"])), 
                      datap[, "pid"])
    datap <- datap[distordp, , drop = F]
    pcoord <- pcoord[distordp, , drop = F]
    netpcoord <- netpcoord[distordp, , drop = F]
    rnames <- NULL
    cnames <- NULL
    if (length(grep("tail", CorModels)) > 0) {
      if (length(grep("taildown", CorModels)) > 1) 
        stop("Cannot have more than 1 tailup model")
      dist.junc.a <- matrix(0, nrow = nobs, ncol = npred)
      dist.junc.b <- matrix(0, nrow = npred, ncol = nobs)
      net.zero <- matrix(0, nrow = nobs, ncol = npred)
      nSoFari <- 0
      nSoFarj <- 0
      nIDs <- sort(unique(c(as.integer(as.character(datao[, 
                                                          "netID"])), as.integer(as.character(datap[, "netID"])))))
      for (i in nIDs) {
        ind.obs <- as.numeric(as.character(datao$netID)) == 
          as.numeric(i)
        site.no <- nrow(datao[ind.obs, ])
        ind.pred <- as.numeric(as.character(datap$netID)) == 
          as.numeric(i)
        pred.no <- nrow(datap[ind.pred, ])
        if (site.no * pred.no != 0) {
          workspace.name.a <- paste("dist.net", i, ".a.RData", 
                                    sep = "")
          workspace.name.b <- paste("dist.net", i, ".b.RData", 
                                    sep = "")
          path.a <- file.path(object$ssn.object@path, 
                              "distance", predpointsID, workspace.name.a)
          path.b <- file.path(object$ssn.object@path, 
                              "distance", predpointsID, workspace.name.b)
          file_handle = file(path.a, open = "rb")
          distmata <- unserialize(file_handle)
          close(file_handle)
          file_handle = file(path.b, open = "rb")
          distmatb <- unserialize(file_handle)
          close(file_handle)
          ordoi <- order(as.numeric(rownames(distmata)))
          ordpi <- order(as.numeric(rownames(distmatb)))
          ni <- length(ordoi)
          nj <- length(ordpi)
          distmata <- distmata[ordoi, ordpi, drop = F]
          distmatb <- distmatb[ordpi, ordoi, drop = F]
          if (any(rownames(distmata) != colnames(distmatb))) 
            stop("rownames of distmata do not match colnames of distmatb")
          if (any(colnames(distmata) != rownames(distmatb))) 
            stop("colnames of distmata do not match rownames of distmatb")
          dist.junc.a[(nSoFari + 1):(nSoFari + ni), (nSoFarj + 
                                                       1):(nSoFarj + nj)] <- distmata
          dist.junc.b[(nSoFarj + 1):(nSoFarj + nj), (nSoFari + 
                                                       1):(nSoFari + ni)] <- distmatb
          net.zero[(nSoFari + 1):(nSoFari + ni), (nSoFarj + 
                                                    1):(nSoFarj + nj)] <- 1
        }
        else {
          ni <- site.no
          nj <- pred.no
        }
        nSoFari <- nSoFari + ni
        nSoFarj <- nSoFarj + nj
      }
      rnames <- datao[distord, "pid"]
      cnames <- datap[, "pid"]
      rownames(dist.junc.a) <- rnames
      colnames(dist.junc.a) <- cnames
      a.mat <- pmax(dist.junc.a, t(dist.junc.b))
      b.mat <- pmin(dist.junc.a, t(dist.junc.b))
      dist.hydro <- as.matrix(dist.junc.a + t(dist.junc.b))
      if (length(grep("tailup", CorModels)) > 0) {
        if (length(grep("tailup", CorModels)) > 1) 
          stop("Cannot have more than 1 tailup model")
        flow.con.mat <- 1 - (b.mat > 0) * 1
        addfunccol <- object$args$addfunccol
        w.matrix <- sqrt(pmin(outer(datao[distord, addfunccol], 
                                    rep(1, times = npred)), t(outer(datap[, addfunccol], 
                                                                    rep(1, times = nobs))))/pmax(outer(datao[distord, 
                                                                                                             addfunccol], rep(1, times = npred)), t(outer(datap[, 
                                                                                                                                                                addfunccol], rep(1, times = nobs))))) * flow.con.mat * 
          net.zero
        if (any(rownames(w.matrix) != rownames(a.mat))) 
          stop("rownames of w.matrix do not match rownames of a.mat")
        if (any(colnames(w.matrix) != colnames(a.mat))) 
          stop("colnames of w.matrix do not match colnames of a.mat")
      }
    }
    xyobs <- ocoord
    x.samp <- ocoord[distord, 1, drop = F]
    y.samp <- ocoord[distord, 2, drop = F]
    xypred <- pcoord
    x.pred <- pcoord[, 1, drop = F]
    y.pred <- pcoord[, 2, drop = F]
    if (any(rownames(x.samp) != rownames(a.mat))) 
      stop("rownames of x.samp do not match rownames of a.mat")
    if (any(rownames(x.pred) != colnames(a.mat))) 
      stop("rownames of x.pred do not match colnames of a.mat")
    REPs <- NULL
    if (!is.null(REs)) {
      REnames <- names(REs)
      for (ii in 1:length(REnames)) if (any(is.na(datap[, 
                                                        REnames[ii]]))) 
        stop("Cannot having missing values when creating random effects")
      REOs <- list()
      REPs <- list()
      ObsSimDF <- datao
      PredSimDF <- datap
      for (ii in 1:length(REnames)) {
        plevels <- unique(c(levels(PredSimDF[, REnames[[ii]]]), 
                            paste("o", levels(ObsSimDF[, REnames[[ii]]]), 
                                  sep = ""), paste("p", levels(PredSimDF[, 
                                                                         REnames[[ii]]]), sep = "")))
        pino <- PredSimDF[, REnames[[ii]]] %in% ObsSimDF[, 
                                                         REnames[[ii]]]
        ObsSimDF[, REnames[[ii]]] <- paste("o", ObsSimDF[, 
                                                         REnames[[ii]]], sep = "")
        ObsSimDF[, REnames[[ii]]] <- as.factor(as.character(ObsSimDF[, 
                                                                     REnames[[ii]]]))
        levels(PredSimDF[, REnames[[ii]]]) <- plevels
        if (any(pino)) 
          PredSimDF[pino, REnames[[ii]]] <- paste("o", 
                                                  PredSimDF[pino, REnames[[ii]]], sep = "")
        if (any(!pino)) 
          PredSimDF[!pino, REnames[[ii]]] <- paste("p", 
                                                   PredSimDF[!pino, REnames[[ii]]], sep = "")
        PredSimDF[, REnames[[ii]]] <- as.factor(as.character(PredSimDF[, 
                                                                       REnames[[ii]]]))
        blevels <- unique(c(levels(ObsSimDF[, REnames[[ii]]]), 
                            levels(PredSimDF[, REnames[[ii]]])))
        ObsSimDF[, REnames[[ii]]] <- factor(ObsSimDF[, 
                                                     REnames[[ii]]], levels = blevels, ordered = FALSE)
        PredSimDF[, REnames[[ii]]] <- factor(PredSimDF[, 
                                                       REnames[[ii]]], levels = blevels, ordered = FALSE)
        REOs[[ii]] <- model.matrix(~ObsSimDF[distord, 
                                             REnames[[ii]]] - 1)
        REPs[[ii]] <- model.matrix(~PredSimDF[, REnames[[ii]]] - 
                                     1)
        rownames(REOs[[ii]]) <- datao[distord, "pid"]
        rownames(REPs[[ii]]) <- datap[, "pid"]
        if (any(rownames(REOs[[ii]]) != rownames(a.mat))) 
          stop("rownames RE for obs do not match rownames of a.mat")
        if (any(rownames(REPs[[ii]]) != colnames(a.mat))) 
          stop("rownames RE for preds do not match colnames of a.mat")
      }
      for (ii in 1:length(REnames)) REPs[[ii]] <- REOs[[ii]] %*% 
        t(REPs[[ii]])
    }
    Vpred <- makeCovMat(theta = theta, dist.hydro = dist.hydro, 
                        a.mat = a.mat, b.mat = b.mat, w.matrix = w.matrix, 
                        net.zero = net.zero, x.row = x.samp, y.row = y.samp, 
                        x.col = x.pred, y.col = y.pred, CorModels = CorModels, 
                        useTailDownWeight = useTailDownWeight, use.nugget = FALSE, 
                        use.anisotropy = object$args$use.anisotropy, REs = REPs)
    Vpred <- Vpred[ind, , drop = F]
  }
  response.col <- object$args$zcol
  mod.names <- as.character(attr(terms(object$args$formula, 
                                       data = object$ssn.object@obspoints@SSNPoints[[1]]@point.data), 
                                 "variables"))
  nc.tmp <- length(mod.names)
  ind.allcov <- rep(TRUE, times = length(datap[, 1]))
  if (nc.tmp > 2) {
    for (i in 3:nc.tmp) ind.allcov <- ind.allcov & !is.na(datap[, 
                                                                mod.names[i]])
  }
  np.allcov <- sum(ind.allcov)
  datap[, response.col] <- NA
  datap[, paste(response.col, ".predSE", sep = "")] <- NA
  datap1 <- datap[ind.allcov, ]
  datap1[, response.col] <- -1
  datap1[, paste(response.col, ".predSE", sep = "")] <- NA
  formula <- object$args$formula
  mf <- model.frame(formula, data = datap1)
  mt <- attr(mf, "terms")
  Xpred <- model.matrix(mt, mf, contrasts)
  Xpred <- Xpred[, object$sampinfo$cutX1toX2, drop = F]
  sumparsil <- sum(theta[attr(theta, "type") == "parsill"])
  Vi <- object$estimates$Vi
  covb <- object$estimates$covb
  Xobs <- object$sampinfo$X
  z <- object$sampinfo$z
  n <- object$sampinfo$obs.sample.size
  p <- object$sampinfo$rankX
  parsilvec <- rep(sumparsil, times = length(Vpred[1, ]))
  if (object$args$family == "poisson") {
    beta.hat <- object$estimates$betahat
    eta.hatp <- Xpred[ind.allcov, ] %*% beta.hat
    eta.hato <- Xobs %*% beta.hat
    Del.ip <- as.vector(1/exp(eta.hatp))
    Del.io <- as.vector(1/exp(eta.hato))
    A.5p <- as.vector(sqrt(exp(eta.hatp)))
    A.5o <- as.vector(sqrt(exp(eta.hato)))
    Vpred <- t((Del.ip * A.5p) * t(Del.io * A.5o * Vpred))
    parsilvec <- sumparsil * (A.5p * Del.ip)^2
  }
  if (object$args$family == "binomial") {
    beta.hat <- object$estimates$betahat
    eta.hatp <- Xpred[ind.allcov, ] %*% beta.hat
    eta.hato <- Xobs %*% beta.hat
    Del.ip <- as.vector((1 + exp(eta.hatp))^2/exp(eta.hatp))
    Del.io <- as.vector((1 + exp(eta.hato))^2/exp(eta.hato))
    A.5p <- as.vector(sqrt(exp(eta.hatp)/(1 + exp(eta.hatp))^2))
    A.5o <- as.vector(sqrt(exp(eta.hato)/(1 + exp(eta.hato))^2/object$sampinfo$trialsvec))
    Vpred <- t((Del.ip * A.5p) * t(Del.io * A.5o * Vpred))
    parsilvec <- sumparsil * (A.5p * Del.ip)^2
  }
  M <- rbind(Vpred, t(Xpred), parsilvec)
  XXSiXi <- Xobs %*% covb
  XSi <- t(Xobs) %*% Vi
  
  vars <- M[(n+1):(n+p),pid]
  vars["PDISRSA_1YR"] <- Xs[1]
  vars["POWNRCA_PRI"] <- Xs[2]
  vars["DAPOPRCA2010"] <- Xs[3]
  
  r1 <- (vars - XSi %*% M[1:n,pid])
  m <- covb %*% r1
  tlam <- t(M[1:n,pid] + XXSiXi %*% r1) %*% Vi
  p1 <- tlam %*% z
  
  p2 <- target
  (p2 - p1)^2

#   pred.out <- t(apply(M, 2, UK4Apply, covb = covb, XXSiXi = XXSiXi, 
#                       XSi = XSi, Vi = Vi, z = z, n = n, p = p))
#   datap1[, response.col] <- pred.out[, 1]
#   datap1[, paste(response.col, ".predSE", sep = "")] <- pred.out[,2]

#  datap[ind.allcov, ] <- datap1
#   for (i in 1:length((object$ssn.object@predpoints@SSNPoints))) if (object$ssn.object@predpoints@ID[i] == 
#                                                                       predpointsID) {
#     object$ssn.object@predpoints@SSNPoints[[i]]@point.data <- datap[order(distordp), 
#                                                                     ]
#   }
#   object$args$predpointsID <- predpointsID
#   class(object) <- "glmssn.predict"
#   object
}

#### ####
  makeCovMat <-
  function(theta, dist.hydro, a.mat, b.mat, w.matrix = NULL,
           net.zero, x.row, y.row, x.col, y.col, useTailDownWeight,
           CorModels, use.nugget, use.anisotropy, REs)
  {
    
    nRow <- length(x.row)
    nCol <- length(x.col)
    
    if(is.null(net.zero)) net.zero <- matrix(1, nrow = nRow, ncol = nCol)
    V <- matrix(0, nrow = nRow, ncol = nCol)
    
    # create covariance matrix component for tailup models
    npar.sofar <- 0
    if(length(grep("tailup",CorModels)) > 0){
      if(length(grep("tailup",CorModels)) > 1)
        stop("Cannot have more than 1 tailup model")
      funname <- tolower(paste(substr(unlist(strsplit(CorModels,".", 
                                                      fixed = T))[(1:length(unlist(strsplit(CorModels,".", 
                                                                                            fixed = T))))[unlist(strsplit(CorModels,".", fixed = T)) == 
                                                                                                            "tailup"] - 1], 1, 3),".tailup", sep = ""))
      tailupmod <- call(funname, dist.hydro=dist.hydro, 
                        weight = w.matrix, parsil = theta[npar.sofar + 1], 
                        range1 = theta[npar.sofar + 2])
      V <- V + eval(tailupmod)*net.zero
      npar.sofar <- npar.sofar + 2
    }
    # create covariance matrix component for taildown models
    if(length(grep("taildown",CorModels)) > 0){
      if(length(grep("taildown",CorModels)) > 1)
        stop("Cannot have more than 1 taildown model")
      funname <- tolower(paste(substr(unlist(strsplit(CorModels,".", 
                                                      fixed = T))[(1:length(unlist(strsplit(CorModels,".", 
                                                                                            fixed = T))))[unlist(strsplit(CorModels,".", fixed = T)) == 
                                                                                                            "taildown"] - 1], 1, 3),".taildown", sep = ""))
      taildnmod <- call(funname, dist.hydro=dist.hydro, 
                        a.mat = a.mat, b.mat = b.mat, parsil = theta[npar.sofar + 1], 
                        useTailDownWeight = useTailDownWeight, weight = w.matrix,
                        range1 = theta[npar.sofar + 2])
      V <- V + eval(taildnmod)*net.zero
      npar.sofar <- npar.sofar + 2
    }
    # create covariance matrix componenet for Euclidean models
    if(length(grep("Euclid",CorModels)) > 0){
      if(length(grep("Euclid",CorModels)) > 1)
        stop("Cannot have more than 1 Euclidean model")
      npar.parsil <- npar.sofar + 1
      if(use.anisotropy == FALSE) {
        dist.mat <- distGeo(x.row, y.row, x.col, y.col, 
                            theta[npar.sofar + 2])
        npar.sofar <- npar.sofar + 2
      }
      else {
        dist.mat <- distGeo(x.row, y.row, x.col, y.col, 
                            theta[npar.sofar + 2], theta[npar.sofar + 3], 
                            theta[npar.sofar + 4])
        npar.sofar <- npar.sofar + 4
      }
      funname <- paste(tolower(substr(unlist(strsplit(CorModels,".", 
                                                      fixed = T))[(1:length(unlist(strsplit(CorModels,".", 
                                                                                            fixed = T))))[unlist(strsplit(CorModels,".", fixed = T)) == 
                                                                                                            "Euclid"] - 1], 1, 3)),".Euclid", sep = "")
      taileumod <- call(funname, distance.matrix = dist.mat,
                        parsil = theta[npar.parsil])
      V <- V + eval(taileumod)
    }
    
    if(length(REs)) {
      for(ii in 1:length(REs)) {
        npar.sofar <- npar.sofar + 1
        V <- V + theta[npar.sofar]*REs[[ii]]
      }
    }
    
    # create diagonal covariance matrix component for nugget effect
    if(use.nugget == TRUE) {
      if(nRow != nCol) stop(		
        "covariancd matrix asymmetric -- cannot use nugget")
      npar.sofar <- npar.sofar + 1
      V <- V + diag(theta[npar.sofar], nrow = nRow, ncol = nCol)
    } else if(nRow == nCol){
      V + diag(1e-6, nrow = nRow, ncol = nCol)
    }
    
    V
    
  }

exp.tailup <- function(dist.hydro, weight, parsil = parsil, range1 = range1)
{
  parsil*exp(-3*dist.hydro/range1)*weight
}

exp.taildown <- function(dist.hydro, a.mat, b.mat, parsil = parsil, 
                         range1 = range1, useTailDownWeight, weight = NULL)
{
  V <- parsil*exp(-3*dist.hydro/range1)
  if(useTailDownWeight == TRUE) V <- V*weight
  V
}

exp.Euclid <-
  function(parsil, distance.matrix)
  {
    parsil*exp(-3*distance.matrix)
  }

distGeo <-
  function(x.row, y.row, x.col, y.col, range = 1, minorp = 1, rotate = 0)
  {
    if(range < 0)
      stop("range parameter less than 0 in distGeo")
    if(rotate < 0 || rotate > 180)
      stop("rotation parameter beyond 0 - 180 range in distGeo")
    if(minorp < 0 || minorp > 1)
      stop("minor range proportion beyond 0 - 1 range")
    # total number of row observations
    n.row <- length(x.row)
    # total number of column observations
    n.col <- length(x.col)
    # expand all x-coordinates for rows
    sxr <- outer(as.vector(x.row), rep(1, times = n.col))
    # expand all x-coordinates for columns
    sxc <- outer(rep(1, times = n.row), as.vector(x.col))
    # find difference in x-coordinates between all pairwise locations
    sxdif <- sxr - sxc
    # expand all x-coordinates for rows
    syr <- outer(as.vector(y.row), rep(1, times = n.col))
    # expand all x-coordinates for columns
    syc <- outer(rep(1, times = n.row), as.vector(y.col))
    # find difference in x-coordinates between all pairwise locations
    sydif <- syr - syc
    # rotate coordinates
    newx <- cos(rotate*.0174533)*sxdif - sin(rotate*.0174533)*sydif
    newy <- sin(rotate*.0174533)*sxdif + cos(rotate*.0174533)*sydif
    # scale coordinates by minor and major axes */
    newx <- newx/(range*minorp)
    newy <- newy/range
    # compute distance for the scaled and rotated coordinates */
    sqrt(newx^2 + newy^2)
  }

UK4Apply <-
  function(vec, covb, XXSiXi, XSi, Vi, z, n, p)
  {
    r1 <- (vec[(n+1):(n+p)] - XSi %*% vec[1:n])
    m <- covb %*% r1
    tlam <- t(vec[1:n] + XXSiXi %*% r1) %*% Vi
    cbind(tlam %*% z,
          sqrt(vec[n+p+1] - tlam %*% vec[1:n] + t(m) %*% vec[(n+1):(n+p)]))
  }


SSNpred.site <- function(vars = NULL, ssn, pid) {
  library(SSN)
  ssn1.glmssn.EEE.5 <- ssn
  
    #get the data frame and build the scenario
  preds.5 <- getSSNdata.frame(ssn1.glmssn.EEE.5, Name = "preds")
  if (length(vars > 0)) {
    preds.5[preds.5$pid == pid,'PDISRSA_1YR'] <- vars["PDISRSA_1YR"]
    preds.5[preds.5$pid == pid,'DAPOPRCA2010'] <- vars["DAPOPRCA2010"]
    preds.5[preds.5$pid == pid,'POWNRCA_PRI'] <- vars["POWNRCA_PRI"]
  }
  preds.5 <- preds.5[match(pid.order,preds.5$pid),]
  row.names(preds.5) <- preds.5$pid
  ssn1.glmssn.EEE.5 <- putSSNdata.frame(preds.5, ssn1.glmssn.EEE.5, Name = "preds")
  
  #Run the prediction
  ssn1.glmssn.EEE.5.preds <- predict(ssn1.glmssn.EEE.5,predpointsID = "preds", newdata = 'preds')
  
  #Check the results
  preds.5.dis <- getSSNdata.frame(ssn1.glmssn.EEE.5.preds, Name = 'preds')
  
  critval <- qnorm(0.975)
  preds.5.dis$uci <- preds.5.dis$log10_FSS_26Aug14 + (critval * preds.5.dis$log10_FSS_26Aug14.predSE)
  preds.5.dis$lci <- preds.5.dis$log10_FSS_26Aug14 - (critval * preds.5.dis$log10_FSS_26Aug14.predSE)
      
  #outputs FSS prediction
  return(c(predict = preds.5.dis[preds.5.dis$pid == pid,'log10_FSS_26Aug14'],
           predict.se = preds.5.dis[preds.5.dis$pid == pid,"log10_FSS_26Aug14.predSE"],
           pred.int.up = preds.5.dis[preds.5.dis$pid == pid,"uci"],
           pred.int.low = preds.5.dis[preds.5.dis$pid == pid,"lci"],
           observed = preds.5[preds.5$pid == imp.pid[i,'pid'],'log10_FSS_26Aug14']))
}
#### ####
#load('predOptim.Rdata')
load('ssn1_glmssn_EEE.Rdata')
load('minmax.Rdata')
load('precip_daily_sum_1095_days.Rdata')

impaired <- data.frame(STATION_KEY = c(21842,34660,21792,33361,26818,33418,33417,34695,26822,33320,33333,30403,34665,26816,25297,26964,29906,33327),
                       TMDL_Target = c(14,14,14,3,7,rep(14,6),8,14,8,14,8,14,14))
impaired$TMDL_Target_Scaled_log <- (log10(impaired$TMDL_Target)-min.max[min.max$variable == 'log10_FSS_26Aug14','min_val'])/(min.max[min.max$variable == 'log10_FSS_26Aug14','max_val']-min.max[min.max$variable == 'log10_FSS_26Aug14','min_val'])


obs <- getSSNdata.frame(ssn1.glmssn.EEE, Name = 'Obs')

#Generate table of quantiles
library(plyr)
qall <- ddply(dfdall,.(STATION_KEY),function(x){quantile(x$sum_1095_days, probs = seq(0,1,0.1))})
qall.scaled <- data.frame(STATION_KEY = qall$STATION_KEY, as.data.frame(lapply(qall[,setdiff(names(qall),'STATION_KEY')],
                     function(x) {(x-min.max[min.max$variable == 'sum_1095_days','min_val'])/
                                    (min.max[min.max$variable == 'sum_1095_days','max_val']-
                                       min.max[min.max$variable == 'sum_1095_days','min_val'])})))

library(dplyr)
max.obs <- data.frame(obs %>% group_by(STATION_KEY) %>% filter(log10_FSS_26Aug14 == max(log10_FSS_26Aug14)))
max.obs <- max.obs[!duplicated(max.obs$STATION_KEY),]
preds <- getSSNdata.frame(ssn1.glmssn.EEE, Name = "preds")
pid.order <- preds$pid
preds <- rename(preds, STATION_KEY = STATION_KE)
preds <- merge(preds, max.obs[,c('STATION_KEY','sum_1095_days','PALITHERODRCA','PDISRSA_1YR',
                                 'PASILTRCA','DAPOPRCA2010','POWNRCA_PRI','log10_FSS_26Aug14')],by = 'STATION_KEY',all.x = TRUE)
preds$STATION_KEY <- as.character(preds$STATION_KEY)
preds <- preds[match(pid.order,preds$pid),]
row.names(preds) <- preds$pid
ssn1.glmssn.EEE <- putSSNdata.frame(preds, ssn1.glmssn.EEE, Name = "preds")

#Fill in sum_1095_days with the 90th percentile of 3 year sum from 1995 through 2012 at each impaired site
preds.p <- preds
preds.p[preds.p$STATION_KEY %in% qall.scaled$STATION_KEY,'sum_1095_days'] <- qall.scaled[order(match(qall.scaled$STATION_KEY,
                                                                                                 preds.p[preds.p$STATION_KEY %in% qall.scaled$STATION_KEY,
                                                                                                       'STATION_KEY'])),'X90.']
preds.p <- preds.p[match(pid.order,preds.p$pid),]
row.names(preds.p) <- preds.p$pid
ssn1.glmssn.EEE.p <- putSSNdata.frame(preds.p, ssn1.glmssn.EEE, Name = "preds")

imp.pid <- merge(preds[,c('STATION_KEY','pid')], impaired, all.y = T, by ="STATION_KEY")

result <- list()

for (i in 1:nrow(imp.pid)) {
  tmp <- optim(preds[preds$pid == as.character(imp.pid[i,"pid"]),c("PDISRSA_1YR","POWNRCA_PRI","DAPOPRCA2010")], 
               predict.PTB,
               target = imp.pid[i,'TMDL_Target_Scaled_log'], 
               pid = as.character(imp.pid[i,"pid"]), 
               object = ssn1.glmssn.EEE.p, 
               predpointsID='preds', 
               lower = 0,upper = 1,method = "L-BFGS-B") #, lower = 0,upper = 1, method = "L-BFGS-B"
  tmp$imp.info <- imp.pid[i,]
  tmp$p.var <- preds.p[preds.p$pid == as.character(imp.pid[i,"pid"]),"sum_1095_days"]
  tmp$p1.SSN <- SSNpred.site(vars = tmp$par, ssn = ssn1.glmssn.EEE.p, pid = imp.pid[i,"pid"])[1:4]
  tmp$orig.vars <- preds[preds$pid == as.character(imp.pid[i,"pid"]),c("PDISRSA_1YR","POWNRCA_PRI","DAPOPRCA2010","sum_1095_days",'PALITHERODRCA','PASILTRCA')]
  tmp$p1.orig <- SSNpred.site(ssn = ssn1.glmssn.EEE, pid = imp.pid[i,"pid"])
  result <- append(result, list(tmp))
}

fss.group <- read.csv('//deqhq1/tmdl/tmdl_wr/midcoast/models/sediment/target_development/final/R_output_bugs_CART_ALL_final_2013-06-15_trans.csv')

tmp.vars <- obs.complete[,c('STATION_KEY','SVN','Ref_Samples','PDISRSA_1YR','POWNRCA_PRI','DAPOPRCA2010','sum_1095_days','PALITHERODRCA','PASILTRCA')]

vars.group <- merge(tmp.vars,fss.group[,c('STATION_KEY','Group')],by='STATION_KEY',all.x = TRUE)

imp <- vars.group[vars.group$SVN %in% rid.cross[rid.cross[,'STATION_KEY'] %in% imp.pid[,'STATION_KEY'],'SVN'],c('SVN','STATION_KEY','Group',"PDISRSA_1YR","POWNRCA_PRI","DAPOPRCA2010",'sum_1095_days','PALITHERODRCA','PASILTRCA')]
imp$group <- 'Impaired'

ref <- vars.group[vars.group$Ref_Samples == '1',c('SVN','STATION_KEY','Group',"PDISRSA_1YR","POWNRCA_PRI","DAPOPRCA2010",'sum_1095_days','PALITHERODRCA','PASILTRCA')]
ref$group <- 'Reference'

toplot <- rbind(imp, ref)
toplot$toplot.name <- paste(toplot$group, toplot$Group)
group = 'Group 4'
toplot <- toplot[toplot$Group == group,]
boxplot(toplot$PDISRSA_1YR~toplot$toplot.name,main='PDISRSA_1YR')
boxplot(toplot$POWNRCA_PRI~toplot$toplot.name,main='POWNRCA_PRI')
boxplot(toplot$DAPOPRCA2010~toplot$toplot.name,main='DAPOPRCA2010')
boxplot(toplot$sum_1095_days~toplot$toplot.name,main='sum_1095_days')
boxplot(toplot$PALITHERODRCA~toplot$toplot.name,main='PALITHERODRCA')
boxplot(toplot$PASILTRCA~toplot$toplot.name,main='PASILTRCA')
