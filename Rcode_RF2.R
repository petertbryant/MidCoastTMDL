
################################################################################################
#http://evansmurphy.wix.com/evansspatial#!random‐forests‐model‐select/cksm

rf.modelSel <- function(xdata, ydata, imp.scale="mir", r=c(0.25, 0.50, 0.75),  
                        final=FALSE, plot.imp=TRUE, parsimony=NULL, ...) 
{
  if (!require (randomForest)) stop("randomForest PACKAGE MISSING")
  rf.ImpScale <- function (x, scale="mir") { 
    if (!inherits(x, "randomForest")) 
      stop(deparse(substitute(x)), " Must be a randomForest object")
    if (x$type == "regression") {
      if (is.null(x$importanceSD) == TRUE | "%IncMSE" %in% 
            names(as.data.frame(x$importance)) == FALSE)
        stop("OBJECT DOES NOT CONTAIN PERMUTATED IMPORTANCE, PLEASE RUN 
             randomForest WITH importance=TRUE")   
      rf.imp <- x$importance[,"%IncMSE"]
      rf.impSD <- x$importanceSD
      rf.impSD[rf.impSD == 0] <- 0.000000001  
      if (scale == "mir") {
        i <- rf.imp / max(rf.imp) 
      }	  
      if (scale == "se") {
        i <- ( rf.imp / rf.impSD ) / sum(rf.imp / rf.impSD, na.rm=TRUE)			 
      }
    }
    if (x$type == "classification" | x$type == "unsupervised") {
      if (is.null(x$importanceSD) == TRUE | "MeanDecreaseAccuracy" %in% 
            names(as.data.frame(x$importance)) == FALSE)
        stop("OBJECT DOES NOT CONTAIN PERMUTATED IMPORTANCE, PLEASE RUN 
             randomForest WITH importance=TRUE") 
      rf.imp <- x$importance[,"MeanDecreaseAccuracy"]
      rf.impSD <- x$importanceSD[,"MeanDecreaseAccuracy"]
      rf.impSD[rf.impSD == 0] <- 0.000000001	
      if (scale == "mir") {
        i <- rf.imp / max(rf.imp) 
      }	  
      if (scale == "se") {
        i <- ( rf.imp / rf.impSD ) / sum(rf.imp / rf.impSD, na.rm=TRUE)			 
      }
    }
    i <- as.data.frame(i)
    names(i) <- "importance" 
    row.names(i) <- names(rf.imp)	
    return( i )            
  }
  
  RFtype <- is.factor(ydata) #TEST FOR FACTOR IN Y 
  ##CLASSIFICATION##
  if (RFtype == "TRUE") {
    model.vars <- list()
    ln <- 0
    rf.all <- randomForest(x=xdata, y=ydata, importance=TRUE, ...) 
    model.vars[[ln <- ln + 1]] <- rownames(rf.all$importance)  
    class.errors <- as.data.frame(rf.all$err.rate)
    class.errors <- na.omit(class.errors)  
    class.errors[class.errors == NaN] <- 0
    class.errors[class.errors == Inf] <- 1         
    i <- vector()
    for ( l in 2:nlevels(as.factor(names(class.errors))) ) {              
      i <- append(i, median(class.errors[,l]))
    }        
    max.error = max(i) 
    imp <- rf.ImpScale(rf.all, scale=imp.scale) 
    results <- as.data.frame(array(0, dim=c( 0, 4 )))
    x <- c(0, (median(rf.all$err.rate[,"OOB"]) * 100), max.error * 100, dim(xdata)[2] )
    results <- rbind(results, x) 	 	 
    for (p in 1:length(r) ) {
      t = quantile(imp[,1], probs=r[p])
      sel.imp <- subset(imp, importance >= t)
      sel.vars <- rownames(sel.imp)
      if (length( sel.vars ) > 1) {                             
        xdata.sub <- xdata[,sel.vars]       
        rf.model <- randomForest(x=xdata.sub, y=ydata, importance=TRUE)          
        class.errors <- as.data.frame(rf.model$err.rate)
        class.errors <- na.omit(class.errors)  
        class.errors[class.errors == NaN] <- 0
        class.errors[class.errors == Inf] <- 1      
        i <- as.vector(array(0, dim=c((0),(1))))
        for ( l in 2:nlevels(as.factor(names(class.errors))) )
        {
          x.bar <- mean(class.errors[,l])              
          i <- as.vector(append(i, x.bar, after=length(i) ))
        }        
        max.error = max(i[2:length(i)] )     
        x <- c(t, median(rf.model$err.rate[,1]) * 100, max.error * 100, length(sel.vars) )
        results <- rbind(results, x)
        model.vars[[ln <- ln + 1]] <- rownames(rf.model$importance)
      }
    }
    names(results) <- c("THRESHOLD", "OOBERROR", "CLASS.ERROR", "NPARAMETERS")
    results <- results[order(results$CLASS.ERROR, results$OOBERROR, results$NPARAMETERS),]
    if (is.null(parsimony) == FALSE) { 
      if(parsimony < 0.00000001 | parsimony > 0.9) stop( "parsomony MUST RANGE 0-1")
      oob <- "TRUE"
      for(i in 2:nrow(results)) {
        if( abs((results[i,][2] - results[1,][2] ) / results[1,][2]) <= parsimony  &
              abs( (results[i,][3] - results[1,][3] ) / results[1,][3] ) <= parsimony ) {
          oob <- append(oob, "TRUE")
        } else {
          oob <- append(oob, "FALSE")
        }
        final <- results[which( oob == "TRUE" ),]
        final <- final[final$NPARAMETERS == min(final$NPARAMETERS) ,]$THRESHOLD
      }
    } else {		
      final <- as.vector(results[,"THRESHOLD"])[1]
    }	
    sel.imp <- subset(imp, importance >= final)    
    sel.vars <- rownames(sel.imp)
    sel.post=which( results$NPARAMETERS == length(sel.vars) ) 
    results <- rbind(results[sel.post,],results[-sel.post,]) 	
  } # END OF CLASSIFICATION
  
  ##REGRESSION## 
  if (RFtype == "FALSE") {
    model.vars <- list()
    ln <- 0      
    rf.all <- randomForest(x=xdata, y=ydata, importance=TRUE, ...) 
    model.vars[[ln <- ln + 1]] <- rownames(rf.all$importance)  
    imp <- rf.ImpScale(rf.all, scale=imp.scale) 
    results <- as.data.frame(array(0, dim=c( 0, 4 )))
    x <- c(0, (median(rf.all$rsq)), mean(rf.all$mse), dim(xdata)[2] )
    results <- rbind(results, x)     
    for (p in 1:length(r) ) {
      t = quantile(imp[,1], probs=r[p])		 
      sel.vars <- rownames(subset(imp, importance >= t))  
      if (length( sel.vars ) > 1) {                             
        xdata.sub <- as.data.frame(xdata[,sel.vars]) 
        rf.model <- randomForest(x=xdata.sub, y=ydata, importance=TRUE, ...)          
        x <- c(t, (median(rf.model$rsq)), mean(rf.model$mse), length(sel.vars) )
        results <- rbind(results, x)
        model.vars[[ln <- ln + 1]] <- rownames(rf.model$importance)
      }
    }
    names(results) <- c("THRESHOLD", "VAREXP", "MSE", "NPARAMETERS")
    results <- results[order(-results$VAREXP, results$MSE, results$NPARAMETERS),]  
    if (!is.null(parsimony)) {
      if(parsimony < 0.00000001 | parsimony > 0.9) stop( "parsomony MUST RANGE 0-1")	
      oob <- "TRUE"
      for(i in 2:nrow(results)) {
        if( abs((results[i,][2] - results[1,][2] ) / results[1,][2]) <= parsimony  &
              abs( (results[i,][3] - results[1,][3] ) / results[1,][3] ) <= parsimony ) {
          oob <- append(oob, "TRUE")
        } else {
          oob <- append(oob, "FALSE")
        }
        final <- results[which( oob == "TRUE" ),]
        final <- final[final$NPARAMETERS == min(final$NPARAMETERS) ,]$THRESHOLD
      }
    } else {		
      final <- as.vector(results[,"THRESHOLD"])[1]
    }	
    sel.imp <- subset(imp, importance >= final)    
    sel.vars <- rownames(sel.imp)
    sel.post=which( results$NPARAMETERS == length(sel.vars) ) 
    results <- rbind(results[sel.post,],results[-sel.post,]) 			
  } # END OF REGRESSION 
  
  if (plot.imp == TRUE) {
    if (imp.scale=="mir") {lable="Row Standardization Variable Importance"} 	
    if (imp.scale=="se") {lable="Standardized Error Variable Importance"}
    p <- as.matrix(subset(imp, importance >= final))    
    ord <- rev(order(p[,1], decreasing=TRUE)[1:dim(p)[1]])  
    dotchart(p[ord,1], main=lable, pch=19)
  }
  
  if (final == TRUE) {
    sub.xdata <- xdata[,sel.vars]  #SUBSET VARIABLES
    rf.final <- randomForest(x=sub.xdata, y=ydata, importance=TRUE, ...)           
    ( list(MODEL=rf.final, SELVARS=sel.vars, TEST=results, IMPORTANCE=sel.imp, PARAMETERS=model.vars) )      
  } else {
    ( list(SELVARS=sel.vars, TEST=results, IMPORTANCE=sel.imp, PARAMETERS=model.vars) ) 
  }     
}


################################################################################################

MultiColinear <- function(x, p=1e-07) {
  if (!inherits(x, "data.frame")) stop("X MUST BE A data.frame")
  if ( (dim(x)[2] < 2) == TRUE) stop("NOT ENOUGH VARIABLES TO TEST")
  xtest <- x
  x.names <- names(xtest)
  qrx <- qr(xtest, tol=p)
  if (length(names(xtest)[qrx$pivot[1:qrx$rank]]) != length(xtest) )  
  {  
    keep <- names(xtest)[qrx$pivot[1:qrx$rank]]
    warning("MULTI-COLINEAR VARIABLES: ", paste(setdiff(x.names, keep), collapse = ","))
    return(paste(setdiff(x.names, keep)))
  } else { print(" NO MULTICOLINEAR VARIABLES IDENTIFIED")
  } 
}
################################################################################################

rf.Permutation <- function(x, y, q=0.99, p=0.05, nperm=999, plot=TRUE, ...) {
  if (!require (randomForest)) stop("randomForest PACKAGE MISSING")
  if (is.factor(y)) stop("y CANNOT BE A FACTOR") 
  rf.test <- randomForest(x=x, y=y, ...)
  test.rsq <- median(rf.test$rsq)
  rand.dist <- vector() 
  for( i in 1:nperm) {	
    rand.y <- sample(y, length(y)) 
    rf.test <- randomForest(x=x, y=rand.y, ...)
    rand.dist <- append(rand.dist, median(rf.test$rsq)) 
  }	
  Pval <- function(x, test, nperm) { 
    if ( length( x[x >= test.rsq] ) < 1 ) { 
      error = 1 
    } else { 
      error = length( x[x >= test.rsq] ) + 1
    }	
    return( error / nperm )
  }				
  if( plot == TRUE) { 
    den=density(rand.dist)
    den$y <- den$y/max(den$y)		
    plot(den, type="n", xlim=c(min(rand.dist), 1), xlab="R-square", ylab="",  
         main="Distribution of randomized models")
    polygon(den, col="blue")
    abline(v=test.rsq, col="black", lwd=1.5, lty=2)
    abline(v=quantile(rand.dist,p=q),lwd=1.5, lty=2, col="red") 
    legend("topright", c("model", "null"), bg="white",  
           col=c("black","red"), lty=c(2,2), lwd=c(1.5,1.5) )				   
  }
  pValue=round(Pval(x=rand.dist, test=test.rsq, nperm=nperm), digits=6)	
  if( pValue <= p ) accept=TRUE else accept=FALSE 
  if (accept == TRUE) accept <- paste("MODEL SIGNIFICANT AT p=", pValue, sep="" ) 
  if (accept == FALSE) accept <- paste("MODEL NOT SIGNIFICANT AT p= ", pValue, sep="" )
  print(accept)		
  return( list(RandRsquare=rand.dist, Rsquare=test.rsq, Accept=accept, TestQuantile=q, 
               pValueThreshold=p, pValue=pValue, nPerm=nperm) )
} 
################################################################################################

library(randomForest)
library(reshape)

options(stringsAsFactors = FALSE)
vars <- read.csv("VarNames_RF.csv")
bugs <- read.csv("ssn_RF_data.csv")

#bugs <- arrange(bugs, STATION_KEY, desc(YEAR))
#bugs <- bugs[!duplicated(bugs$STATION_KEY),]

vars.fss2.s1 <- vars[vars$fss2.rf_keep == 1,]
fss2.s1 <- bugs[,colnames(bugs) %in% vars.fss2.s1$var]

colnames(fss2.s1)

fss2.s1$PSUSCEP_DE <- fss2.s1$PSUSCEP4_DE + fss2.s1$PSUSCEP5_DE
fss2.s1$PSUSCEP_LI <- fss2.s1$PSUSCEP4_LI + fss2.s1$PSUSCEP5_LI
fss2.s1$PASUSCEP_DE <- fss2.s1$PASUSCEP4_DE + fss2.s1$PASUSCEP5_DE
fss2.s1$PASUSCEP_LI <- fss2.s1$PASUSCEP4_LI + fss2.s1$PASUSCEP5_LI

fss2.s1$PSUSCEP4_DE <- NULL
fss2.s1$PSUSCEP5_DE <- NULL
fss2.s1$PSUSCEP4_LI <- NULL
fss2.s1$PSUSCEP5_LI <- NULL
fss2.s1$PASUSCEP4_DE <- NULL
fss2.s1$PASUSCEP5_DE <- NULL
fss2.s1$PASUSCEP4_LI <- NULL
fss2.s1$PASUSCEP5_LI <- NULL

# remove NAs in response variable
fss2.s1 <- fss2.s1[(!is.na(fss2.s1$FSS_26Aug14)),]

colnames(fss2.s1)

# impute the NAs
set.seed(111)
fss2.s1.imputed <- rfImpute(FSS_26Aug14 ~ ., fss2.s1, ntree=2000, iter=3)

colnames(fss2.s1)

# impute the NAs
set.seed(111)
fss1.s1.imputed <- rfImpute(FSS_26Aug14 ~ ., fss1.s1, ntree=2000, iter=3)

colnames(fss1.s1.imputed)

# initialize the variable importance df
fss2.s1.vi <- data.frame(matrix(, nrow = ncol(fss2.s1.imputed)-1, ncol = 50))
fss2.s1.visd <- data.frame(matrix(, nrow = ncol(fss2.s1.imputed)-1, ncol = 50))

fss2.s1.col <- colnames(fss2.s1.imputed)
fss2.s1.col <- fss2.s1.col[!(fss2.s1.col == "FSS_26Aug14")]

# WARNING - Takes about 1 hour
beg <- Sys.time()
set.seed(42)
for (i in 1:50) {
  fss2.s1.rf <- randomForest(FSS_26Aug14 ~ ., 
                             data = fss2.s1.imputed, 
                             ntree = 2000, 
                             keep.forest = TRUE, 
                             importance = TRUE)
  fss2.s1.vi[,i] <- fss2.s1.rf$importance[,1]
  fss2.s1.visd[,i] <- fss2.s1.rf$importanceSD
}
print(Sys.time() - beg)

# Add var names and index
fss1.s1.vi[,51]<- fss1.s1.col
fss1.s1.vi[,52]<-c(1:length(fss1.s1.col))
colnames(fss1.s1.vi)[51] <- "var_name"
colnames(fss1.s1.vi)[52] <- "var_index"
fss1.s1.visd[,51]<- fss1.s1.col
fss1.s1.visd[,52]<-c(1:length(fss1.s1.col))
colnames(fss1.s1.visd)[51] <- "var_name"
colnames(fss1.s1.visd)[52] <- "var_index"

# ----------
# Save the df with a timestamp so we don't accidently overwrite it.
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
save(fss2.s1.vi, file=paste0("fss2_s1_vi_",timestamp,".RData"))
save(fss2.s1.visd, file=paste0("fss2_s1_visd_",timestamp,".RData"))

#-----------
# Test for multicolinearity, remove them
cl <- MultiColinear(fss2.s1.imputed[,c(1:49,51:ncol(fss2.s1.imputed))], p=0.05)
xdata <- fss2.s1.imputed[,c(1:49,51:ncol(fss2.s1.imputed))] 
for(l in cl) {
  cl.test <- xdata[,-which(names(xdata)==l)]
  print(paste("REMOVE VARIABLE", l, sep=": "))
  MultiColinear(cl.test, p=0.05) 
}
for(l in cl) { fss2.s1.imputed <- fss2.s1.imputed[,-which(names(fss2.s1.imputed)==l)] }

#-----------
colnames(fss2.s1.imputed)

set.seed(42)
fss2.s1.rf <- rf.modelSel(xdata=fss2.s1.imputed[,2:ncol(fss2.s1.imputed)], 
                        ydata=fss2.s1.imputed[,"FSS_26Aug14"], 
                        imp.scale="se", r=c(0.1,0.25, 0.50, 0.75, 0.9),  
                        final=TRUE, plot.imp=TRUE, parsimony=NULL, ntree=2000) 


#-----------
fss2.s2.rf <- rf.Permutation(x=fss2.s1.imputed[,fss2.s1.rf$PARAMETERS[[5]]], 
               y=fss2.s1.imputed[,"FSS_26Aug14"],
               q=0.99, p=0.05, nperm=999, plot=TRUE)

# ----------
# Save the df with a timestamp so we don't accidently overwrite it.
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
save(fss2.s1.rf, file=paste0("fss2_s1_rf_varsel_",timestamp,".RData"))
save(fss2.s2.rf, file=paste0("fss2_s1_rf_nullsig_",timestamp,".RData"))

#-----------