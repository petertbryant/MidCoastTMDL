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