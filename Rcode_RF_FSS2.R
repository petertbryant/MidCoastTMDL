library(randomForest)
library(reshape)
library(plyr)

options(stringsAsFactors = FALSE)
vars <- read.csv("VarNames_RF.csv")
bugs <- read.csv("ssn_RF_data.csv")

#### Correlation Matrix Functions ####
# source
# http://stackoverflow.com/questions/15271103/how-to-modify-this-correlation-matrix-plot

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
#### Random forest step 1 #### 
# ----------------------------------------------------------- #
# FSS2 - Random forests excluding the physical habitat data
# ----------------------------------------------------------- #
# Step 1, per Genuer et al (2000), use brute force to identify important variables. 
# Since we have a large amount of variables random forest is run 50 times with a 
# very high number of trees per forest (ntree). This yields a distribution of importance scores.
# Removal is based on these distributions. We keep the variables with the highest scores.

vars.fss2.s1 <- vars[vars$fss2.rf_keep == 1,]
fss2.s1 <- bugs[,colnames(bugs) %in% vars.fss2.s1$var]

# keep a copy with all the data
fss2.s1.na <- fss2.s1

fss2.s1 <- fss2.s1[,!(colnames(fss2.s1) %in% c("PSUSCEP4_LI","PSUSCEP5_LI",
                                               "PASUSCEP4_LI", "PASUSCEP5_LI"))]
colnames(fss2.s1)

# remove NAs in response variable
fss2.s1 <- fss2.s1[(!is.na(fss2.s1$FSS_26Aug14)),]
# remove any NAs
fss2.s1 <- data.frame(na.omit(fss2.s1))

# mtry and ntree values 
mtry.fss2.s1 <- as.integer(((ncol(fss2.s1)-1) / 3),0)

# initialize the variable importance df
fss2.s1.vi <- data.frame(matrix(, nrow = ncol(fss2.s1)-1, ncol = 50))
fss2.s1.visd <- data.frame(matrix(, nrow = ncol(fss2.s1)-1, ncol = 50))

fss2.s1.col <- colnames(fss2.s1)
fss2.s1.col <- fss2.s1.col[!(fss2.s1.col == "FSS_26Aug14")]

# WARNING - Takes about 30 min
beg <- Sys.time()
set.seed(42)
for (i in 1:50) {
  fss2.s1.rf <- randomForest(FSS_26Aug14 ~ ., 
                             data = fss2.s1, 
                             ntree = 2000, 
                             keep.forest = TRUE, 
                             importance = TRUE)
  fss2.s1.vi[,i] <- fss2.s1.rf$importance[,1]
  fss2.s1.visd[,i] <- fss2.s1.rf$importanceSD
}
print(Sys.time() - beg)

# Add var names and index
fss2.s1.vi[,51]<- fss2.s1.col
fss2.s1.vi[,52]<-c(1:length(fss2.s1.col))
colnames(fss2.s1.vi)[51] <- "var_name"
colnames(fss2.s1.vi)[52] <- "var_index"
fss2.s1.visd[,51]<- fss2.s1.col
fss2.s1.visd[,52]<-c(1:length(fss2.s1.col))
colnames(fss2.s1.visd)[51] <- "var_name"
colnames(fss2.s1.visd)[52] <- "var_index"

# Save the df with a timestamp so we don't accidently overwrite it.
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
save(fss2.s1.vi, file=paste0("fss2_s1_vi_",timestamp,".RData"))
save(fss2.s1.visd, file=paste0("fss2_s1_visd_",timestamp,".RData"))
timestamp
load("fss2_s1_vi_20141105_1138.RData")
load("fss2_s1_visd_20141105_1138.RData")

#### s1 Boxplots ####

fss2.s1.vi.l <- melt(fss2.s1.vi, id=c("var_name","var_index"))

png('varImpALL.png', width = 960, height = 960)
bymedian <- with(fss2.s1.vi.l, reorder(var_index, value, median))
boxplot(value ~ bymedian, data = fss2.s1.vi.l,
        ylab = "Variable index", xlab = "% Increase MSE", 
        varwidth = TRUE,
        col = "lightgray", horizontal = TRUE)
dev.off()

fss2.s1.vi.median <- cast(fss2.s1.vi.l,var_name + var_index ~ ., value ='value', median)
colnames(fss2.s1.vi.median )[3] <- "median"

# sort the data so largest at the top
fss2.s1.vi.median <- fss2.s1.vi.median[with(fss2.s1.vi.median, order(-median)), ]

# save the interim dataframes so re-rerunning is easier
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
save(fss2.s1.vi.median, file=paste0("fss2_s1_vi_median_",timestamp,".RData"))
save(fss2.s1, file=paste0("fss2_s1_",timestamp,".RData"))
timestamp
load("fss2_s1_vi_median_20141105_1138.RData")
load("fss2_s1_20141105_1138.RData")

# Variable removal
# After the first 24 largest median values starts to flatten out so 
# we will take the top 40 to step 2. (33%)

# grab all variable names with median values > 0.8 = 33% of the data
fss2.s2.col <- fss2.s1.vi.median[fss2.s1.vi.median$median >= 0.85272628,][,1]
fss2.s2.col <- c("FSS_26Aug14",fss2.s2.col)
fss2.s2 <- fss2.s1[,colnames(fss2.s1) %in% fss2.s2.col]

#### Correlation plots ####
fss2.s2.col

# Precip "sum_1095_days","PPT_1981_2010","sum_365_days","sum_60_days","sum_180_days"
png('precip_cor.png')
pairs(fss2.s2[,fss2.s2.col[c(3,4,7,36)]],
      lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist)
dev.off()
# keep "sum_1095_days"

# Disturb [1] "PADISRSA_1YR" "PDISRCA_1YR"  "PDISRSA_1YR"  "PDISRCA_3YR"
pairs(fss2.s2[,fss2.s2.col[c(8,11,12,17)]],
      lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist)
# everything is coorelated
# keep "PDISRSA_1YR",

# Lithology/soils
# [1] "PALITHERODRCA"  "PALITHERODRSA"  "PASILTRCA"      "PACLAYRCA"      "PASILT_CLAYRCA" "PASANDRCA"      "MAKFACTRCA"     "PSILTRCA"      
# [9] "PCLAYRCA"       "PLITHERODRSA"   "PSANDRCA"
pairs(fss2.s2[,fss2.s2.col[c(2,5,14,15,18,20,21,23,24,31,41)]],
      lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist)
# almost everything is coorelated
# Keep "PALITHERODRCA", "PASILTRCA"

# Ownership "APOPRCA2010"  "POWNRCA_PRI"  
pairs(fss2.s2[,fss2.s2.col[c(13,35)]],
      lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist)
#When including the susceptibility breakouts we only get "POWNRCA_PRI" as a top variable

#Susceptibility
#[1] "PCONVEXRSA"   "PASLOPERCA"   "PACONVEXRCA"  "PASUSCEP4_DE" "PSUSCEP5_DE"  "PCONVEXRCA"
pairs(fss2.s2[,fss2.s2.col[c(10,16,25,26,28,30)]],
      lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist)
#almost everything is correlated
#Keep "PCONVEXRSA"

# Others
#  [1] "STRMPWR"     "XSLOPE_MAP"  "MIN_Z"       "upDist"      "LAT_RAW"     "afvArea"     "LONG_RAW"    "ARCASQM"     "ARSASQM"    
# [11] "PATYPEF"   
pairs(fss2.s2[,fss2.s2.col[c(6,9,22,29,32,33,34,38,39,40)]],
      lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist)
# udist and Long coorelated. arca and arsa correlated with each other.
# keeep "STRMPWR", "XSLOPE_MAP", "MIN_Z", "upDist", "LAT_RAW", "ARCASQM","PATYPEF"

keeps.s2 <- c("FSS_26Aug14",
              "sum_1095_days", 
              "PDISRSA_1YR","PALITHERODRCA", "PASILTRCA",
              "APOPRCA2010","POWNRCA_PRI",
              "PCONVEXRSA",
              "STRMPWR", "XSLOPE_MAP", "MIN_Z", "upDist", "LAT_RAW", "ARCASQM","PATYPEF")

#Further remove variables to reduce the influence of correlation on raising variable importance
fss2.s2 <- fss2.s2[,colnames(fss2.s2) %in% keeps.s2]
colnames(fss2.s2)

# remove any NAs
fss2.s2 <- data.frame(na.omit(fss2.s2))
#write.csv(fss2.s2, 'fss2_s2_data.csv')

#### Random Forest Step 2 ####
colnames(fss2.s2)

#DO we want to run this here?
# set.seed(42)
# fss2.s2.rf <- rf.modelSel(xdata=fss2.s2[,c(1:7,9:ncol(fss2.s2))], 
#                           ydata=fss2.s2[,"FSS_26Aug14"], 
#                           imp.scale="mir", r=c(0.5,0.10, 0.15,0.20,0.25,0.30,0.35,0.40,0.45, 0.5,0.55,0.60,0.75,0.80,0.85,0.90, 0.95),  
#                           final=TRUE, plot.imp=TRUE, parsimony=0.03, ntree=2000) 

# mtry value for use in the random forest function
mtry.fss2.s2 <- as.integer(((ncol(fss2.s2)-1) / 3),0)

# initialize the variable importance df
fss2.s2.vi <- data.frame(matrix(, nrow = ncol(fss2.s2)-1, ncol = 50))
fss2.s2.visd <- data.frame(matrix(, nrow = ncol(fss2.s2)-1, ncol = 50))

fss2.s2.col <- colnames(fss2.s2)
fss2.s2.col <- fss2.s2.col[!(fss2.s2.col == "FSS_26Aug14")]

#Run the random forest on the reduced set of variables
beg <- Sys.time()
set.seed(100)
for (i in 1:50) {
  fss2.s2.rf <- randomForest(FSS_26Aug14 ~ ., 
                             data = fss2.s2, 
                             ntree = 1000, 
                             keep.forest = TRUE, 
                             importance = TRUE)
  fss2.s2.vi[i,] <- importance(fss2.s2.rf, conditional = TRUE)
  fss2.s2.visd[,i] <- fss2.s2.rf$importanceSD
}
print(Sys.time() - beg)

# Add var names and index
fss2.s2.vi[,51]<- fss2.s2.col
fss2.s2.vi[,52]<-c(1:length(fss2.s2.col))
colnames(fss2.s2.vi)[51] <- "var_name"
colnames(fss2.s2.vi)[52] <- "var_index"
fss2.s2.visd[,51]<- fss2.s2.col
fss2.s2.visd[,52]<-c(1:length(fss2.s2.col))
colnames(fss2.s2.visd)[51] <- "var_name"
colnames(fss2.s2.visd)[52] <- "var_index"

#Save the dfs with a timestamp so we don't accidentally overwrite them
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
save(fss2.s2.vi, file=paste0("fss2_s2_vi_",timestamp,".RData"))
save(fss2.s2.visd, file=paste0("fss2_s2_visd_",timestamp,".RData"))
timestamp
load("fss2_s2_vi_20141106_0858.RData")
load("fss2_s2_visd_20141106_0858.RData")

#### s2 Boxplot ####

fss2.s2.vi.l <- melt(fss2.s2.vi, id=c("var_name","var_index"))

png('varImpALL_s2.png', width = 960, height = 960)
bymedian <- with(fss2.s2.vi.l, reorder(var_name, value, median))
par(yaxt="n",mar=c(5, 8, 4, 5))
boxplot(value ~ bymedian, data = fss2.s2.vi.l,
        ylab = "var name", xlab = "%InMSE", 
        varwidth = TRUE,
        col = "lightgray", horizontal = TRUE)
lablist.y<-levels(bymedian)
axis(2, labels = FALSE)
text(y = 1:15, par("usr")[1], labels = lablist.y, pos = 2, xpd = TRUE)
dev.off()

#### Partial dependence plots ####

#Partial dependence plots
for (i in 1:length(fss2.s2.vi$var_name)) {
  filename <- paste("partialPlot_", fss2.s2.vi$var_name[i], ".png",sep = "")
  png(filename, width = 960, height = 960)
  partialPlot(fss2.s2.rf, 
              fss2.s2, 
              x.var = fss2.s2.vi$var_name[i], 
              ylab = 'Mean FSS', 
              ylim = c(9,18))
  dev.off()  
}