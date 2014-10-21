# socentor subunit stream temperature modeling and prediction
# Aug 2014

library(SSN)
library(foreign)

# Function to standardize variables
stand <- function(x) { (x-mean(x))/(2*sd(x))}

# Function to standardize a prediction dataset based on the fitting dataset 
stdpreds <- function(newset,originalset) {
	xnames <- colnames(newset)
	sx <- matrix(rep(NA,ncol(newset)*nrow(newset)),nrow=nrow(newset))
	for(i in 1:ncol(newset)) {
		var <- with(originalset,get(xnames[i]))
		sx[,i] <- (newset[,i]-mean(var))/(2*sd(var))
		}
	colnames(sx) <- colnames(newset)
	return(sx)
}

# Get fixed effects and SEs from glmssn
ests <- function(x) {
	means <- round(x$estimates$betahat,3)
	ses <- round(sqrt(diag(x$estimates$covb)),3)
	output <- cbind(means,ses)
	colnames(output) <- c("Estimate","SE")
	return(output)
}

setwd("D:\\Cutthroat_Climate_Project\\GNLCC-temperature\\socentor")

# Load the ssn and all sets of prediction points
soc <- importSSN("socentor.ssn",predpts="preds")

# Create distance matrices
createDistMat(soc,o.write=T,predpts="preds")

# Create raw torgegram
soctg <- Torgegram(soc,"STREAM_AUG",nlag=20)
jpeg("socrawtorg.jpg")
plot(soctg, main = "Raw Data Torgegram") 
dev.off()

# Extract dataframe and fit basic aspatial model with un-standardized predictors
socdf <- getSSNdata.frame(soc)

soc.lm <- lm(STREAM_AUG ~ ELEV + CANOPY + SLOPE + PRECIP + CUMDRAINAG + Y_COORD + NLCD11PC + BFI + Air_Aug + Flow_Aug, data=socdf)
summary(soc.lm)

# Standardize continuous covariates and add factor for year
continuous <- socdf[,c(11:21)]
cont.s <- apply(continuous,2,stand)
colnames(cont.s) <- c("elev","canopy","slope","precip","drainage","lat","water","glacier","bfi","airtemp","flow")
socdf.s <- data.frame(socdf,cont.s)
socdf.s$yearf <- factor(socdf.s$SAMPLEYEAR)
socs <- putSSNdata.frame(socdf.s,soc,"Obs")

soc.lm2 <- lm(STREAM_AUG ~ elev + canopy + slope + precip + drainage + lat + water + bfi + airtemp + flow, data=socdf.s)

# Extract and save pure aspatial results
socAe <- summary(soc.lm2)$coefficients
backtrans <- socAe[-1,1:2]/(2*sapply(continuous[,-8],sd))
esttable <- cbind(rbind(socAe[1,1:2],backtrans),socAe)
rownames(esttable) <- rownames(socAe)
write.csv(esttable,"aspatialestimates.csv")

# Aspatial performance
predictA <- predict(soc.lm)
sqrt(mean((predictA-socdf.s$STREAM_AUG)^2))

# Examine correlations
library(ellipse)
cor(cont.s)
jpeg("corrplot.jpg")
plotcorr(cor(cont.s),type="lower")
dev.off()

# Get VIFs
library(car)
vif(soc.lm2)

# Standardize preds based on obs
socpreddf <- getSSNdata.frame(soc, "preds")
contpred <- socpreddf[,c(11:21)]
contpred.s <- stdpreds(contpred,continuous)
colnames(contpred.s) <- c("elev","canopy","slope","precip","drainage","lat","water","glacier","bfi","airtemp","flow")
socpreddf.s <- data.frame(socpreddf,contpred.s)
socpreddf.s$yearf <- factor(socpreddf.s$SAMPLEYEAR)
socs <- putSSNdata.frame(socpreddf.s,socs,"preds")

###########
# Model1  #
###########

starttime <- Sys.time()
soc1 <- glmssn(STREAM_AUG ~ elev + canopy + slope + precip + drainage + lat + water + bfi + airtemp + flow, socs, EstMeth= "ML", family="gaussian",CorModels = c("locID","yearf","Exponential.tailup","Exponential.taildown","Exponential.Euclid"), addfunccol = "afvArea")
elapsed <- Sys.time()-starttime
elapsed

# Extract predictions/ LOO CV predictions
soc1r <- residuals(soc1)
soc1rdf <- getSSNdata.frame(soc1r)

# Torgegram of fitted model
soc1t <- Torgegram(soc1r,"_resid.crossv_",nlag=20)
jpeg("soc1residtorg.jpg")
plot(soc1t, main= "Model1 Residuals Torgegram") 
dev.off()

#RMSPE of cross-validated predictions
sqrt(mean((soc1rdf[,"_CrossValPred_"]-soc1rdf$obsval)^2))
# 0.949

#RMSPE of fixed effects only
sqrt(mean((soc1rdf[,"_fit_"]-soc1rdf$obsval)^2))
# 2.236

# Null RMSPE
sqrt(mean((soc1rdf$obsval-mean(soc1rdf$obsval))^2))
# 3.545

#r2 of cross-validated predictions. 
cor(soc1rdf$obsval,soc1rdf[,"_CrossValPred_"])^2
# 0.928


# Summarize estimates and save
soc1e <- summary(soc1)$fixed.effects.estimates
backtrans <- soc1e[-1,2:3]/(2*sapply(continuous[,-8],sd))
esttable <- cbind(soc1e[,1],rbind(soc1e[1,2:3],backtrans),soc1e[,2:5])
write.csv(esttable,"soc1estimates.csv")

# Outlier analysis
resids <- soc1rdf[,"_resid.student_"]
outlie <- soc1rdf[abs(resids)>5,c("pid","obsval","_CrossValPred_","_resid.student_")]
outlie
# Only 2 records have abs studentized residuals >5 (5.3 and 5.1).

varcomp(soc1)
#               VarComp  Proportion
#1    Covariates (R-sq) 0.185721545
#2   Exponential.tailup 0.464684138
#3 Exponential.taildown 0.018729707
#4   Exponential.Euclid 0.188076898
#5                locID 0.055287992
#6                yearf 0.008398156
#7               Nugget 0.079101563

save.image(file = "D:/Cutthroat_Climate_Project/GNLCC-temperature/socentor/socSSNfit1.RData")


###############
# Predictions #
###############
soc1p1 <- predict(soc1,"preds")

pred1df <- getSSNdata.frame(soc1p1,"preds")

allpreds <- pred1df[,c("OBSPREDID","STREAM_AUG","STREAM_AUG.predSE")]
colnames(allpreds) <- c("OBSPREDID","predtemp","predtempse")

write.csv(allpreds,"soc1preds.csv",row.names=F)
write.dbf(allpreds,"soc1preds.dbf")

# Predicted and observed at fitting sites
predobs <- data.frame(soc1rdf$OBSPREDID,soc1rdf[,"_CrossValPred_"],soc1rdf$obsval)
colnames(predobs) <- c("obspredid","predicted","observed")
write.csv(predobs,"predobs.csv",row.names=F)

save.image(file = "D:/Cutthroat_Climate_Project/GNLCC-temperature/socentor/soc-pred.RData")


###########
# Model2 #
###########

# Aspatial model 

starttime <- Sys.time()
soc2 <- glmssn(STREAM_AUG ~ elev + canopy + slope + precip + drainage + lat + water + bfi + airtemp + flow, socs, EstMeth= "ML", family="gaussian",CorModels = c("locID","yearf"), addfunccol = "afvArea")
elapsed <- Sys.time()-starttime
elapsed

# Extract predictions/ LOO CV predictions
soc2r <- residuals(soc2)
soc2rdf <- getSSNdata.frame(soc2r)

soc2t <- Torgegram(soc2r,"_resid.crossv_",nlag=20)
jpeg("soc2residtorg.jpg")
plot(soc2t, main= "Model2 Residuals Torgegram") 
dev.off()

#RMSPE of cross-validated predictions
sqrt(mean((soc2rdf[,"_CrossValPred_"]-soc2rdf$obsval)^2))
# 1.092

#RMSPE of fixed effects only
sqrt(mean((soc2rdf[,"_fit_"]-soc2rdf$obsval)^2))
# 2.134

# Null RMSPEs
sqrt(mean((soc2rdf$obsval-mean(soc2rdf$obsval))^2))
# 3.545

#r2 of cross-validated predictions. 
cor(soc2rdf$obsval,soc2rdf[,"_CrossValPred_"])^2
# r2 = 0.905

soc2e <- summary(soc2)$fixed.effects.estimates
backtrans <- soc1e[-1,2:3]/(2*sapply(continuous[,-8],sd))
esttable <- cbind(soc1e[,1],rbind(soc1e[1,2:3],backtrans),soc1e[,2:5])

write.csv(esttable,"soc2estimates.csv")

varcomp(soc2)
#           VarComp  Proportion
#1 Covariates (R-sq) 0.387212009
#2             locID 0.543223150
#3             yearf 0.006954102
#4            Nugget 0.062610738
save.image(file = "D:/Cutthroat_Climate_Project/GNLCC-temperature/socentor/socSSNfit2.RData")


