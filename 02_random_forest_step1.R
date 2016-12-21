#This file needs to be run in 64 bit R

library(randomForest)
library(reshape)
library(plyr)

options(stringsAsFactors = FALSE)

ssn_RF_data <- read.csv("ssn_RF_data.csv")

# ----------------------------------------------------------- #
# FSS2 - Random forests excluding the physical habitat data ####
# ----------------------------------------------------------- #
# Step 1, per Genuer et al (2000), use brute force to identify important variables. 
# Since we have a large amount of variables random forest is run 50 times with a 
# very high number of trees per forest (ntree). This yields a distribution of importance scores.
# Removal is based on these distributions. We keep the variables with the highest scores.

bsti <- ssn_RF_data[, !colnames(ssn_RF_data) %in% c('SVN','STATION_KEY',
                                                       'rid','rid_LSN04')]

bsti$DATE <- as.POSIXct(bsti$DATE)

# #The below steps were necessary during initial processing and are not necessary now that we are 
# #using the modified dataset
# #remove NAs in response variable
# bsti <- bsti[(!is.na(bsti$BSTI)),]
# # remove any NAs - this has been taken care of prior to this but for good measure we'll leave it in
# #This removes the whole row where there is an NA. There are columns that are all NA due to scaling 
# bsti <- data.frame(na.omit(bsti))

# #Normalize the data to the same scale - Completed on 12/17 at 1618/1620. Did not affect rank order
# Will not normalize to same scale at this point in the data processing since it affects the 
# value of the %inc MSE because the scale of the response is so small. 
# bsti$LONG_RAW <- abs(bsti$LONG_RAW)
# bsti <- as.data.frame(lapply(bsti,function(x) {(x-min(x))/(max(x)-min(x))}))
# 
#This removes those variables where all the values are 0
bsti <- bsti[, setdiff(names(bsti), c("X2year_count_60_days",
                                               "X10year_count_60_days", 
                                               "X25year_count_60_days", 
                                               "X50year_count_60_days",
                                               "X100year_count_60_days", 
                                               "X10year_count_180_days", 
                                               "X25year_count_180_days",
                                               "X50year_count_180_days", 
                                               "X100year_count_180_days"))]
#Need to convert Inf values to 0
bsti[, grep("X", names(bsti))] <- as.data.frame(sapply(
  bsti[, grep("X",names(bsti))], function(x) {
    replace(x, is.infinite(x),0)
    }
  ))

#Soil characteristics are known complements of each other. We are going
#to select the size class representative of fine sediment <0.05mm
bsti <- bsti[, -grep('^SAND|^CLAY|^SILT_P', names(bsti))]

#COMP is a complement of EROD and has high correlation. Although EROD is 
#a component of the derivation of SUSCEP it maintains low correlation < 0.4
bsti <- bsti[, -grep('COMP', names(bsti))]

#SLOPE and Q0001E_adj are precursors to the calculated STRMPWR
#based on previous runs of the model with these variables used independently
#SLOPE and Q0001E_adj arrived at a model with lower AIC suggesting they
#produce a more likely model. Both may be run but for now we will use the
#component variables instead of the calculated variavle
for (RUN in 1:2) {
if (RUN == 1) {
  bsti.run <- within(bsti, rm(STRMPWR))
} else if (RUN == 2) {
  bsti.run <- within(bsti, rm(Q0001E_adj, XSLOPE_MAP))
}

#Normalize by maximum range
# melted <- melt(bsti[,names(bsti[,-c(grep('_P',names(bsti)),
#                                           which(names(bsti) %in%
#                                                   c('DATE')))])])
# 
# min.max <- ddply(melted, .(variable),
#                  summarize,
#                  min_val = min(value),
#                  max_val = max(value))
# bsti[,-c(grep('_P',names(bsti)),
#             which(names(bsti) %in%
#                     c('DATE')))] <- as.data.frame(
#                       lapply(bsti[,-c(grep('_P', names(bsti)),
#                                          which(names(bsti) %in% c('DATE')))],
#                              function(x) {((x) / (max(x)))*100}))
#OR standardize using mean and two times standard deviation to aid in coefficient interpretation
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

bsti.run <- (as.data.frame(stdpreds(bsti.run, bsti.run)))
#look <- lapply(bsti.run_std, hist)

# mtry value
mtry.bsti <- as.integer(((ncol(bsti.run) - 1) / 3), 0)

# initialize the variable importance df to store the importance scores from each run
bsti.vi <- data.frame(matrix(nrow = ncol(bsti.run) - 1, ncol = 50))
bsti.visd <- data.frame(matrix(nrow = ncol(bsti.run) - 1, ncol = 50))

bsti.col <- colnames(bsti.run)
bsti.col <- bsti.col[!(bsti.col == "BSTI")]


#Run the randomForest
# WARNING - Takes about 30 min
beg <- Sys.time()
print(beg)
for (i in 1:50) {
  set.seed(i)
  bsti.rf <- randomForest(BSTI ~ ., 
                             data = bsti.run,
                             mtry = mtry.bsti,
                             ntree = 2000, 
                             keep.forest = TRUE, 
                             importance = TRUE)
  bsti.vi[, i] <- bsti.rf$importance[, 1]
  bsti.visd[, i] <- bsti.rf$importanceSD
}
print(Sys.time() - beg)

# Add var names and index
bsti.vi[, 51]<- bsti.col
bsti.vi[, 52]<-c(1:length(bsti.col))
colnames(bsti.vi)[51] <- "var_name"
colnames(bsti.vi)[52] <- "var_index"
bsti.visd[, 51]<- bsti.col
bsti.visd[, 52]<- c(1:length(bsti.col))
colnames(bsti.visd)[51] <- "var_name"
colnames(bsti.visd)[52] <- "var_index"

#### sort the variables by median importance ####
bsti.vi.l <- melt(bsti.vi, id = c("var_name", "var_index"))

bsti.vi.median <- cast(bsti.vi.l, var_name + var_index ~ ., 
                          value ='value', median)
colnames(bsti.vi.median )[3] <- "median"

#Take the median sd as well
bsti.visd.l <- melt(bsti.visd, id = c("var_name", "var_index"))

bsti.visd.median <- cast(bsti.visd.l, var_name + var_index ~ ., 
                          value ='value', median)
colnames(bsti.visd.median )[3] <- "median_sd"

# sort the data so largest at the top
bsti.vi.median <- bsti.vi.median[with(bsti.vi.median, order(-median)), ]

## Save the precursors to s2 with a timestamp so we don't accidently overwrite it.
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
vi_median_name <- paste0("bsti_vi_median_", timestamp, ".RData")
bsti_name <- paste0("bsti_", timestamp, ".RData")
save(bsti.vi.median, file = vi_median_name)
save(bsti.run, file = bsti_name)
timestamp

}