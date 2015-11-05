#This file needs to be run in 64 bit R

library(randomForest)
library(reshape)
library(plyr)

source('funCorrelationPlots.R')

options(stringsAsFactors = FALSE)

vars <- read.csv("VarNames_RF_v2.csv")
bugs <- read.csv("ssn_RF_data.csv")

# ----------------------------------------------------------- #
# FSS2 - Random forests excluding the physical habitat data ####
# ----------------------------------------------------------- #
# Step 1, per Genuer et al (2000), use brute force to identify important variables. 
# Since we have a large amount of variables random forest is run 50 times with a 
# very high number of trees per forest (ntree). This yields a distribution of importance scores.
# Removal is based on these distributions. We keep the variables with the highest scores.

vars.fss2.s1 <- vars[vars$keep == 1, ]

fss2.s1 <- bugs[, colnames(bugs) %in% vars.fss2.s1$var]

fss2.s1$DATE <- as.POSIXct(fss2.s1$DATE)

# #The below steps were necessary during initial processing and are not necessary now that we are 
# #using the modified dataset
# #remove NAs in response variable
# fss2.s1 <- fss2.s1[(!is.na(fss2.s1$FSS_26Aug14)),]
# # remove any NAs - this has been taken care of prior to this but for good measure we'll leave it in
# #This removes the whole row where there is an NA. There are columns that are all NA due to scaling 
# fss2.s1 <- data.frame(na.omit(fss2.s1))

# #Normalize the data to the same scale - Completed on 12/17 at 1618/1620. Did not affect rank order
# Will not normalize to same scale at this point in the data processing since it affects the 
# value of the %inc MSE because the scale of the response is so small. 
# fss2.s1$LONG_RAW <- abs(fss2.s1$LONG_RAW)
# fss2.s1 <- as.data.frame(lapply(fss2.s1,function(x) {(x-min(x))/(max(x)-min(x))}))
# 
#This removes those variables where all the values are 0
fss2.s1 <- fss2.s1[, setdiff(names(fss2.s1), c("X2year_count_60_days",
                                               "X10year_count_60_days", 
                                               "X25year_count_60_days", 
                                               "X50year_count_60_days",
                                               "X100year_count_60_days", 
                                               "X10year_count_180_days", 
                                               "X25year_count_180_days",
                                               "X50year_count_180_days", 
                                               "X100year_count_180_days"))]
#Need to convert Inf values to 0
fss2.s1[, grep("X", names(fss2.s1))] <- as.data.frame(sapply(
  fss2.s1[, grep("X",names(fss2.s1))], function(x) {
    replace(x, is.infinite(x),0)
    }
  ))

#Normalize by maximum range
melted <- melt(fss2.s1[,names(fss2.s1[,-c(grep('_P',names(fss2.s1)),
                                          which(names(fss2.s1) %in% 
                                                  c('DATE')))])])
min.max <- ddply(melted, .(variable), 
                 summarize, 
                 min_val = min(value), 
                 max_val = max(value))
fss2.s1[,-c(grep('_P',names(fss2.s1)),
            which(names(fss2.s1) %in% 
                    c('DATE')))] <- as.data.frame(
                      lapply(fss2.s1[,-c(grep('_P', names(fss2.s1)), 
                                         which(names(fss2.s1) %in% c('DATE')))], 
                             function(x) {((x - min(x)) / (max(x) - 
                                                             min(x)))*100}))

# mtry value
mtry.fss2.s1 <- as.integer(((ncol(fss2.s1) - 1) / 3), 0)

# initialize the variable importance df to store the importance scores from each run
fss2.s1.vi <- data.frame(matrix(nrow = ncol(fss2.s1) - 1, ncol = 50))
fss2.s1.visd <- data.frame(matrix(nrow = ncol(fss2.s1) - 1, ncol = 50))

fss2.s1.col <- colnames(fss2.s1)
fss2.s1.col <- fss2.s1.col[!(fss2.s1.col == "FSS_26Aug14")]

#Run the randomForest
# WARNING - Takes about 30 min
beg <- Sys.time()
print(beg)
set.seed(100)
for (i in 1:50) {
  fss2.s1.rf <- randomForest(FSS_26Aug14 ~ ., 
                             data = fss2.s1,
                             mtry = mtry.fss2.s1,
                             ntree = 2000, 
                             keep.forest = TRUE, 
                             importance = TRUE)
  fss2.s1.vi[, i] <- fss2.s1.rf$importance[, 1]
  fss2.s1.visd[, i] <- fss2.s1.rf$importanceSD
}
print(Sys.time() - beg)

# Add var names and index
fss2.s1.vi[, 51]<- fss2.s1.col
fss2.s1.vi[, 52]<-c(1:length(fss2.s1.col))
colnames(fss2.s1.vi)[51] <- "var_name"
colnames(fss2.s1.vi)[52] <- "var_index"
fss2.s1.visd[, 51]<- fss2.s1.col
fss2.s1.visd[, 52]<- c(1:length(fss2.s1.col))
colnames(fss2.s1.visd)[51] <- "var_name"
colnames(fss2.s1.visd)[52] <- "var_index"

# Save the df with a timestamp so we don't accidently overwrite it.
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
save(fss2.s1.vi, file = paste0("fss2_s1_vi_", timestamp, ".RData"))
save(fss2.s1.visd, file = paste0("fss2_s1_visd_", timestamp, ".RData"))
timestamp

#### sort the variables by median importance ####
fss2.s1.vi.l <- melt(fss2.s1.vi, id = c("var_name", "var_index"))

fss2.s1.vi.median <- cast(fss2.s1.vi.l, var_name + var_index ~ ., 
                          value ='value', median)
colnames(fss2.s1.vi.median )[3] <- "median"

#Take the median sd as well
fss2.s1.visd.l <- melt(fss2.s1.visd, id = c("var_name", "var_index"))

fss2.s1.visd.median <- cast(fss2.s1.visd.l, var_name + var_index ~ ., 
                          value ='value', median)
colnames(fss2.s1.visd.median )[3] <- "median_sd"

# sort the data so largest at the top
fss2.s1.vi.median <- fss2.s1.vi.median[with(fss2.s1.vi.median, order(-median)), ]

## Save the precursors to s2 with a timestamp so we don't accidently overwrite it.
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
save(fss2.s1.vi.median, file=paste0("fss2_s1_vi_median_", timestamp, ".RData"))
save(fss2.s1, file=paste0("fss2_s1_", timestamp, ".RData"))
timestamp

