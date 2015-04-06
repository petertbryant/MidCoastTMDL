library(reshape)
library(stringr)
library(plyr)
library(xlsx)
options(scipen = 100, stringsAsFactors = FALSE)

#### Bring in source data ####

#Bring in the cART variables for the CTSI sites
source('Rcode_CART_vars_for_CTSI.R')

CTSI_chars$FSS_May05_trans<-log10(as.numeric(CTSI_chars$FSS_May05)+1)
CTSI_chars$LONG_NHD_trans<-as.numeric(CTSI_chars$LONG_NHD)*-1   #can't do calculations on negative values for transformation code
CTSI_chars$PRECIPSITE_MM_trans<-log10(as.numeric(CTSI_chars$PRECIPSITE_MM) +1)
CTSI_chars$STRMPOWER_trans<-log10(as.numeric(CTSI_chars$STRMPOWER) +1)
CTSI_chars$MAFLOWU_trans<-log10(as.numeric(CTSI_chars$MAFLOWU) +1)
CTSI_chars$ELEV_FT_trans<-sqrt(as.numeric(CTSI_chars$ELEV_FT))
CTSI_chars$SLOPE_trans<-sqrt(as.numeric(CTSI_chars$SLOPE))
CTSI_chars$PCT_FSP3_trans<-asin(sqrt(round(as.numeric(CTSI_chars$PCT_FSP3),2)))

#Bring in data used as input to random forest 
#These data were used as input to Rcode_transformations.R
obs.complete <- read.csv("ssn_RF_data.csv")
obs.complete$SVN <- str_trim(obs.complete$SVN)

#we name this bugs.all instead of obs.vars as in Rcode_transformations.R
bugs.all <- obs.complete[,c('SVN','STATION_KEY','SITE_NAME','PREDATOR_Nov05_model','PREDATOR_Nov05_score','PREDATOR_Nov05_Condition',
                            'YEAR')]
bugs.all <- rename(bugs.all, c('YEAR' = 'Year_Sampled'))

#The biocriteria status as previously determined
bc_old <- read.csv('//deqhq1/tmdl/tmdl_wr/midcoast/models/sediment/target_development/R_output_biocriteria_status.csv')

#The original target values as output from the first time Ryan ran the CART model
targets <- read.csv('//deqhq1/tmdl/tmdl_wr/midcoast/models/sediment/target_development/final/R_output_bugs_CART_ALL_final_2013-06-15_trans.csv')

#Data submitted to EPA for 2012 Integrated Report - Has updated FSS values
FSS <- read.xlsx2('C:/users/pbryant/desktop/OE_Stress_Abunds_data quality_25feb15.xlsx', sheetName = 'all')

#These are the input data for the CART calculated for ALL samples (except CTSI sites)
chars <- read.csv('//deqhq1/tmdl/tmdl_wr/midcoast/data/benthicmacros/stationwork/r_inputs/r_input_cart_2013_06_15.csv')
chars$SVN <- str_trim(chars$SVN)

#data transformations so the calculated data will fit with the model predict function
chars$FSS_May05_trans<-log10(as.numeric(chars$FSS_May05)+1)
chars$LONG_NHD_trans<-as.numeric(chars$LONG_NHD)*-1   #can't do calculations on negative values for transformation code
chars$PRECIPSITE_MM_trans<-log10(as.numeric(chars$PRECIPSITE_MM) +1)
chars$STRMPOWER_trans<-log10(as.numeric(chars$STRMPOWER) +1)
chars$MAFLOWU_trans<-log10(as.numeric(chars$MAFLOWU) +1)
chars$ELEV_FT_trans<-sqrt(as.numeric(chars$ELEV_FT))
chars$SLOPE_trans<-sqrt(as.numeric(chars$SLOPE))
chars$PCT_FSP3_trans<-asin(sqrt(round(as.numeric(chars$PCT_FSP3),2)))

#The CART model as fit to the reference population
load("//DEQHQ1/TMDL/TMDL_WR/MidCoast/Models/Sediment/Target_Development/R_output_Ref_CART_pruned_model_2013-06-10_trans.RData")

#The source data used to build the original targets
bugs<- read.csv("//Deqhq1/tmdl/TMDL_WR/MidCoast/Documents/TMDL/BiocriteriaMethodsCART_DATA.csv", na.strings = c("",-9999,-9998), stringsAsFactors = FALSE)
colnames(bugs)[3]<-paste("LONG_NHD")
colnames(bugs)[6]<-paste("PREDATOR_Nov05_model")
colnames(bugs)[4]<-paste("LAT_NHD")
colnames(bugs)[17]<-paste("PCT_FSP3")
colnames(bugs)[30]<-paste("FSS_May05")
colnames(bugs)[31]<-paste("PREDATOR_Nov05_score")

# Data transformations to the original target data
bugs$FSS_May05_trans<-log10(bugs$FSS_May05+1)
bugs$LONG_NHD_trans<-bugs$LONG_NHD*-1   #can't do calculations on negative values for transformation code
bugs$PRECIPSITE_MM_trans<-log10(bugs$PRECIPSITE_MM +1)
bugs$STRMPOWER_trans<-log10(bugs$STRMPOWER +1)
bugs$MAFLOWU_trans<-log10(bugs$MAFLOWU +1)
bugs$ELEV_FT_trans<-sqrt(bugs$ELEV_FT)
bugs$SLOPE_trans<-sqrt(bugs$SLOPE)
bugs$PCT_FSP3_trans<-asin(sqrt(bugs$PCT_FSP3))

#Getting the reference data statewide from the original target data
ref<-subset(bugs, REF_SAMPLES == "1")
ref$BASINNAME<-factor(ref$BASINNAME,levels=unique(ref$BASINNAME))
ref$ERODRES_CAT<-factor(ref$ERODRES_CAT,levels=unique(ref$ERODRES_CAT))
ref<-ref[,c("FSS_May05","FSS_May05_trans","AREAWTMAP","MAFLOWU_trans","SLOPE_trans","PRECIPSITE_MM_trans","STRMPOWER_trans","LONG_NHD_trans","LAT_NHD","ELEV_FT_trans","PCT_FSP3_trans","ERODRES_CAT")]

#### Continue data input build out ####
#Combine CART characteristics to the samples used inlcuding new and old
bugs.all <- merge(bugs.all,chars[,c("SVN","Ref_Samples","NewSample","FSS_May05","FSS_May05_trans","AREAWTMAP","MAFLOWU_trans","SLOPE_trans","PRECIPSITE_MM_trans",
                                    "STRMPOWER_trans","LONG_NHD_trans","LAT_NHD","ELEV_FT_trans","PCT_FSP3_trans","ERODRES_CAT")],
                  by = 'SVN',all.x=TRUE)
#Pull in the CTSI CART predictors
bugs.all <- merge(bugs.all, CTSI_chars[,c("STATION_KEY","FSS_May05","FSS_May05_trans","AREAWTMAP","MAFLOWU_trans","SLOPE_trans","PRECIPSITE_MM_trans",
                                          "STRMPOWER_trans","LONG_NHD_trans","LAT_NHD","ELEV_FT_trans","PCT_FSP3_trans","ERODRES_CAT")], by = 'STATION_KEY', 
                  all.x = TRUE, suffixes = c('','.y'))
bugs.all[grep('CTSI',bugs.all$STATION_KEY),c("FSS_May05_trans","AREAWTMAP","MAFLOWU_trans","SLOPE_trans","PRECIPSITE_MM_trans",
                                          "STRMPOWER_trans","LONG_NHD_trans","LAT_NHD","ELEV_FT_trans",
                                          "PCT_FSP3_trans","ERODRES_CAT")] <- bugs.all[grep('CTSI',bugs.all$STATION_KEY),
                                                                                       c("FSS_May05_trans.y","AREAWTMAP.y","MAFLOWU_trans.y","SLOPE_trans.y",
                                                                                         "PRECIPSITE_MM_trans.y", "STRMPOWER_trans.y","LONG_NHD_trans.y",
                                                                                         "LAT_NHD.y","ELEV_FT_trans.y","PCT_FSP3_trans.y","ERODRES_CAT.y")]
#One CTSI station (CTSI_57) was assigned to Station 29898 - SVN: SLTZ0812IRM0004
bugs.all[bugs.all$SVN == 'SLTZ0812IRM0004',c("AREAWTMAP","MAFLOWU_trans","SLOPE_trans","PRECIPSITE_MM_trans",
                                             "STRMPOWER_trans","LONG_NHD_trans","LAT_NHD","ELEV_FT_trans",
                                             "PCT_FSP3_trans","ERODRES_CAT")] <- bugs.all[bugs.all$STATION_KEY == 29898,c("AREAWTMAP","MAFLOWU_trans",
                                                                                                                          "SLOPE_trans","PRECIPSITE_MM_trans",
                                                                                                                          "STRMPOWER_trans","LONG_NHD_trans","LAT_NHD",
                                                                                                                          "ELEV_FT_trans", "PCT_FSP3_trans","ERODRES_CAT")][2,]

#### Evalute biocriteria status ####
#Bring in Data Quality objective determinations
bugs.all <- merge(bugs.all, FSS[,c('Sample','Use.303d')], by.x = 'SVN', by.y = 'Sample', all.x = TRUE)
bugs.all[grep('SLTZ',bugs.all$SVN),'Use.303d'] <- 'Yes'
bugs.dqo <- bugs.all[bugs.all$Use.303d == 'Yes',]

#Average scores for same year
bugs.s <- aggregate(PREDATOR_Nov05_score ~ STATION_KEY + Year_Sampled + PREDATOR_Nov05_model, data=bugs.dqo, mean)
bugs.s <- bugs.s[bugs.s$Year_Sampled >= 2002,]

#Evaluate MWCF model samples
MWCF <- subset(bugs.s, PREDATOR_Nov05_model == "MWCF")
MWCF$sample_status <- ifelse(MWCF$PREDATOR_Nov05_score <= 0.85, "Impaired", 
                             ifelse(MWCF$PREDATOR_Nov05_score > 1.24, "Potential Concern", "Attaining"))

#Evaluate WCCP model samples
WCCP <- subset(bugs.s, PREDATOR_Nov05_model == "WCCP")
WCCP$sample_status <- ifelse(WCCP$PREDATOR_Nov05_score <= 0.78, "Impaired", 
                             ifelse(WCCP$PREDATOR_Nov05_score > 1.23, "Potential Concern","Attaining"))
bugs.s <- rbind(MWCF,WCCP)
bugs.s$SY <- paste(bugs.s$STATION_KEY, bugs.s$Year_Sampled)
bugs.all$SY <- paste(bugs.all$STATION_KEY, bugs.all$Year_Sampled)
bugs.all <- merge(bugs.all,bugs.s[,c('SY','sample_status')],all.x=TRUE)

#Prepare for status assignment
bugs.s$samples_tot <- 1
bugs.s <- merge(bugs.s,aggregate(samples_tot ~ STATION_KEY, data=bugs.s,  sum), by ="STATION_KEY")
bugs.s <- bugs.s[with(bugs.s, order(STATION_KEY, -Year_Sampled)), ]
bugs.s$samples_total <- bugs.s$samples_tot.y
bugs.s$samples_tot.y <- NULL
bugs.s$samples_tot.x <- NULL

# Rank observations starting with the most recent
bugs.s <- ddply(bugs.s, "STATION_KEY", transform, obs = seq_along(STATION_KEY))

# Status of sites with just one sample
samples1 <- subset(bugs.s, samples_total == 1)
samples1$biocriteria_status <- samples1$sample_status

# more than one sample and not attaining status
samples2 <- subset(bugs.s, samples_total > 1) #select statons w/ more than one sample
samples2.1 <- subset(samples2, obs == 1 & sample_status != "Attaining") # select most recent observation that is not attaining
samples2.1$biocriteria_status <- samples2.1$sample_status

# more than one sample and attaining
dummy <- subset(samples2, obs == 1 & sample_status == "Attaining")
samples2.2_pc <- subset(samples2, STATION_KEY %in% unique(dummy$STATION_KEY) & obs == 2 & sample_status != "Attaining") # select observation #2 that is not attaining
samples2.2_pc$biocriteria_status <- "Potential Concern"
samples2.2_att <- subset(samples2, STATION_KEY %in% unique(dummy$STATION_KEY) & obs == 2 & sample_status == "Attaining")# select observation #2 that is attaining

# Evaluate if the attaining samples outnumber earlier impaired samples
dummy2 <- subset(bugs.s, STATION_KEY %in% unique(samples2.2_att$STATION_KEY)) # select all samples from stations that we are making this evaulation for
dummy2$samples <- 1
dummy3  <- cast(dummy2, STATION_KEY~sample_status, value="samples", sum) #sum the sample status by station
dummy3$biocriteria_status <- ifelse(dummy3$Attaining > dummy3$Impaired,"Attaining","Impaired") # See if attaining samples outnumber impaired samples
dummy3 <-  dummy3[,c("STATION_KEY","biocriteria_status")]
samples2.2_att <- merge(samples2.2_att,dummy3, by = "STATION_KEY", all.x = TRUE)

#put it together
bugs.s <- rbind(samples1,samples2.1,samples2.2_pc,samples2.2_att) 
bc_status <-  bugs.s[,c("STATION_KEY","biocriteria_status")]

#Check against previous status determinations
bc_status <- merge(bc_status, bc_old, all.x = TRUE, by = 'STATION_KEY', suffixes = c("",".OLD"))

#Stations that now show impairment
newImp <- bc_status[which(bc_status$biocriteria_status != bc_status$biocriteria_status.OLD & bc_status$biocriteria_status == 'Impaired'),]

#Stations that now show attainment
newAtt <- bc_status[which(bc_status$biocriteria_status != bc_status$biocriteria_status.OLD & bc_status$biocriteria_status.OLD == 'Impaired' & bc_status$biocriteria_status == 'Attaining'),]

#bring in status, old status, correct FSS values (hopefully), original targets where calculated
bugs.all <- merge(bugs.all, bc_status[,c('STATION_KEY','biocriteria_status')], all.x = TRUE, by = 'STATION_KEY')
bugs.all <- merge(bugs.all, bc_old, all.x = TRUE, by = 'STATION_KEY',suffixes = c("",".OLD"))
bugs.all <- merge(bugs.all, FSS[,c('Sample','FSS')], by.x = 'SVN', by.y = 'Sample', all.x = TRUE)
bugs.all <- merge(bugs.all, targets[,c('SVN','Q75TH')], by = 'SVN', all.x = TRUE)
bugs.all[is.na(bugs.all$biocriteria_status),'biocriteria_status'] <- bugs.all[is.na(bugs.all$biocriteria_status),'biocriteria_status.OLD']

#the fss values for CTSI need to be added in
bugs.all[grep('SLTZ',bugs.all$SVN),'FSS'] <- bugs.all[grep('SLTZ',bugs.all$SVN),'FSS_May05.y']

#with the correct FSS values we now need to match column names
bugs.all$FSS_May05 <- as.numeric(bugs.all$FSS)
bugs.all$FSS_May05_trans<-log10(as.numeric(bugs.all$FSS_May05)+1)

#### Apply CART model ####
bugs.all$fss_pred_trans <- predict(ref.F.cart.prune, newdata=bugs.all)
bugs.all$fss_pred <- (10^bugs.all$fss_pred_trans)-1
bugs.all$fss_pred_resids_trans <- bugs.all$FSS_May05_trans - bugs.all$fss_pred_trans
bugs.all$fss_pred_resids <- bugs.all$FSS_May05 - bugs.all$fss_pred

ref$fss_pred_trans <- predict(ref.F.cart.prune, newdata=ref)
ref$fss_pred <- (10^ref$fss_pred_trans)-1
ref$fss_pred_resids_trans <- ref$FSS_May05_trans - ref$fss_pred_trans
ref$fss_pred_resids <- ref$FSS_May05 - ref$fss_pred

#### Assign CART Groups ####
cartgroups <- unique(ref$fss_pred)
cartgroups <- (as.data.frame(cartgroups))
cartgroups <- cartgroups[with(cartgroups, order(cartgroups)), ]
cartgroups

#create dummy vector
bugs.all$Group <- bugs.all$fss_pred

## Group 1 = 1.757594 = cartgroups[1]
cartg1 <- subset(bugs.all, fss_pred == cartgroups[1])
cartg1$Group <- "Group 1"

## Group 2 =  4.218690 = cartgroups[3]
cartg2 <- subset(bugs.all, fss_pred == cartgroups[3])
cartg2$Group <- "Group 2"

## Group 3 = 4.95204 = cartgroups[2]
cartg3 <- subset(bugs.all, fss_pred == cartgroups[2])
cartg3$Group <- "Group 3"

## Group 4 = 9.357024 = cartgroups[4]
cartg4 <- subset(bugs.all, fss_pred == cartgroups[4])
cartg4$Group <- "Group 4"

## Group 5 = 16.493012 = cartgroups[5]
cartg5 <- subset(bugs.all, fss_pred == cartgroups[5])
cartg5$Group <- "Group 5"

bugs.all <- rbind(cartg1,cartg2,cartg3,cartg4,cartg5)
rm(cartg1,cartg2,cartg3,cartg4,cartg5)

table(bugs.all$Group)

#bugs.all$FSS <- bugs.all$FSS_May05

#### Determine target benchmarks ####
quants<-c(.5,.75,.9,.95)
10^quantile(ref$fss_pred_resids_trans,probs=quants)
#   50%         75%         90%         95% 
#   1.018466 1.541863 2.175810 2.708002 

# (any site where the Observed FSS is 2.2X greater than the Predicted = Poor FSS condition)
boxplot(10^ref$fss_pred_resids_trans)

bugs.all$Q75TH <- round(bugs.all$fss_pred * 10^(quantile(ref$fss_pred_resids_trans,probs=.75)), digits = 0)
bugs.all$Q90TH <- round(bugs.all$fss_pred * 10^(quantile(ref$fss_pred_resids_trans,probs=.90)), digits = 0)

#### Determine sediment stressor status ####
#set up residual conditions 
bugs.all$fss_resid_cond <- ifelse(as.numeric(bugs.all$FSS) >= as.numeric(bugs.all$Q90TH), "Poor", ifelse(bugs.all$FSS <= bugs.all$Q75TH, "Good", "Fair"))

# Not impaired or Potential Concern
sed01 <- subset(bugs.all, biocriteria_status == "Attaining" | biocriteria_status == "Potential Concern")
sed01$sediment_resid_status <- sed01$biocriteria_status

# Impaired - Unknown Stressor
dummy01 <- subset(bugs.all, biocriteria_status == "Impaired" & fss_resid_cond == "Good")
dummy02 <- subset(bugs.all, biocriteria_status == "Impaired" & fss_resid_cond == "Fair")
sed02 <- rbind(dummy01, dummy02)
sed02$sediment_resid_status <- "Impaired - Unknown Stressor"
rm(dummy01,dummy02)

# Impaired - Sediment Stressor
sed03 <- subset(bugs.all, biocriteria_status == "Impaired" & fss_resid_cond == "Poor")
sed03$sediment_resid_status <- "Impaired - Sediment Stressor"

#rm(bugs.all)
sed_resid_stat <- rbind(sed01[,c('SVN','sediment_resid_status')],sed02[,c('SVN','sediment_resid_status')],sed03[,c('SVN','sediment_resid_status')])
bugs.all <- merge(bugs.all, sed_resid_stat, by = 'SVN')
rm(sed_resid_stat,sed01,sed02,sed03)

bugs.all.any.sed <- ddply(bugs.all, .(STATION_KEY), function(x) {ifelse(any(x$sediment_resid_status == "Impaired - Sediment Stressor"),x$anysed <- 1, x$anysed <- 0)})
bugs.all.any.sed <- rename(bugs.all.any.sed, c("V1" = "ansysed"))
bugs.all <- merge(bugs.all, bugs.all.any.sed, by = 'STATION_KEY', all.x = TRUE)

table(bugs.all$sediment_resid_status)

#### compare midcoast stations only ####
mc <- read.csv("E:/Github/KML-R-tool/midcoast_stns.csv")
mc2 <- bugs.all[bugs.all$STATION_KEY %in% mc$STATION_KEY,]
#mc3 <- merge(mc,mc2,by='SVN')
#mc3[mc3$FSS_26Aug14 != mc3$FSS,c('FSS_26Aug14','FSS')]

#mc4<-merge(mc3,bc_old,by.x='STATION_KEY.x',by.y='STATION_KEY',all.x=TRUE, suffixes = c('.NEW','.OLD'))
#View(arrange(mc4[which(mc4$biocriteria_status.NEW != mc4$biocriteria_status.OLD),],STATION_KEY.X,Year_Sampled))

#### Determine which stations have new data ####
bugs.F <- aggregate(PREDATOR_Nov05_score ~ STATION_KEY + Year_sampled + PREDATOR_Nov05_model, data=bugs, mean)
bugs.F$samples_tot <- 1
bugs.F <- merge(bugs.F,aggregate(samples_tot ~ STATION_KEY, data=bugs.F,  sum), by ="STATION_KEY")
bugs.F <- bugs.F[with(bugs.F, order(STATION_KEY, -Year_sampled)), ]
bugs.F.mc <- bugs.F[bugs.F$STATION_KEY %in% mc$STATION_KEY,]

bugs.all[,'NEW_STATION'] <- ifelse(!bugs.all$STATION_KEY %in% bugs.F$STATION_KEY,'NEW','OLD')
bugs.all[,'NEW_SAMPLE'] <- ifelse(!bugs.all$SVN %in% bugs$SVN,'NEW','OLD')
bugs.all[,'NEW'] <- ifelse(bugs.all$NEW_STATION == 'NEW' | bugs.all$NEW_SAMPLE == 'NEW','NEW','OLD')

bugs.all <- bugs.all[,!grepl('\\.y',names(bugs.all))]

write.csv(bugs.all,'allstns_new_status.csv')

mc2 <- bugs.all[bugs.all$STATION_KEY %in% mc$STATION_KEY,]

write.csv(mc2,'midcoast_new_status.csv')
