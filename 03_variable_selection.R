library(randomForest)
library(reshape)
library(plyr)

options(stringsAsFactors = FALSE)

ssns <- list()
obs.bsti.lst <- list()
obs.vars.lst <- list()

for (j in 1:2) {
  if (j == 1) {
      #RUN1
  vi_median_name <- "bsti_vi_median_20161215_1347.Rdata"
  bsti_name <- "bsti_20161215_1347.Rdata"

  load(vi_median_name)
  load(bsti_name)
  } else {
      #RUN2
  vi_median_name <- "bsti_vi_median_20161215_1412.Rdata"
  bsti_name <- "bsti_20161215_1412.Rdata"
  
  load(vi_median_name)
  load(bsti_name)
  }

  bsti <- bsti.run
  #### Variable selection ####
  # Values drop off and then level out. Arbitrarily going with 50% of the variables.
  # grab all variable names with median values > 1.004880e-04 = 50% of the data
  # This 50% of the data reflects 50% of the original list of variables prior to scaling
  # Scaling had the effect of dropping variables that were all 0s anyway.
  bsti.s2.col <- bsti.vi.median[1:ceiling(nrow(bsti.vi.median)) / 2, ][, 1]
  bsti.vi.median <- bsti.vi.median[1:ceiling(nrow(bsti.vi.median)) / 2, ]
  #bsti.s2.col <- c("FSS_26Aug14",(bsti.vi.median[,'var_name']))
  bsti.s2 <- bsti[, colnames(bsti) %in% c(bsti.s2.col,'HDWTR')]
  
  # bsti.s2.col <- vars[vars$var %in% names(bsti.s2),]
  # #bsti.s2.col <- bsti.s2.col[bsti.s2.col$var != 'FSS_26Aug14',]
  # bsti.s2.col <- merge(bsti.s2.col, bsti.vi.median[, c('var_name','median')], 
  #                      by.x = 'var', by.y = 'var_name', all.x = TRUE)
  # bsti.s2.col <- arrange(bsti.s2.col, desc(median))
  
  
  # # #All together
  correlation_threshold <- 0.4
  pcor <- cor(bsti[,setdiff(bsti.vi.median$var_name,"DATE")])
  pnames <- attr(pcor, "dimnames")[[1]]
  pkeep <- pnames
  tracking <- vector("list", length(pnames))
  names(tracking) <- pnames
  for (i in length(pnames):1) {
    tracking[[i]] <- pcor[order(abs(pcor[,i])),i]#[apply(pcor, MARGIN = 1, FUN = function(x) abs(x) >= correlation_threshold)[,i]]
    if (any(round(abs(pcor[i,][-i]),2) >= correlation_threshold)) {
      if (i != 1) {
        pkeep <- pkeep[-i]
        pcor <- pcor[-i,-i,drop=FALSE]
      }
    } 
  }
  
  #Further remove variables to reduce the influence of correlation on raising variable importance
  #bsti.s2 <- bsti.s2[,colnames(bsti.s2) %in% keeps.s2]
  bsti.s2 <- bsti.s2[, c(pkeep,'HDWTR')]
  colnames(bsti.s2)
  
  # remove any NAs
  bsti.s2 <- data.frame(na.omit(bsti.s2))
  
  # -----------------------------------------------------------
  # SSN
  library(SSN)
  library(stringr)
  library(MASS)
  library(plyr)
  library(reshape2)
  obs.complete <- read.csv("ssn_RF_data.csv")
  obs.complete$SVN <- str_trim(obs.complete$SVN)
  #obs.bsti <- read.csv('fss2_s2_data_testing.csv')
  #obs.bsti <- within(obs.bsti, rm(X))
  #load('fss2_s2_scaled.Rdata')
  obs.bsti <- bsti.s2
  
  #in order to separate the obs data frame we have to build a new ssn object
  #this is done in arcGIS
  #but first we need to know which obs to put in which data frame
  #Split the observations up to generate an estimate and prediction set
  # lng_rng <- seq_along(obs.complete$SVN)
  # set.seed(100)
  # lng_est <- sample(lng_rng, 500)
  # set.seed(100)
  # lng_prd <- sample(lng_rng[-1 * lng_est], 260)
  # #Generate character vector for use in subsetting in arcgis context
  # paste(obs.complete[lng_est, 'SVN'], collapse = "','")
  # paste(obs.complete[lng_prd, 'SVN'], collapse = "','")
  # #Actually subset the data itself too
  # obs_est <- obs.complete[lng_est, ]
  # obs_prd <- obs.complete[lng_prd, ]
  
  ssn1 <- importSSN("//deqhq1/TMDL/TMDL_WR/MidCoast/Models/Sediment/SSN/LSN06/lsn.ssn",
                    predpts = "preds", o.write = TRUE)
  # ssn1 <- importSSN("//deqhq1/TMDL/TMDL_WR/MidCoast/Models/Sediment/SSN/LSN07/lsn.ssn",
  #                   predpts = "preds", o.write = TRUE)
  #only needs to be run once - RAN 07-01-2016
  #createDistMat(ssn1, o.write = TRUE)
  # ssn2 <- importSSN("//deqhq1/TMDL/TMDL_WR/MidCoast/Models/Sediment/SSN/LSN06/lsn.ssn",
  #                   predpts = "preds", o.write = TRUE)
  #ssn1 <- importSSN('C:/users/pbryant/desktop/midcoasttmdl-gis/revisedssn/lsn05/lsn.ssn', o.write = TRUE)
  obs<- getSSNdata.frame(ssn1, Name = "Obs")
  
  #Put together a data frame of the selected predictors and identifying information
  vars <- c("STATION_KEY", "SVN", "DATE","YEAR",'BSTI')
  obs.complete.vars <- obs.complete[,vars]
  obs.complete.vars <- cbind(obs.complete.vars, obs.bsti)
  
  #obs.vars <- obs.vars[!is.na(obs.vars$STRMPWR),]
  
  #RID crosswalk by SVN
  #rid.cross <- merge(obs[,c('STATION_KE','SVN','rid','pid')],
  #       obs.complete[,c('SVN','rid','pid')],by='SVN',suffixes=c('.NEW','.OLD'))
  
  #Output for GIS display
  #write.csv(obs.vars[,c("SVN","rid","sum_1095_days","FSS_26Aug14", 
  #                      "PDISRSA_1YR", "POWNRCA_PRI", "PALITHERODRCA", 
  #                      "PASILTRCA", "DAPOPRCA2010" )], 'final_model_vars.csv')
  
  #don't run when going through. only for data exploration
  #obs.vars <- arrange(obs.vars, STATION_KEY, desc(YEAR))
  #obs.vars.sub <- obs.vars[!duplicated(obs.vars$STATION_KEY),]
  
  #### Response TRANSFORMATION####
  names(obs.complete.vars)
  ####FSS_26Aug14####
  # hist(obs.vars.sub$FSS_26Aug14)
  # plot(density(obs.vars.sub$FSS_26Aug14))
  # hist(log10(obs.vars.sub$FSS_26Aug14))
  # plot(density(log10(obs.vars.sub$FSS_26Aug14)))
  # shapiro.test(log10(obs.vars.sub$FSS_26Aug14))
  # ks.test(log10(obs.vars.sub$FSS_26Aug14), "pnorm")
  #log should be just fine
  obs.complete.vars$log10_BSTI <- log10(obs.complete.vars$BSTI)
  
  #### Variable Scaling #### 
  # melted <- melt(obs.vars[,c(names(obs.bsti),'log10_FSS_26Aug14')])
  # min.max <- ddply(melted, .(variable), 
  #                  summarize, 
  #                  min_val = min(value), 
  #                  max_val = max(value))
  # obs.vars[,c(names(obs.bsti),'log10_FSS_26Aug14')] <- as.data.frame(lapply(
  #   obs.vars[,c(names(obs.bsti),'log10_FSS_26Aug14')], function(x) {(x-min(x))/
  #       (max(x)-min(x))}))
  # save(min.max, file = 'minmax.Rdata')
  # max_log10_bsti <- max(obs.complete.vars$log10_BSTI, na.rm = TRUE)
  # 
  # obs.complete.vars$log10_BSTI <- obs.complete.vars$log10_BSTI / max_log10_bsti * 100
  
  #Merge the selected predictors with the critical columns by SVN so they match
  #to the correct sample
  obs.vars <- merge(obs[,c("SVN","rid", "ratio", "locID", "netID", "pid", "upDist",
                           "afvArea", "HU_6_NAME", "HU_8_NAME", "HU_10_NAME", 
                           "HU_12_NAME", "HU_08", "NHDHigh", "LONG_RAW", 
                           "NHDh_Reach", "NHDP21_Rea", "NHDP12_COM", "HU_10", 
                           "HU_12")],
                    obs.complete.vars, 
                    by = "SVN",
                    all.x = TRUE)
  
  #### Put the data back into the ssn ####
  #Now that we have the transformed variables we put them back in the 
  #SSN object
  obs.vars$HDWTR <- as.factor(obs.vars$HDWTR)
  levels(obs.vars$HDWTR) <- c("0","1")
  obs.vars <- obs.vars[match(obs$pid, obs.vars$pid),]
  row.names(obs.vars) <- obs.vars$pid
  ssn1 <- putSSNdata.frame(obs.vars, ssn1, Name = 'Obs')
  
  obs <- getSSNdata.frame(ssn1, Name = 'Obs')
  preds <- getSSNdata.frame(ssn1, Name = "preds")
  # preds <- getSSNdata.frame(ssn1, Name = "obs_prd")
  
  pid.order <- preds$pid
  preds <- rename(preds, c('STATION_KE' = "STATION_KEY"))
  obs_sub <- obs[,c('STATION_KEY',pkeep,'log10_BSTI','HDWTR')]
  obs_sub <- obs_sub[order(obs_sub$STATION_KEY, obs_sub$log10_BSTI, decreasing = TRUE),]
  obs_sub <- obs_sub[!duplicated(obs_sub$STATION_KEY),]
  preds.vars <- merge(preds, obs_sub, by = 'STATION_KEY', all.x = TRUE)
  # preds.vars <- merge(preds[,c("SVN","rid", "ratio", "locID", "netID", "pid", "upDist",
  #                          "afvArea", "HU_6_NAME", "HU_8_NAME", "HU_10_NAME", 
  #                          "HU_12_NAME", "HU_08", "LONG_RAW", "NHDHigh",
  #                          "NHDh_Reach", "NHDP21_Rea", "NHDP12_COM", "HU_10", 
  #                          "HU_12")],
  #                   obs.complete.vars, 
  #                   by = "SVN",
  #                   all.x = TRUE)
  preds.vars <- preds.vars[match(pid.order,preds.vars$pid),]
  row.names(preds.vars) <- preds.vars$pid
  # ssn1 <- putSSNdata.frame(preds.vars, ssn1, Name = "obs_prd")
  ssn1 <- putSSNdata.frame(preds.vars, ssn1, Name = "preds")
  #save the ssn object to the github folder
  #writeSSN(ssn1, filename = 'bugs.ssn')
  
  ssns[[j]] <- ssn1
  obs.bsti.lst[[j]] <- obs.bsti
  obs.vars.lst[[j]] <- obs.vars
}

rm(list = setdiff(ls(), c("ssns", "obs.bsti.lst", "obs.vars.lst")))

