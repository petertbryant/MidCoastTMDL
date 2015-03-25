diff.SVN <- c('05007MIX','05009MIX','05010MIX','12001MCB','12002MCB','12003MCB','12004MCB','12005MCB','12006MCB','12007MCB','12008MCB','12010MCB','12013MCB','12014MCB','12016MCB','12018MCB','12019MCB','12020MCB','12021MCB','12022MCB','12023MCB','12024MCB','12025MCB','12026MCB','12027MCB','12028MCB','12029MCB','12031MCB','12032MCB','12033MCB','12035MCB','12036MCB','12037MCB','12038MCB','12039MCB','12040MCB')

targets <- read.csv('//deqhq1/tmdl/tmdl_wr/midcoast/models/sediment/target_development/final/R_output_bugs_CART_ALL_final_2013-06-15_trans.csv')

library(xlsx)
sheet1 <- read.xlsx('//Deqhq1/tmdl/TMDL_WR/MidCoast/Data/BenthicMacros/Raw_From_Shannon/Stressor ID_recaclulation QA Check_FSS output_26aug14_RM.xlsx', sheetName='Sheet1')
orig <- read.xlsx('//Deqhq1/tmdl/TMDL_WR/MidCoast/Data/BenthicMacros/2012_Macro_collection/Mid_Coast_2012_MacroSampling_2013_0422.xlsx', sheetName = 'RESULTS')

x <- obs.complete[obs.complete$SVN %in% diff.SVN,]
x <- merge(x, impaired, by = 'STATION_KEY', all.x = TRUE)
x <- merge(x, targets[,c('SVN','Q75TH',"biocriteria_status")], by = 'SVN', all.x = TRUE)
x <- merge(x, sheet1[,c('SVN','FSS')], by = 'SVN', all.x = TRUE)
x <- merge(x, orig[,c('SVN','FSS')], by = 'SVN', all.x = TRUE, suffixes = c('.sheet1','.orig'))

write.csv(x[,c('SVN','FSS_26Aug14','STATION_KEY','Q75TH','biocriteria_status','FSS.sheet1','FSS.orig')],'FSSreconcile.csv')

#Are there other samples
y <- x[which(x$biocriteria_status == 'Impaired' & x$FSS.sheet1 < x$Q75TH),c('SVN','STATION_KEY','SITE_NAME','FSS_26Aug14','Q75TH','biocriteria_status','FSS.sheet1','FSS.orig')]

obs.complete[obs.complete$STATION_KEY %in% y$STATION_KEY,c('SVN','STATION_KEY','SITE_NAME','FSS_26Aug14')]
