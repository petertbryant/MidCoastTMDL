#Control script to run random forest, correlations, model selection and model averaging

#source('R_code_combine_data.R')
#Inputs: Calculated watershed characteristics. Mostly stored in LSN04/Tables.mdb. Some stored in various .csv
#Outputs: VarNames_RF_V2.csv, ssn_RF_data.csv

#source('R_code_RF_FSS2.R')
#Inputs: VarNames_RF_V2.csv, ssn_RF_data.csv
#Outputs: fss2_s1_vi_20150720_1627.RData = First run variable importance scores for all 50 runs 
#         fss2_s1_visd_20150720_1627.RData = First run variable importance standard deviations
#         fss2_s1_vi_median_20150722_1035.RData = Median importance values as based on the 50 runs
#         fss2_s1_20150722_1035.RData = The edited data used as input to randomForest function

#Placeholder for breakout of correlation variable selection step
#Inputs: fss2_s1_vi_median_20150722_1035.RData, fss2_s1_20150722_1035.RData 
#Outputs: fss2_s2_data.csv

#source('Rcode_transformations.R')
#Inputs: ssn_RF_data.csv, fss_s2_data.csv, ssn object, 
#Outputs: ssn object with transformed and scaled variables, minmax.Rdata = Min and max used for scaling. To be used for descaling.

#source('funDredgeGLMSSN.R')
#Inputs: ssn object, funMuMInhelpers.R
#Output: rval.out = model.selection object

#source('funModelAvgGLMSSN.R')
#Inputs: rval.out
#Output: ret = averaging object with averaged regression coefficients


#deprecated source('Rcode_ssn.R')
#Inputs: ssn object, minmax.Rdata
#OUtput: final model object