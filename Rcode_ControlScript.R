#Control script to run random forest, correlations, model selection and model averaging

#Step 1:
#source('01_combine_data.R')
#Inputs: Calculated watershed characteristics. Mostly stored in LSN04/Tables.mdb. Some stored in various .csv
#Outputs: VarNames_RF_V2.csv, ssn_RF_data.csv

#### START AT STEP 2 FOR CURRENT MODEL RUNS ####

#RUN1: No STRMPWR, using Q0001E_adj and XSLOPE_MAP instead since STRMPWR = Q*SLOPE
#RUN2: No Q0001E_adj or XSLOPE_MAP, using STRMPWR instead to test for better predictor
RUN <- 1

#Step 2:
#source('02_random_forest_step1.R')
#Inputs: VarNames_RF_V2.csv, ssn_RF_data.csv
#Outputs: fss2_s1_vi_20150720_1627.RData = First run variable importance scores for all 50 runs 
#         fss2_s1_visd_20150720_1627.RData = First run variable importance standard deviations
#         fss2_s1_vi_median_20150722_1035.RData = Median importance values as based on the 50 runs
#         fss2_s1_20150722_1035.RData = The edited data used as input to randomForest function

#Step 2a:
#source('02a_random_forest_step1_plots.R')
#Inputs: fss2.s1.vi.l
#Outputs: Plot of variable importance by median with standard deviation

#Step 3:
#source('03_variable_selection.R')
#Inputs: fss2_s1_vi_median_20150722_1035.RData, fss2_s1_20150722_1035.RData 
#Outputs: fss2_s2_data.csv

#Step 3a:
#source('03a_variable_selection_plots.R')
#Inputs: fss.s2, fss.s2.col
#Outputs: Correlation plots to R console (unless png device is activated)

#Step 3b:
#source('03b_random_forest_step2.R')
#Inputs: fss2.s2
#Outputs: Really it's for the plot of the variables to advance to regression

#Step 4:
#source('04_transformations.R')
#Inputs: ssn_RF_data.csv, fss_s2_data.csv, ssn object, 
#Outputs: ssn object with transformed and scaled variables, minmax.Rdata = Min and max used for scaling. To be used for descaling.

#Step 5:
#source('05_backward_stepwise_regression.R')
#Inputs: ssn object, obs.fss2
#Output: saves regression objects to folder

#source('05a_forward_regression.R')
#Inputs: ssn object, vars from pkeep
#Output: saves regression objects to folder

#source('05b_model_selection.R)
#Inputs: ssn regression objects in folder
#Outputs: Results as csv file with dAIC



#source('funModelAvgGLMSSN.R')
#Inputs: rval.out
#Output: ret = averaging object with averaged regression coefficients


#deprecated source('Rcode_ssn.R')
#Inputs: ssn object, minmax.Rdata
#OUtput: final model object