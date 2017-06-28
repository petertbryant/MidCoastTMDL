MidCoastTMDL
============

This folder contains all the code for model development, target/Loading Capacity
calculations and predictions. 

# Model development

## 01_combine_data.R
This file takes all the input data that has been calculated by Ryan Michie and
myself, Peter Bryant. It also handles consistent name formatting and some record specific
data corrections. The PPRCA and PPRSA Zonal Stats are the result of calculations
to make as many of the land area summarizations calculated to pour point specific
watersheds for each observation point (hence the PP for pour point). Because
many of the calculations were already completed to the catchments I do a subtraction
of the lowest RCA and then add in the pour point catchment area in order to complete
the calculations. 

The output of this file is ssn_RF_data.csv which is used in many subsequent scripts.

## 02_random_forest_step1.R
This file takes the input data and fits the random forest models to derive the
variable importance scores. Several different options were considered to filter
input variables ahead of time to reduce potential and known representation of
the same information from different variables. The standardization of the variables
also occurs in this script so important to take note so back-transformation can be
completed.

The outputs of this file are time-stamped files for the variable importance and
the standardized data.

### 02a_random_forest_step1_plots.R
This file was used to make the plot of the random forest variable importance
scores that is included in the source anlaysis methods document.

## 03_variable selection
Because I was considering the use of Stream power versus Slope and Discharge I ran
two versions from Script 02 and this script is set up to run the variable selection 
algorithm on both outputs. This removes all correlated variables save for the 
correlated variables with the highest importance. The intent here is to reduce as
much known colinearity between input variables to the SSN multiple regression. The
output of this script is a list with SSN objects with each of the potential sets
of variables.

### 03a_variable_selection_plots.R
This is code used to generate plots for inclusion in the source anlaysis methods
document.

#### 03b_random_forest_step2.R
This is old code that was used to make partial dependence plots and calcuate the 
R2 for the random forest model.

#### 04_transformations.R
This file was actually circumvented and included in Script 03 to streamline the 
evaluation of both candidate sets of variables.

## 05_backward_stepwise_regression.R
This script fits several models implementing a backward stepwise regression
algorithm and saves each one as a separate object.

#### 05a_forward_regression.R
This was an attempt to also include forward regression as a model fitting technique.
This method was not ultimately used.

## 05b_model_selection.R
This step calculates model fitting statistics including RMSE and AIC and Cross
Validation stats. It uses the functions_custom.R script which includes functions
I built around implementing a k-fold cross validation that instead of randomly selecting
which observations to subset or selecting one at a time, will drop all repeated 
observations at each sampling location and predict at those for calculating the 
error statistics.

#### 05c_Refit_with_factor.R
This was an attempt at improving model fit. 
This method was not ultimately used.

## 05d_Autocorrelation_selection.R
This file comprehensively fit all possible combinations of autocorrelation models

### 05e_ParameterUncertaintyAnalysis.R
This script was an attempt at explicitly incorporating uncertainty in parameter
estimations into the model outputs but ultimately is not being reported out at
this time.

## 06_scenarios.R
This script takes the best model fit and creates the reference (or loading capacity)
scenario that is then used to identify which biocriteria impaired sites have sediment
as a stressor. I am also using this script to rearrange the equation to allow
rainfall to vary in order to generate plots and equations for use as an implementation
tool.

#### 06a_Target_test_RUN2.R
OLD. Could probably just be deleted

#### 06a_Target_test_RUN3.R
OLD. Could probably just be deleted.

### 06_AutoMapGen.R
This was an attempt at making it easy to output maps for each of the impaired
watersheds of all of the anthropogenic activities for which we had data. This
is not currently being used but could potentially be. It's just that the map
formatting is not great.

## 07_Model_performance_evaluation.R
This script runs all the model checking functions that were included in the 
Source Analysis Model Development document.

## 08_Calculate_precip_new_sites.R
This is old code that I got from Ryan that I'm not actually using. This attempt
eventually morphed into Script 09.

## 09_impaired_ws_preds.R
This script uses the impaired watershed SSN (LSN07_IndianCreek) and generates
predictions using the dense prediction network that was populated with the model
variable data using the python scripts in the Watershed Characteristics/Scripts
folder. The output is a csv that can be used in ArcGIS to visualize the predictions
of where we may expect sediment to be a stressor.

### autocorfunPlots.R
This script was used to generate a chart for an earlier version of the Source 
Analysis Methods document.

#### biocriteria_add_stns.R
This is an old file and is no longer applicable to current efforts.

### FSS_correction_check_and_track.R
This script reconciled changes based on updates Shannon made to the FSS calculation.
This is old and the changes have been incorporated into current calculations.

#### PRECIP_rollSum_ImpairedStations.R
This script was an early development script and is now not used.

### Rcode_Disturbance.R
This code implemented the additional calculations necessary to make the 
pour point specific calculations.

### Rcode_Disturbance_01202016_revisions.R
This is code to update the watershed characteristic database with updated disturbance values

### Rcode_pprca_combine.R
This code builds the table to insert into the waterhsed characteristic database.

#### station_export_for_google_earth.R
At one point I wanted to provide a google earth layer for the TWG to be able
to look at station locations. Ultimately this was not used.

## functions_custom.R
This script compiles all the modifications I made to the SSN cross validation
functions to work with our unique case of multiple samples at a location

# dma_unique.R
This script identifies the DMAs that are covered in each impaired watershed.