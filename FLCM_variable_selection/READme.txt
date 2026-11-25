The following funtions are used.
########## function for preprocessing original covariates with measurement error #####
function:preprocess 
funtioncall: preprocess(df)
input: df #dataframe in long format with following columns
id   : id of the subjects.
time : timepoints vectorized.
rest of the  columns are covariates
######################################################################################
### function for performing Variable Selection ###
function:FLCM.select 
funtioncall: FLCM.select(y,mydata,nbasis1,nbasis2,cvl)  
input: (y,mydata,nbasis1,nbasis2,cvl)  
y: the response vectorized (as in long format)
mydata: dataframe in long format column 1=id, column 2 = time, rest covariates
nbasis1: Number of basis functions for intercept
nbasis2: Number of basis functions for regression functions 
cvl    : Length of crossvalidation grid for tuning parameter \psi
######################################################################################
"FLCM_selection.html" illustrates appication in calcium data
######################################################################################
"varselect.R" contains the main functions and needs to be sourced 
######################################################################################
calcium dataset included : "calcium.csv"        
######################################################################################



