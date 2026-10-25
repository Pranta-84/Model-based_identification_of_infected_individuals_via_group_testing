# Cleaning the R environment
rm(list = ls())
#### All the packages needed ####
library(glmnet)
library(reticulate)
library(R.matlab)
library(tibble)

############### Set up to run MATLAB functions in R ############################
# Telling which version of python to use in the conda environment "r-matlab311"
Sys.setenv(RETICULATE_PYTHON = "/home/bilder/pranta/.conda/envs/r-matlab311/bin/python")
# Removing any active virtual environment that may conflict with the conda environment
Sys.unsetenv("VIRTUAL_ENV")
# Seeing which version of python is used (we would need python 3.11)
py_config()
# Importing numpy
np <- reticulate::import("numpy")
# Importing "matlab.engine" module of python; this will allow R to communicate with MATLAB
matlab <- import("matlab.engine")
# Starting a matlab session in the background (-nosplash supresses the MATLAB splash screen and -nodesktop runs MATLAB without GUI)
eng <- matlab$start_matlab("-nosplash -nodesktop")
# Adding the folder that contains MATLAB functions
eng$addpath("~/PBestMatlabCode/", nargout = as.integer(0))
# Importing the "matlab" package of python
matlab_types <- import("matlab")
# Explicitly telling which python version to use
#use_python("/home/bilder/pranta/.conda/envs/r-matlab311/bin/python", required = TRUE)
################################################################################

############################## Importing Tapestry functions ####################
sbl <- import("sbl")
l1ls <- import("l1ls")
nnompcv <-  import('nnompcv')

####################### Reading pooling matrix #########################
#M <- readMat("PBestMatlabCode/poolingMatrix.mat")
#M <- M$poolingMatrix
#M <- read.table("optimized_M_93_1240_kirkman.txt", header = F)
#M <- M[, 1:961]
#M <- matrix(as.numeric(as.matrix(M)), nrow = nrow(M))


################## 188x282 pooling matrix ################################
# M4.RDS is same as "P-BEST188x282.csv".
## The original function uses M4.RDS so I keeping it here.
## I am not replacing it by "PBEST188x282.csv" to avoid any kind of errors.
## This R program calls MATLAB and python; changing the file may cause conversion error
## which need to be fixed to completely run the program. 
# M <- readRDS("M4.RDS")
#M <- read.csv("PBEST188x282.csv")
#M <- as.matrix(M)
# write.csv(M, "PBEST188x282.csv", row.names = F)
# as.matrix(read.csv("PBEST188x282.csv"))
###########################################################################

################## 45x105 pooling matrix ##################################
# 45Times105kirkman.csv is same as "tapestry45x105.csv". 
## The original function uses 45Times105kirkman so I keeping it here.
## I am not replacing it by "tapestry45x105.csv" to avoid any kind of errors.
# This R program calls MATLAB and python; changing the file may cause conversion error
# which need to be fixed to completely run the program. 
M <- read.csv("45Times105kirkman.csv", header = T)
M <- as.matrix(M)
# write.csv(M, "tapestry45x105.csv", row.names = F)
# as.matrix(read.csv("tapestry45x105.csv"))
###########################################################################

####################### Pooling matrix for Dorfman testing ####################
# This is when 282 individuals are in the matrix
#M_Dorf <- M[1:48, ]

# This is when 105 individuals are in the matrix
M_Dorf <- read.csv("DorfMatrix_Spe105.csv") 
M_Dorf <- as.matrix(M_Dorf)
####################### Pooling matrix for optimal Dorfman testing ############
# The optimal group at p = 0.08 with Se = Sp = 1 or Se = 0.95, Sp = 1 is same
# So, depending on number of specimens, we either use "DorfMatrix_p0.08Spe282Se0.95Sp1.csv" or 
## "DorfMatrix_p0.08Spe105Se0.95Sp1.csv" for both _Se0.95 and _Se1 cases

# Pooling matrix for optimal Dorfman testing when testing 282 specimens
#M_Dorf_opt_Se1 <- read.csv("DorfMatrix_p0.08Spe282Se0.95Sp1.csv")

# Pooling matrix for optimal Dorfman testing when testing 105 specimens
M_Dorf_opt_Se1 <- read.csv("DorfMatrix_p0.08Spe105Se0.95Sp1.csv")
M_Dorf_opt_Se1 <- as.matrix(M_Dorf_opt_Se1)

# Pooling matrix for optimal Dorfman testing when testing 282 specimens
#M_Dorf_opt_Se0.95 <- read.csv("DorfMatrix_p0.08Spe282Se0.95Sp1.csv")
# Pooling matrix for optimal Dorfman testing when testing 282 specimens
M_Dorf_opt_Se0.95 <- read.csv("DorfMatrix_p0.08Spe105Se0.95Sp1.csv")
M_Dorf_opt_Se0.95 <- as.matrix(M_Dorf_opt_Se0.95)

source("FunctionsForApproximatedREDDILabDataAnalysis.R")
## ReddiLab filtered data
work_data <- read.csv("ReddiLabLognormalData.csv")
names(work_data)
#work_data2 <- work_data[, c(3, 17, 20)]
#names(work_data2)
#length(unique(work_data2$ID))

#work_data2$VL[work_data2$Ct >= 33.00] <- 0

#work_data2$ID <- 1:72016

#write.csv(work_data2, "ReddiLab_data_filterred2.csv", row.names = FALSE)


# prevalence
sum(work_data$VL >= 0.011)/length(work_data$VL) # prevalence is approx 0.08


#summary(work_data[work_data$Ct >= 32 & work_data$Ct < 33, 3])

# Total specimens
Total_spe <- dim(work_data)[1]
# Number of sample
Num_Sample <- dim(M)[2]

# Setting seed for reproducibility
set.seed(1234)
# ID's of specimens to test
Total_ID <- 1:Total_spe
active_ID <- sample(Total_ID, size = Num_Sample, replace = F)
# IDs that are tested
tested_ID <- c()
# Number of batches of specimens needed to test all the IDs
Num_batch <- ceiling(Total_spe/Num_Sample)

# Place to store all the results
store_results <- array(NA, dim = c(12, Num_batch, 19))
# Naming the dimensions of the array
dimnames(store_results) <- list(
  algorithm = c("Lasso_VL", "SBL_VL", "LS_VL", "NNOMP_VL", "PBest-2021_reduced",
                "PBest-2021_full", "PBest-2024",
                "Regular_testing_(retest_upon_ambiguity)",
                "Regular_testing_(always_retest)",
                "Dorfman_testing", "Optimal_Dorfman_testing_Se1",
                "Optimal_Dorfman_testing_Se0.95"),
  iteration = paste0("iter_", 1:Num_batch),
  quantity = c("Tested_neg_out_of_true", "Num_true_neg",
               "Tested_pos_out_of_true", "Num_true_pos",
               "False_pos", "False_neg", "True_neg_out_of_tested", "Num_test_neg",
               "True_pos_out_of_tested", "Num_test_pos", "Num_spec_need_model_prediction",
               "Num_spec_predicted_pos", "Num_spec_predicted_neg", "lambda",
               "Num_retest", "Total_test", "Model_Used", "Length_pred_pos_PBest2024",
               "Correct_pred_pos_PBest2024")
)


Rem <- 0

for(num_iter in 1:Num_batch){
  
  # Data of specimens to test in this batch
  active_work_data <- work_data$VL[active_ID]
  
  if(Rem == 0){
    
    # Results from all the algorithms
    res_AA <- res_All_Algo_ReddiLab(true_ind_VL = active_work_data,
                                    Num_Sample = 105, Num_Group = 45, Num_Group_Each_Sample = 3,
                                    M = M, M_Dorf = M_Dorf, M_Dorf_opt_Se1 = M_Dorf_opt_Se1,
                                    M_Dorf_opt_Se0.95 = M_Dorf_opt_Se0.95, sigma = 1,
                                    cut_off_VL = 0.011, Spe_Amnt_microL = 8)
    
  } else {
    # Number of samples in the remainder
    Num_Sample_rem <- length(active_ID)
    
    # M for the remainder samples
    M <- M[, 1:Num_Sample_rem]
    # Omitting the groups of M_rem that does not contain any specimens
    M <- M[which(apply(M, 1, sum) >= 1), ]
    # Number of groups in the M for remainder samples
    Num_Group_rem <- dim(M)[1]
    
    # M_Dorf for the remainder samples
    M_Dorf <- M_Dorf[, 1:Num_Sample_rem]
    # Omitting the groups of M_Dorf_rem that does not contain any specimens
    M_Dorf <- M_Dorf[which(apply(M_Dorf, 1, sum) >= 1), ]
    
    # M_Dorf_opt_Se1 for the remainder samples
    M_Dorf_opt_Se1 <- M_Dorf_opt_Se1[, 1:Num_Sample_rem]
    # Omitting the groups of M_Dorf_Se1_rem that does not contain any specimens
    M_Dorf_opt_Se1 <- M_Dorf_opt_Se1[which(apply(M_Dorf_opt_Se1, 1, sum) >= 1), ]
    
    # M_Dorf_opt_Se0.95 for the remainder samples
    M_Dorf_opt_Se0.95 <- M_Dorf_opt_Se0.95[, 1:Num_Sample_rem]
    # Omitting the groups of M_Dorf_Se0.95_rem that does not contain any specimens
    M_Dorf_opt_Se0.95 <- M_Dorf_opt_Se0.95[which(apply(M_Dorf_opt_Se0.95, 1, sum) >= 1), ]
    
    # Results from all the algorithms
    res_AA <- res_All_Algo_ReddiLab(true_ind_VL = active_work_data,
                                    Num_Sample = Num_Sample_rem, Num_Group = dim(M)[1],
                                    Num_Group_Each_Sample = 3,
                                    M = M, M_Dorf = M_Dorf, M_Dorf_opt_Se1 = M_Dorf_opt_Se1,
                                    M_Dorf_opt_Se0.95 = M_Dorf_opt_Se0.95, sigma = 1,
                                    cut_off_VL = 0.011, Spe_Amnt_microL = 8)
  }
  
  ######################## Lasso VL ###########################
  # Getting quantities to calculate accuracy of Lasso
  accuracy_Lasso_VL <- PSe(status = res_AA$true_status, result = res_AA$test_result_Lasso_VL)
  predictive_value_Lasso_VL <- PSe(status = res_AA$test_result_Lasso_VL, result = res_AA$true_status)
  
  # Storing results of Lasso when applied with VL
  store_results["Lasso_VL", num_iter, ] <- c(unlist(accuracy_Lasso_VL), unlist(predictive_value_Lasso_VL)[1:4],
                                             res_AA$num_specimens_needed_model_prediction,
                                             res_AA$num_specimens_Lasso_predicted_positive_VL,
                                             res_AA$num_specimens_Lasso_predicted_negative_VL,
                                             res_AA$lmd_Lasso_VL, res_AA$num_retest,
                                             dim(M)[1] + res_AA$num_retest,
                                             res_AA$Model_Used, NA, NA)
  
  ######################## SBL VL ################################
  # Getting quantities to calculate accuracy of SBL
  accuracy_SBL_VL <- PSe(status = res_AA$true_status, result = res_AA$test_result_SBL_VL)
  predictive_value_SBL_VL <- PSe(status = res_AA$test_result_SBL_VL, result = res_AA$true_status)
  
  # Storing results of SBL when applied with VL
  store_results["SBL_VL", num_iter, ] <- c(unlist(accuracy_SBL_VL), unlist(predictive_value_SBL_VL)[1:4], 
                                           res_AA$num_specimens_needed_model_prediction,
                                           res_AA$num_specimens_SBL_predicted_positive_VL,
                                           res_AA$num_specimens_SBL_predicted_negative_VL, NA,
                                           res_AA$num_retest, dim(M)[1] + res_AA$num_retest,
                                           res_AA$Model_Used, NA, NA)
  
  
  
  ##################### LS VL ##################################
  # Getting quantities to calculate accuracy of LS
  accuracy_LS_VL <- PSe(status = res_AA$true_status, result = res_AA$test_result_LS_VL)
  predictive_value_LS_VL <- PSe(status = res_AA$test_result_LS_VL, result = res_AA$true_status)
  
  # Storing results of LS when applied with VL
  store_results["LS_VL", num_iter, ] <- c(unlist(accuracy_LS_VL), unlist(predictive_value_LS_VL)[1:4], 
                                          res_AA$num_specimens_needed_model_prediction,
                                          res_AA$num_specimens_LS_predicted_positive_VL,
                                          res_AA$num_specimens_LS_predicted_negative_VL, NA,
                                          res_AA$num_retest, dim(M)[1] + res_AA$num_retest,
                                          res_AA$Model_Used, NA, NA)
  
  #################### NNOMP VL ################################
  # Getting quantities to calculate accuracy of NNOMP
  accuracy_NNOMP_VL <- PSe(status = res_AA$true_status, result = res_AA$test_result_NNOMP_VL)
  predictive_value_NNOMP_VL <- PSe(status = res_AA$test_result_NNOMP_VL, result = res_AA$true_status)
  
  # Storing results of NNOMP when applied with VL
  store_results["NNOMP_VL", num_iter, ] <- c(unlist(accuracy_NNOMP_VL), unlist(predictive_value_NNOMP_VL)[1:4],
                                             res_AA$num_specimens_needed_model_prediction,
                                             res_AA$num_specimens_NNOMP_predicted_positive_VL,
                                             res_AA$num_specimens_NNOMP_predicted_negative_VL,
                                             NA, res_AA$num_retest, dim(M)[1] + res_AA$num_retest,
                                             res_AA$Model_Used, NA, NA)
  
  
  ####################### P-Best reduced ##################
  # Getting quantities to calculate accuracy of P-Best with reduced pooling matrix
  accuracy_PBest_reduced <- PSe(status = res_AA$true_status, result = res_AA$test_result_PBest_reduced)
  predictive_value_PBest_reduced <- PSe(status = res_AA$test_result_PBest_reduced, result = res_AA$true_status)
  
  # Storing results of P-Best when applied with VL and reduced pooling matrix
  store_results["PBest-2021_reduced", num_iter, ] <- c(unlist(accuracy_PBest_reduced), unlist(predictive_value_PBest_reduced)[1:4],
                                                       res_AA$num_specimens_needed_model_prediction, 
                                                       res_AA$num_specimens_PBest_reduced_predicted_positive,
                                                       res_AA$num_specimens_PBest_reduced_predicted_negative,
                                                       res_AA$PBest_reduced_lambda,
                                                       res_AA$num_retest, dim(M)[1] + res_AA$num_retest,
                                                       res_AA$Model_Used, NA, NA)
  
  ####################### P-Best full ####################
  # Getting quantities to calculate accuracy of P-Best with full pooling matrix
  accuracy_PBest_full <- PSe(status = res_AA$true_status, result = res_AA$test_result_PBest_full)
  predicte_value_PBest_full <- PSe(status = res_AA$test_result_PBest_full, result = res_AA$true_status)
  
  # Storing results of P-Best when applied with full pooling matrix
  store_results["PBest-2021_full", num_iter, ] <- c(unlist(accuracy_PBest_full),
                                                    unlist(predicte_value_PBest_full)[1:4],
                                                    res_AA$num_specimens_needed_model_prediction,
                                                    res_AA$num_specimens_PBest_full_predicted_positive,
                                                    res_AA$num_specimens_PBest_full_predicted_negative,
                                                    res_AA$PBest_full_lambda, res_AA$num_retest,
                                                    dim(M)[1] + res_AA$num_retest,
                                                    NA, NA, NA)
  
  
  ####################### P-Best 2024 ####################
  # Getting quantities to calculate accuracy of P-Best 2024
  accuracy_PBest_2024 <- PSe(status = res_AA$true_status, result = res_AA$res_PBest_2024$test_result)
  predicte_value_PBest_2024 <- PSe(status = res_AA$res_PBest_2024$test_result, result = res_AA$true_status)
  
  # Storing results of P-Best when applied with full pooling matrix
  store_results["PBest-2024", num_iter, ] <- c(unlist(accuracy_PBest_2024),
                                               unlist(predicte_value_PBest_2024)[1:4],
                                               NA,
                                               NA,
                                               NA,
                                               res_AA$res_PBest_2024$lambda,
                                               res_AA$res_PBest_2024$Num_retest,
                                               dim(M)[1] + res_AA$res_PBest_2024$Num_retest,
                                               NA, res_AA$res_PBest_2024$len_pred_pos,
                                               res_AA$res_PBest_2024$correct_pos)
  
  
  
  ####################### Regular testing (retest upon ambiguity) ###################
  # Getting quantities to calculate accuracy of regular testing (retest upon ambiguity)
  accuracy_regular_retest_upon_ambiguity <- PSe(status = res_AA$true_status,
                                                result = res_AA$test_result_regular)
  predictive_value_regular_retest_upon_ambiguity <- PSe(status = res_AA$test_result_regular,
                                                        result = res_AA$true_status)
  
  # Storing results of regular testing (retest upon ambiguity) 
  store_results["Regular_testing_(retest_upon_ambiguity)", num_iter, ] <- c(unlist(accuracy_regular_retest_upon_ambiguity),
                                                                            unlist(predictive_value_regular_retest_upon_ambiguity)[1:4],
                                                                            res_AA$num_specimens_needed_model_prediction,
                                                                            res_AA$num_specimens_tested_positive_regular,
                                                                            res_AA$num_specimens_tested_negative_regular,
                                                                            NA, res_AA$num_retest_regular,
                                                                            dim(M)[1] + res_AA$num_retest_regular,
                                                                            res_AA$Model_Used,
                                                                            NA, NA)
  
  ####################### Regular testing (always retest) ########################
  # Getting quantities to calculate accuracy of regular testing (always retest)
  accuracy_regular_always_retest <- PSe(status = res_AA$true_status,
                                        result = res_AA$test_result_regular_AR)
  predictive_value_regular_always_retest <- PSe(status = res_AA$test_result_regular_AR,
                                                result = res_AA$true_status)
  # Storing results of regular testing (always retest)
  store_results["Regular_testing_(always_retest)", num_iter, ] <- c(unlist(accuracy_regular_always_retest),
                                                                    unlist(predictive_value_regular_always_retest)[1:4],
                                                                    NA, NA, NA, NA, res_AA$num_retest_regular_AR,
                                                                    dim(M)[1] + res_AA$num_retest_regular_AR, NA, NA, NA)
  
  
  ###################### Dorfman testing #########################
  # Getting quantities to calculate accuracy of Dorfman testing
  accuracy_Dorfman <- PSe(status = res_AA$true_status, result = res_AA$res_Dorfman$test_result)
  predictive_value_Dorfman <- PSe(status = res_AA$res_Dorfman$test_result,
                                  result = res_AA$true_status)
  # Storing results of Dorfman testing
  store_results["Dorfman_testing", num_iter, ] <- c(unlist(accuracy_Dorfman),
                                                    unlist(predictive_value_Dorfman)[1:4],
                                                    NA, NA, NA, NA,
                                                    res_AA$res_Dorfman$Num_retest,
                                                    dim(M_Dorf)[1] + res_AA$res_Dorfman$Num_retest, NA,
                                                    NA, NA)
  
  
  ###################### Optimal Dorfman testing at Se = Sp = 1#########################
  # Getting quantities to calculate accuracy of optimal Dorfman testing at Se = Sp = 1
  accuracy_Dorfman_opt_Se1 <- PSe(status = res_AA$true_status, result = res_AA$res_Dorfman_opt_Se1$test_result)
  predictive_value_Dorfman_opt_Se1 <- PSe(status = res_AA$res_Dorfman_opt_Se1$test_result,
                                          result = res_AA$true_status)
  
  # Storing results of optimal Dorfman testing at Se = Sp = 1
  store_results["Optimal_Dorfman_testing_Se1", num_iter, ] <- c(unlist(accuracy_Dorfman_opt_Se1),
                                                                unlist(predictive_value_Dorfman_opt_Se1)[1:4],
                                                                NA, NA, NA, NA,
                                                                res_AA$res_Dorfman_opt_Se1$Num_retest,
                                                                dim(M_Dorf_opt_Se1)[1] + res_AA$res_Dorfman_opt_Se1$Num_retest, NA,
                                                                NA, NA)
  
  ###################### Optimal Dorfman testing at Se = 0.95, Sp = 1#########################
  # Getting quantities to calculate accuracy of optimal Dorfman testing at Se = 0.95, Sp = 1
  accuracy_Dorfman_opt_Se0.95 <- PSe(status = res_AA$true_status, result = res_AA$res_Dorfman_opt_Se0.95$test_result)
  predictive_value_Dorfman_opt_Se0.95 <- PSe(status = res_AA$res_Dorfman_opt_Se0.95$test_result,
                                             result = res_AA$true_status)
  
  # Storing results of optimal Dorfman testing at Se = Sp = 1
  store_results["Optimal_Dorfman_testing_Se0.95", num_iter, ] <- c(unlist(accuracy_Dorfman_opt_Se0.95),
                                                                   unlist(predictive_value_Dorfman_opt_Se0.95)[1:4],
                                                                   NA, NA, NA, NA,
                                                                   res_AA$res_Dorfman_opt_Se0.95$Num_retest,
                                                                   dim(M_Dorf_opt_Se0.95)[1] + res_AA$res_Dorfman_opt_Se0.95$Num_retest,
                                                                   NA, NA, NA)
  # IDs that are tested
  tested_ID <- c(tested_ID, active_ID)
  # IDs that are available to test
  candidate_ID <- setdiff(Total_ID, tested_ID)
  
  if(length(candidate_ID) >= Num_Sample){
    # IDs to test in the next batch
    active_ID <- sample(candidate_ID, size = Num_Sample, replace = F)
    # Indicator of remainder group
    Rem <- 0
  } else{
    # IDs to test will be same as candidate IDs
    active_ID <- candidate_ID
    # Indicator of testing the remainder group
    Rem <- 1
  }
  
  cat("Num batch = ", num_iter, "\n")
}

res_across_all_iter <- apply(store_results[, , c(1:4, 7:10)], c(1, 3), sum, na.rm = TRUE)
# Calculating PPA
PPA <- res_across_all_iter[ , "Tested_pos_out_of_true"]/res_across_all_iter[, "Num_true_pos"]
# Calculating NPA
NPA <- res_across_all_iter[ , "Tested_neg_out_of_true"]/res_across_all_iter[ , "Num_true_neg"]
# Calculating PPV
PPV <- res_across_all_iter[ , "True_pos_out_of_tested"]/res_across_all_iter[ , "Num_test_pos"]
# Calculating NPV
NPV <- res_across_all_iter[ , "True_neg_out_of_tested"]/res_across_all_iter[ , "Num_test_neg"]
# Expected number of retest
mean_retest <- apply(store_results[, 1:(Num_batch - 1), "Num_retest"], 1, mean)
# Expected number of tests
ET <- apply(store_results[, 1:Num_batch, "Total_test"], 1, mean)


# Expected number of tests per individual
ET.per <- ET/Num_Sample
# Calculating proportion of ambiguity
Ambi_Prop <- apply(store_results[, , "Model_Used"], 1, mean, na.rm = TRUE)

# Assigning ambiguity proportion equal NA for P-Best with full pooling matrix,
## regular testing (always retest), Dorfman testing as it does not make sense
## to calculate ambiguity proportion as they always use the full pooling matrix
Ambi_Prop[c(6, 7, 9, 10, 11, 12)] <- NA

# Results summary
res_summary <- as.data.frame(cbind(Ambiguity_Proportion = Ambi_Prop, ET = ET, ET.per, PPA, NPA, PPV, NPV))
# Taking rownames to a column name
res_summary <- rownames_to_column(res_summary, var = "Algorithm")

################### Few names change in the paper ###########
# SBL is spare Bayesian Learning
# LS is least-squares
# P-Best_2021_reduced is now "P-BEST 2020*"
# P-Best_2021_full is now "P-BEST 2020"
# Regular testing is now "standard retesting"
# Dorfman testing is now "regular Dorfman"
# ET.per is now E
# PPA is now PSe
# NPA is now PSp
# PPV is now PPPV
# NPV is now PNPV
#############################################################

# Results summary
res_summary
# Writing the summary file
write.csv(res_summary, "DistBasedReddiLab_DataAnalysis_results188Times282PoolingMatrix8microLSigma1.csv",
          row.names = F)

# Making the file 2-dimensional
flat_raw_results <- as.data.frame.table(store_results, responseName = "value")
write.csv(flat_raw_results, file = "DistBasedReddiLabRawResultsOn188Times282PoolingMatrix8microLSigma1.csv", row.names = F)

