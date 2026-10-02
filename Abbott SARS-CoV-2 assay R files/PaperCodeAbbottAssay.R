# Cleaning the R environment
rm(list = ls())
#### All the packages needed ####
library(glmnet)
library(reticulate)
library(R.matlab)
library(tibble)
################################
# Data used in simulation from Hirschhorn et al. (2021)
ct <- read.csv("HirchhornProb.csv")
# CT values
CT <- ct$CT
# Probability for each CT value
CT_prob <- ct$Prob

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
# use_python("/home/bilder/pranta/.conda/envs/r-matlab311/bin/python", required = TRUE)
################################################################################

############################## Importing Tapestry functions ####################
sbl <- import("sbl")
l1ls <- import("l1ls")
nnompcv <-  import('nnompcv')

####################### Reading pooling matrix #########################
#M <- readMat("PBestMatlabCode/poolingMatrix.mat")
#M <- M$poolingMatrix

# The "tapestry93x961.csv" is created by taking first 961 columns of the matrix obtained from "optimized_M_93_1240_kirkman.txt" as done by Ghosh et al. (2021)
# I am directly reading "tapestry93x961.csv" to avoid any kind of errors. This R program calls MATLAB and python; changing the file may cause conversion error
# which need to be fixed to completely run the program.
M <- read.table("optimized_M_93_1240_kirkman.txt", header = F)
M <- M[, 1:961]
M <- matrix(as.numeric(as.matrix(M)), nrow = nrow(M))

 
# write.csv(M, "tapestry93x961.csv", row.names = F)
# apply(read.csv("tapestry93x961.csv"), 1, sum)
# apply(read.csv("tapestry93x961.csv"), 2, sum)
####################### Pooling matrix for Dorfman testing ####################
#M_Dorf <- M[1:8, ] # Group size is 48 with PBest's 48 by 384 matrix

# Group size is 31 with Tapestry's 93 by 961 matrix. So, we create a 31x961 matrix where each group has a size of 31 and 
## the specimens are non-overlaping.
M_Dorf <- read.csv("Dorf_Pool_matrix.csv") 

####################### Pooling matrix for optimal Dorfman testing ############
# Depending on p, we change the pooling matrix used for optimal Dorfman testing. In our most recent run, we used p = 0.02.
# Note that at p = 0.005, 0.01, 0.015 and 0.02 with Se = 0.95, Sp = 1 or Se = Sp = 1, the optimal group size does not vary.
## Accordingly, we use the file with extension "_Se0.95Sp1.csv" for every p except p = 0.015.

# Matrix for optimal Dorfman testing at Se = Sp = 1 
# The optimal group size at p = 0.005 at Se = Sp = 1 is 15
# This group size is same as at p = 5/961 with Se = 0.95, Sp = 1. So, we are using the file with "_p5by961_" in the name to keep the number of files limited.
#M_Dorf_opt_Se1 <- read.csv("DorfMatrix_p5by961_Se0.95Sp1.csv")

# The optimal group size at p = 0.01 at Se = Sp = 1 is 11
#M_Dorf_opt_Se1 <- read.csv("DorfMatrix_p10by961_Se0.95Sp1.csv")
# The optimal group size at p = 15/961 at Se = Sp = 1 is 9
# The optimal group size at p = 0.015 at Se = Sp = 1 is 9
#M_Dorf_opt_Se1 <- read.csv("DorfMatrix_p15by961_SeSp1.csv")
# The optimal group size at p = 0.02 at Se = Sp = 1 is 8
M_Dorf_opt_Se1 <- read.csv("DorfMatrix_p20by961_Se0.95Sp1.csv")

# Matrix for optimal Dorfman testing at Se = 0.95, Sp = 1
# The optimal group size at p = 5/961 at Se = 0.95 and Sp = 1 is 15
# The optimal group size at p = 0.005 at Se = 0.95, Sp = 1 is 15
# M_Dorf_opt_Se0.95 <- read.csv("DorfMatrix_p5by961_Se0.95Sp1.csv")
# The optimal group size at p = 10/961 at Se = 0.95 and Sp = 1 is 11
# The optimal group size at p = 0.01 at Se = 0.95 and Sp = 1 is 11
#M_Dorf_opt_Se0.95 <- read.csv("DorfMatrix_p10by961_Se0.95Sp1.csv")
# The optimal group size at p = 15/961 at Se = 0.95 and Sp = 1 is 9
# The optimal group size at p = 0.015 at Se = 0.95 and Sp = 1 is 9
#M_Dorf_opt_Se0.95 <- read.csv("DorfMatrix_p15by961_SeSp1.csv")
# The optimal group size at p = 20/961 at Se = 0.95 and Sp = 1 is 8
# The optimal group size at p = 0.02 at Se = 0.95 and Sp = 1 is 8
M_Dorf_opt_Se0.95 <- read.csv("DorfMatrix_p20by961_Se0.95Sp1.csv")


# Function to get summary of the results of each algorithm
get_res_summary <- function(B, p, Num_Sample, Num_Group, Num_Group_Each_Sample, 
         M, M_Dorf, M_Dorf_opt_Se1, M_Dorf_opt_Se0.95,
         CT, CT_prob, sigma = 1,
         cut_off_CT = 31.5,
         Spe_Amnt_mL = 0.8,
         LoD = 26.73){

  # Place to store all the results
  store_res_All_Algo <- array(NA, dim = c(16, B, 18))
  # Naming the dimensions of the array
  dimnames(store_res_All_Algo) <- list(
    algorithm = c("Lasso_VL", "SBL_VL", "LS_VL", "NNOMP_VL", "PBest-2021_reduced", "PBest-2021_full", "PBest-2024",
                  "Regular_testing_(retest_upon_ambiguity)", "Regular_testing_(always_retest)",
                  "Dorfman_testing", "Optimal_Dorfman_testing_Se1", "Optimal_Dorfman_testing_Se0.95",
                  "Lasso_CT", "SBL_CT", "LS_CT", "NNOMP_CT"),
    iteration = paste0("iter_", 1:B),
    quantity = c("Tested_neg_out_of_true", "Num_true_neg",
      "Tested_pos_out_of_true", "Num_true_pos",
      "False_pos", "False_neg", "True_neg_out_of_tested", "Num_test_neg",
      "True_pos_out_of_tested", "Num_test_pos", "Num_spec_need_model_prediction",
      "Num_spec_predicted_pos", "Num_spec_predicted_neg", "lambda",
      "Num_retest", "Model_Used", "Length_pred_pos_PBest2024",
      "Correct_pred_pos_PBest2024")
  )

  for(num_iter in 1:B){  
  
    # Results from all the algorithms
    res_AA <- res_All_Algo(p = p, Num_Sample = Num_Sample, Num_Group = Num_Group,
                          Num_Group_Each_Sample = Num_Group_Each_Sample, M = M,
                          M_Dorf = M_Dorf, M_Dorf_opt_Se1 = M_Dorf_opt_Se1,
                          M_Dorf_opt_Se0.95 = M_Dorf_opt_Se0.95,
                          CT = CT, CT_prob = CT_prob,
                          sigma = sigma, cut_off_CT = cut_off_CT,
                          Spe_Amnt_mL = Spe_Amnt_mL, LoD = LoD)
  
    # Using PSe() to get the quantities required to calculate PPA, NPA, PPV, NPV
    # Such as number of tested positives, true positives, false positives, false negatives
    ######################## Lasso VL ###########################
    # Getting quantities to calculate accuracy of Lasso
    accuracy_Lasso_VL <- PSe(status = res_AA$true_status, result = res_AA$test_result_Lasso_VL)
    predictive_value_Lasso_VL <- PSe(status = res_AA$test_result_Lasso_VL, result = res_AA$true_status)
  
    # Storing results of Lasso when applied with VL
    store_res_All_Algo["Lasso_VL", num_iter, ] <- c(unlist(accuracy_Lasso_VL), unlist(predictive_value_Lasso_VL)[1:4],
                                                    res_AA$num_specimens_needed_model_prediction,
                                                    res_AA$num_specimens_Lasso_predicted_positive_VL,
                                                    res_AA$num_specimens_Lasso_predicted_negative_VL,
                                                    res_AA$lmd_Lasso_VL, res_AA$num_retest,
                                                    res_AA$Model_Used, NA, NA)
  
  ######################## SBL VL ################################
  # Getting quantities to calculate accuracy of SBL
  accuracy_SBL_VL <- PSe(status = res_AA$true_status, result = res_AA$test_result_SBL_VL)
  predictive_value_SBL_VL <- PSe(status = res_AA$test_result_SBL_VL, result = res_AA$true_status)
  
  # Storing results of SBL when applied with VL
  store_res_All_Algo["SBL_VL", num_iter, ] <- c(unlist(accuracy_SBL_VL), unlist(predictive_value_SBL_VL)[1:4], 
                                                res_AA$num_specimens_needed_model_prediction,
                                                res_AA$num_specimens_SBL_predicted_positive_VL,
                                                res_AA$num_specimens_SBL_predicted_negative_VL, NA,
                                                res_AA$num_retest, res_AA$Model_Used, NA, NA)
  
  
  
  ##################### LS VL ##################################
  # Getting quantities to calculate accuracy of LS
  accuracy_LS_VL <- PSe(status = res_AA$true_status, result = res_AA$test_result_LS_VL)
  predictive_value_LS_VL <- PSe(status = res_AA$test_result_LS_VL, result = res_AA$true_status)
  
  # Storing results of LS when applied with VL
  store_res_All_Algo["LS_VL", num_iter, ] <- c(unlist(accuracy_LS_VL), unlist(predictive_value_LS_VL)[1:4], 
                                               res_AA$num_specimens_needed_model_prediction,
                                               res_AA$num_specimens_LS_predicted_positive_VL,
                                               res_AA$num_specimens_LS_predicted_negative_VL, NA,
                                               res_AA$num_retest, res_AA$Model_Used, NA, NA)
  
  #################### NNOMP VL ################################
  # Getting quantities to calculate accuracy of NNOMP
  accuracy_NNOMP_VL <- PSe(status = res_AA$true_status, result = res_AA$test_result_NNOMP_VL)
  predictive_value_NNOMP_VL <- PSe(status = res_AA$test_result_NNOMP_VL, result = res_AA$true_status)
  
  # Storing results of NNOMP when applied with VL
  store_res_All_Algo["NNOMP_VL", num_iter, ] <- c(unlist(accuracy_NNOMP_VL), unlist(predictive_value_NNOMP_VL)[1:4],
                                                  res_AA$num_specimens_needed_model_prediction,
                                                  res_AA$num_specimens_NNOMP_predicted_positive_VL,
                                                  res_AA$num_specimens_NNOMP_predicted_negative_VL,
                                                  NA, res_AA$num_retest, res_AA$Model_Used, NA, NA)
  
  
  ####################### P-Best reduced ##################
  # Getting quantities to calculate accuracy of P-Best with reduced pooling matrix
  accuracy_PBest_reduced <- PSe(status = res_AA$true_status, result = res_AA$test_result_PBest_reduced)
  predictive_value_PBest_reduced <- PSe(status = res_AA$test_result_PBest_reduced, result = res_AA$true_status)
  
  # Storing results of P-Best when applied with VL and reduced pooling matrix
  store_res_All_Algo["PBest-2021_reduced", num_iter, ] <- c(unlist(accuracy_PBest_reduced), unlist(predictive_value_PBest_reduced)[1:4],
                                                       res_AA$num_specimens_needed_model_prediction, 
                                                       res_AA$num_specimens_PBest_reduced_predicted_positive,
                                                       res_AA$num_specimens_PBest_reduced_predicted_negative,
                                                       res_AA$PBest_reduced_lambda,
                                                       res_AA$num_retest, res_AA$Model_Used, NA, NA)
  
  ####################### P-Best full ####################
  # Getting quantities to calculate accuracy of P-Best with full pooling matrix
  accuracy_PBest_full <- PSe(status = res_AA$true_status, result = res_AA$test_result_PBest_full)
  predicte_value_PBest_full <- PSe(status = res_AA$test_result_PBest_full, result = res_AA$true_status)
  
  # Storing results of P-Best when applied with full pooling matrix
  store_res_All_Algo["PBest-2021_full", num_iter, ] <- c(unlist(accuracy_PBest_full),
                                                    unlist(predicte_value_PBest_full)[1:4],
                                                    res_AA$num_specimens_needed_model_prediction,
                                                    res_AA$num_specimens_PBest_full_predicted_positive,
                                                    res_AA$num_specimens_PBest_full_predicted_negative,
                                                    res_AA$PBest_full_lambda, res_AA$num_retest, NA, NA, NA)
  
  
  ####################### P-Best 2024 ####################
  # Getting quantities to calculate accuracy of P-Best 2024
  accuracy_PBest_2024 <- PSe(status = res_AA$true_status, result = res_AA$res_PBest_2024$test_result)
  predicte_value_PBest_2024 <- PSe(status = res_AA$res_PBest_2024$test_result, result = res_AA$true_status)
  
  # Storing results of P-Best when applied with full pooling matrix
  store_res_All_Algo["PBest-2024", num_iter, ] <- c(unlist(accuracy_PBest_2024),
                                                    unlist(predicte_value_PBest_2024)[1:4],
                                                    NA,
                                                    NA,
                                                    NA,
                                                    res_AA$res_PBest_2024$lambda,
                                                    res_AA$res_PBest_2024$Num_retest,
                                                    NA, res_AA$res_PBest_2024$len_pred_pos,
                                                    res_AA$res_PBest_2024$correct_pos)
  
  
  
  
  ####################### Regular testing (retest upon ambiguity) ###################
  # Getting quantities to calculate accuracy of regular testing (retest upon ambiguity)
  accuracy_regular_retest_upon_ambiguity <- PSe(status = res_AA$true_status,
                                                result = res_AA$test_result_regular)
  predictive_value_regular_retest_upon_ambiguity <- PSe(status = res_AA$test_result_regular,
                                                        result = res_AA$true_status)
  
  # Storing results of regular testing (retest upon ambiguity) 
  store_res_All_Algo["Regular_testing_(retest_upon_ambiguity)", num_iter, ] <- c(unlist(accuracy_regular_retest_upon_ambiguity),
                                                                                 unlist(predictive_value_regular_retest_upon_ambiguity)[1:4],
                                                                                 res_AA$num_specimens_needed_model_prediction,
                                                                                 res_AA$num_specimens_tested_positive_regular,
                                                                                 res_AA$num_specimens_tested_negative_regular,
                                                                                 NA, res_AA$num_retest_regular, res_AA$Model_Used,
                                                                                 NA, NA)
  
  ####################### Regular testing (always retest) ########################
  # Getting quantities to calculate accuracy of regular testing (always retest)
  accuracy_regular_always_retest <- PSe(status = res_AA$true_status,
                                        result = res_AA$test_result_regular_AR)
  predictive_value_regular_always_retest <- PSe(status = res_AA$test_result_regular_AR,
                                              result = res_AA$true_status)
  # Storing results of regular testing (always retest)
  store_res_All_Algo["Regular_testing_(always_retest)", num_iter, ] <- c(unlist(accuracy_regular_always_retest),
                                                                         unlist(predictive_value_regular_always_retest)[1:4],
                                                                         NA, NA, NA, NA, res_AA$num_retest_regular_AR, NA, NA, NA)
  
  
  ###################### Dorfman testing #########################
  # Getting quantities to calculate accuracy of Dorfman testing
  accuracy_Dorfman <- PSe(status = res_AA$true_status, result = res_AA$res_Dorfman$test_result)
  predictive_value_Dorfman <- PSe(status = res_AA$res_Dorfman$test_result,
                                  result = res_AA$true_status)
  # Storing results of Dorfman testing
  store_res_All_Algo["Dorfman_testing", num_iter, ] <- c(unlist(accuracy_Dorfman),
                                                         unlist(predictive_value_Dorfman)[1:4],
                                                         NA, NA, NA, NA,
                                                         res_AA$res_Dorfman$Num_retest, NA,
                                                         NA, NA)
  
  
  ###################### Optimal Dorfman testing at Se = Sp = 1#########################
  # Getting quantities to calculate accuracy of optimal Dorfman testing at Se = Sp = 1
  accuracy_Dorfman_opt_Se1 <- PSe(status = res_AA$true_status, result = res_AA$res_Dorfman_opt_Se1$test_result)
  predictive_value_Dorfman_opt_Se1 <- PSe(status = res_AA$res_Dorfman_opt_Se1$test_result,
                                  result = res_AA$true_status)

  # Storing results of optimal Dorfman testing at Se = Sp = 1
  store_res_All_Algo["Optimal_Dorfman_testing_Se1", num_iter, ] <- c(unlist(accuracy_Dorfman_opt_Se1),
                                                         unlist(predictive_value_Dorfman_opt_Se1)[1:4],
                                                         NA, NA, NA, NA,
                                                         res_AA$res_Dorfman_opt_Se1$Num_retest, NA,
                                                         NA, NA)
  
  ###################### Optimal Dorfman testing at Se = 0.95, Sp = 1#########################
  # Getting quantities to calculate accuracy of optimal Dorfman testing at Se = 0.95, Sp = 1
  accuracy_Dorfman_opt_Se0.95 <- PSe(status = res_AA$true_status, result = res_AA$res_Dorfman_opt_Se0.95$test_result)
  predictive_value_Dorfman_opt_Se0.95 <- PSe(status = res_AA$res_Dorfman_opt_Se0.95$test_result,
                                          result = res_AA$true_status)
  
  # Storing results of optimal Dorfman testing at Se = Sp = 1
  store_res_All_Algo["Optimal_Dorfman_testing_Se0.95", num_iter, ] <- c(unlist(accuracy_Dorfman_opt_Se0.95),
                                                                     unlist(predictive_value_Dorfman_opt_Se0.95)[1:4],
                                                                     NA, NA, NA, NA,
                                                                     res_AA$res_Dorfman_opt_Se0.95$Num_retest, NA,
                                                                     NA, NA)
  
  
  ##################### Lasso with scaled CT #####################
  # Getting quantities to calculate accuracy of Lasso with scaled CT
  accuracy_Lasso_CT <- PSe(status = res_AA$true_status, result = res_AA$test_result_Lasso_CT)
  predictive_value_Lasso_CT <- PSe(status = res_AA$test_result_Lasso_CT,
                                   result = res_AA$true_status)
  
  # Storing results of Lasso when applied with scaled CT
  store_res_All_Algo["Lasso_CT", num_iter, ] <- c(unlist(accuracy_Lasso_CT), unlist(predictive_value_Lasso_CT)[1:4],
                                                  res_AA$num_specimens_needed_model_prediction,
                                                  res_AA$num_specimens_Lasso_predicted_positive_CT,
                                                  res_AA$num_specimens_Lasso_predicted_negative_CT,
                                                  res_AA$lmd_Lasso_CT, res_AA$num_retest, res_AA$Model_Used, NA, NA)
  
  ##################### SBL with scaled CT #####################
  # Getting quantities to calculate accuracy of SBL with scaled CT
  accuracy_SBL_CT <- PSe(status = res_AA$true_status, result = res_AA$test_result_SBL_CT)
  predictive_value_SBL_CT <- PSe(status = res_AA$test_result_SBL_CT,
                                   result = res_AA$true_status)
  
  # Storing results of SBL when applied with scaled CT
  store_res_All_Algo["SBL_CT", num_iter, ] <- c(unlist(accuracy_SBL_CT), unlist(predictive_value_SBL_CT)[1:4],
                                                  res_AA$num_specimens_needed_model_prediction,
                                                  res_AA$num_specimens_SBL_predicted_positive_CT,
                                                  res_AA$num_specimens_SBL_predicted_negative_CT,
                                                  NA, res_AA$num_retest, res_AA$Model_Used, NA, NA)
  
  ##################### Least squares (LS) with scaled CT #####################
  # Getting quantities to calculate accuracy of LS with scaled CT
  accuracy_LS_CT <- PSe(status = res_AA$true_status, result = res_AA$test_result_LS_CT)
  predictive_value_LS_CT <- PSe(status = res_AA$test_result_LS_CT,
                                result = res_AA$true_status)
  
  # Storing results of LS when applied with scaled CT
  store_res_All_Algo["LS_CT", num_iter, ] <- c(unlist(accuracy_LS_CT), unlist(predictive_value_LS_CT)[1:4],
                                                res_AA$num_specimens_needed_model_prediction,
                                                res_AA$num_specimens_LS_predicted_positive_CT,
                                                res_AA$num_specimens_LS_predicted_negative_CT,
                                                0.10, res_AA$num_retest, res_AA$Model_Used, NA, NA)
  
  ##################### NNOMP with scaled CT #####################
  # Getting quantities to calculate accuracy of NNOMP with scaled CT
  accuracy_NNOMP_CT <- PSe(status = res_AA$true_status, result = res_AA$test_result_NNOMP_CT)
  predictive_value_NNOMP_CT <- PSe(status = res_AA$test_result_NNOMP_CT,
                                   result = res_AA$true_status)
  
  # Storing results of NNOMP when applied with scaled CT
  store_res_All_Algo["NNOMP_CT", num_iter, ] <- c(unlist(accuracy_NNOMP_CT), unlist(predictive_value_NNOMP_CT)[1:4],
                                               res_AA$num_specimens_needed_model_prediction,
                                               res_AA$num_specimens_NNOMP_predicted_positive_CT,
                                               res_AA$num_specimens_NNOMP_predicted_negative_CT,
                                               NA, res_AA$num_retest, res_AA$Model_Used, NA, NA)
  cat("p = ", p, "Iteration = ", num_iter, "\n")

}

res_across_all_iter <- apply(store_res_All_Algo[, , c(1:4, 7:10)], c(1, 3), sum, na.rm = TRUE)
# Calculating PPA
PPA <- res_across_all_iter[ , "Tested_pos_out_of_true"]/res_across_all_iter[, "Num_true_pos"]
# Calculating NPA
NPA <- res_across_all_iter[ , "Tested_neg_out_of_true"]/res_across_all_iter[ , "Num_true_neg"]
# Calculating PPV
PPV <- res_across_all_iter[ , "True_pos_out_of_tested"]/res_across_all_iter[ , "Num_test_pos"]
# Calculating NPV
NPV <- res_across_all_iter[ , "True_neg_out_of_tested"]/res_across_all_iter[ , "Num_test_neg"]
# Expected number of retest
mean_retest <- apply(store_res_All_Algo[, , "Num_retest"], 1, mean, na.rm = TRUE)
# Expected number of tests
ET <- c(rep(Num_Group, 9), dim(M_Dorf)[1], dim(M_Dorf_opt_Se1)[1], dim(M_Dorf_opt_Se0.95)[1],
        rep(Num_Group, 4)) + mean_retest

# Expected number of tests per individual
ET.per <- ET/Num_Sample
# Calculating proportion of ambiguity
Ambi_Prop <- apply(store_res_All_Algo[, , "Model_Used"], 1, mean, na.rm = TRUE)

# Assigning ambiguity proportion equal NA for P-Best with full pooling matrix,
## regular testing (always retest), Dorfman testing as it does not make sense
## to calculate ambiguity proportion as they always use the full pooling matrix
Ambi_Prop[c(6, 7, 9, 10, 11, 12)] <- NA

# Results summary
res_summary <- as.data.frame(cbind(Ambiguity_Proportion = Ambi_Prop, ET, ET.per, PPA, NPA, PPV, NPV))
# Taking rownames to a column name
res_summary <- rownames_to_column(res_summary, var = "Algorithm")
# Listing output
list(res_summary = res_summary, raw_results = store_res_All_Algo)
}

# Loading all other necessary functions
source("NecessaryFunctions.R")
# Setting seed for reproducibility
set.seed(1234)
start_time <- proc.time()
# Running get_res_summary() for get the measures needed
My_results <- get_res_summary(B = 10000, p = 0.02,
                              Num_Sample = 961, Num_Group = 93, 
                              Num_Group_Each_Sample = 3,
                              M = M, M_Dorf = M_Dorf, 
                              M_Dorf_opt_Se1 = M_Dorf_opt_Se1, 
                              M_Dorf_opt_Se0.95 = M_Dorf_opt_Se0.95, 
                              CT = CT, CT_prob = CT_prob, sigma = 1, 
                              cut_off_CT = 31.5, 
                              Spe_Amnt_mL = 0.8)
end_time <- proc.time()

################### Few names change in the paper ###########
# SBL is spare Bayesian Learning
# LS is least-squares
# P-Best_2021_reduced is now "P-BEST 2021*"
# P-Best_2021_full is now "P-BEST 2021"
# Regular testing or Array is now "standard retesting"
# Dorfman testing is now "regular Dorfman"
# ET.per is now E
# PPA is now PSe
# NPA is now PSp
# PPV is now PPPV
# NPV is now PNPV
#############################################################

# Results summary
My_results$res_summary[1:12, ]


# Flatting raw results
flat_raw_results <- as.data.frame.table(My_results$raw_results, responseName = "value")
# Writing the files in csv format
#write.csv(My_results$res_summary, file = "ResultsSummaryOn93Times961PoolingMatrix.Sigma1.p0.01.csv", row.names = F)
#write.csv(flat_raw_results, file = "RawResultsOn93Times961PoolingMatrix.Sigma1.p0.01.csv", row.names = F)
write.csv(My_results$res_summary, file = "ResultsSummaryOn93Times961PoolingMatrix.Sigma1.p0.02.csv", row.names = F)
write.csv(flat_raw_results, file = "RawResultsOn93Times961PoolingMatrix.Sigma1.p0.02.csv", row.names = F)

# Total time taken in hours
(end_time[3] - start_time[3])/3600
