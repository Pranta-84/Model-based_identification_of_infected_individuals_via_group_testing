# Function for simulating data
sim_data_ReddiLab <- function(data, M, M_Dorf, M_Dorf_opt_Se1, M_Dorf_opt_Se0.95,
                              Spe_Amnt_microL = 8, sigma = 1, cut_off_CT = 33){
  
  # True individual CT
  true_ind_CT <- data$Ct
  
  # true status of the specimens 
  true_status <- as.numeric(true_ind_CT < cut_off_CT)
  
  # Viral loads of the specimens
  VL <- data$VL*true_status
  
  # Group size
  Gr_Sz <- apply(M, 1, sum)
  
  # Number of groups
  Num_group <- dim(M)[1]
  
  # Place to store measured group viral loads
  measured_group_VL <- rep(NA, Num_group)
  
  # Error in group viral loads
  error_gr <- rnorm(Num_group, mean = 0, sd = sigma)
  
  for(i in 1:Num_group){
    
    # Specimen contribution
    Spe_contr <- Spe_Amnt_microL/Gr_Sz[i]
    
    # Measured group viral loads per micro liter of a specimen
    total_VL <- ((Spe_contr*M[i, ]%*%VL)/Spe_Amnt_microL)
    # Viral loads can not be negative 
    measured_group_VL[i] <- ifelse(total_VL == 0, 0,
                                   total_VL + error_gr[i])
  }
  
  # Number of rows in Dorfman pooling matrix
  Num_Group_Dorf <- nrow(M_Dorf)
  
  # Number of rows in Dorfman pooling matrix at Se = Sp = 1
  Num_Group_Dorf_opt_Se1 <- nrow(M_Dorf_opt_Se1)
  
  # Number of rows in Dorfman pooling matrix at Se = 0.95, Sp = 1
  Num_Group_Dorf_opt_Se0.95 <- nrow(M_Dorf_opt_Se0.95)
  
  # Place to store measured CT of groups for Dorfman
  measured_group_VL_Dorf <- rep(NA, Num_Group_Dorf)
  
  # Place to store measured CT of groups for optimal Dorfman at Se = Sp = 1
  measured_group_VL_Dorf_opt_Se1 <- rep(NA, Num_Group_Dorf_opt_Se1)
  
  # Place to store measured CT of groups for optimal Dorfman at Se = 0.95, Sp = 1
  measured_group_VL_Dorf_opt_Se0.95 <- rep(NA, Num_Group_Dorf_opt_Se0.95)
  
  # Group sizes in Pooling matrices for Dorfman
  all_Gr_Sz_Dorf <- apply(M_Dorf, 1, sum)
  
  # Group sizes in Pooling matrices for optimal Dorfman at Se = Sp = 1 
  all_Gr_Sz_Dorf_opt_Se1 <- apply(M_Dorf_opt_Se1, 1, sum)
  
  # Group sizes in Pooling matrices for optimal Dorfman at Se = 0.95, Sp = 1 
  all_Gr_Sz_Dorf_opt_Se0.95 <- apply(M_Dorf_opt_Se0.95, 1, sum)
  
  # Error with opt Dorf Se = Sp = 1
  error_gr_opt_Dorf_Se1 <- rnorm(n = Num_Group_Dorf_opt_Se1, mean = 0, sd = sigma)
  
  # Error with opt Dorf Se = 0.95, Sp = 1
  #error_gr_opt_Dorf_Se0.95 <- rnorm(n = Num_Group_Dorf_opt_Se0.95, mean = 0, sd = sigma)
  
  # Measured group CT with Dorfman pooling matrix with same group size 
  for(d in 1:Num_Group_Dorf){
    
    # Specimen contribution
    Spe_contr_Dorf <- Spe_Amnt_microL/all_Gr_Sz_Dorf[d]
    
    # Row of Dorfman pooling matrix
    M_Dorf_part <- M_Dorf[d, ]
    
    # Measured group viral loads per micro liter of a specimen
    total_VL_Dorf <- (Spe_contr_Dorf*M_Dorf_part%*%VL)/Spe_Amnt_microL
    
    # Viral loads can not be negative 
    measured_group_VL_Dorf[d] <- ifelse(total_VL_Dorf == 0, 0,
                                   total_VL_Dorf + error_gr[d])
  }
  
  for(d in 1:Num_Group_Dorf_opt_Se1){
    
    # Specimen contribution
    Spe_contr_Dorf_Se1 <- Spe_Amnt_microL/all_Gr_Sz_Dorf_opt_Se1[d]
    
    # Row of Dorfman pooling matrix at Se = Sp = 1
    M_Dorf_opt_Se1_part <- unlist(M_Dorf_opt_Se1[d, ])
    
    # Measured group viral loads per micro liter of a specimen
    total_VL_Dorf_Se1 <- (Spe_contr_Dorf_Se1*M_Dorf_opt_Se1_part%*%VL)/Spe_Amnt_microL
    
    # Viral loads can not be negative 
    measured_group_VL_Dorf_opt_Se1[d] <- ifelse(total_VL_Dorf_Se1 == 0, 0,
                                        total_VL_Dorf_Se1 + error_gr_opt_Dorf_Se1[d])
    
  }
  
  for(d in 1:Num_Group_Dorf_opt_Se0.95){
    
    # Specimen contribution
    Spe_contr_Dorf_Se0.95 <- Spe_Amnt_microL/all_Gr_Sz_Dorf_opt_Se0.95[d]
    
    # Row of Dorfman pooling matrix at Se = 0.95, Sp = 1
    M_Dorf_opt_Se0.95_part <- unlist(M_Dorf_opt_Se0.95[d, ])
    
    # Measured group viral loads per micro liter of a specimen
    total_VL_Dorf_Se0.95 <- (Spe_contr_Dorf_Se0.95*M_Dorf_opt_Se0.95_part%*%VL)/Spe_Amnt_microL
    
    # Viral loads can not be negative 
    measured_group_VL_Dorf_opt_Se0.95[d] <- ifelse(total_VL_Dorf_Se0.95 == 0, 0,
                                                total_VL_Dorf_Se0.95 + error_gr_opt_Dorf_Se1[d])
  }
  
  
  # Listings output
  list(measured_group_VL = measured_group_VL, measured_ind_VL = VL,
       measured_group_VL_Dorf = measured_group_VL_Dorf,
       measured_group_VL_Dorf_opt_Se1 = measured_group_VL_Dorf_opt_Se1,
       measured_group_VL_Dorf_opt_Se0.95 = measured_group_VL_Dorf_opt_Se0.95,
       true_status = true_status, true_ind_CT = true_ind_CT)
}

# Function for fitting model with VL
Model_VL <- function(measured_group_VL, Mpos_update, NoF, cut_off_ind){
  
  # Doing cross-validation with viral loads
  cv_model_VL <- cv.glmnet(as.matrix(Mpos_update), measured_group_VL, alpha = 1,
                           intercept = FALSE, grouped = FALSE, standardize = FALSE,
                           family = 'gaussian',
                           lower.limits = 0, nfolds = NoF,
                           type.measure = 'mse')
  
  # Lambda giving minimum MSE
  lmd_VL = cv_model_VL$lambda.min
  
  # Estimated individual viral load
  estimated_ind_VL <- round(coef(cv_model_VL, s = lmd_VL)[-1, ], 2)
  
  # test result from the model
  test_result_model_VL <- as.numeric(estimated_ind_VL >= cut_off_ind)
  
  # Listing outputs
  list(lmd_VL = lmd_VL, estimated_ind_VL = estimated_ind_VL,
       test_result_model_VL = test_result_model_VL)
}

# Function for Dorfman testing
Dorfman_test_ReddiLab <- function(measured_group_VL, cut_off_VL, cut_off_CT, 
                                  Num_Sample, true_ind_CT, M_Dorf, sigma){
  
  # Group test results
  group_test_Dorfman <- as.numeric(measured_group_VL >= cut_off_VL)
  
  # Number of groups tested positive
  num_group_pos_Dorf <- sum(group_test_Dorfman)
  
  # Case 1: When all the groups test negative
  if(num_group_pos_Dorf == 0){
    
    # Number of retest
    num_retest_D <- 0
    
    # Test result
    test_result_D <- rep(0, Num_Sample)
    
    # Retest needed
    retest_needed <- FALSE
    
  } else{
    
    # Place to store test results
    test_result_D <- rep(NA, Num_Sample)
    
    # Index of positive groups
    positive_Group <- which(group_test_Dorfman == 1)
    
    # M with only positive groups
    M_new_Dorfman <- M_Dorf[positive_Group, ]
    
    if(length(positive_Group) == 1){
      
      # Index of potential positive specimens if only one group is positive
      spe_num_pos_D <- which(M_new_Dorfman == 1)
      
    } else{
      
      # Index of potential positive specimens if more than one group is positive
      spe_num_pos_D <- which(colSums(M_new_Dorfman) == 1)
      
    }
    
    # Index of negative specimens
    spe_num_neg_D <- setdiff(c(1:Num_Sample), spe_num_pos_D)
    
    # True CT of potential positive specimens
    true_ind_CT_retest_D <- true_ind_CT[spe_num_pos_D]
    
    # Number of retest
    num_retest_D <- length(true_ind_CT_retest_D)
    
    # Error in CT upon retest
    error_retest <- rnorm(n = num_retest_D, mean = 0, sd = sigma)
    
    # CT of potential positives upon retest
    ind_CT_retest_D <- ifelse(true_ind_CT_retest_D >= cut_off_CT, 34,
                              true_ind_CT_retest_D + error_retest)
    
    # Test result of retested specimens 
    test_result_D_WithRetest <- as.numeric(ind_CT_retest_D < cut_off_CT)
    
    # Storing result of retested specimens
    test_result_D[spe_num_pos_D] <- test_result_D_WithRetest
    
    # Storing test results to sure negatives
    test_result_D[spe_num_neg_D] <- 0
    
    # Retest needed
    retest_needed <- TRUE
  }
  
  # listing output
  list(Num_retest = num_retest_D,
       test_result = test_result_D,
       retest_needed = retest_needed)
}

# Function for standard retesting
Regular_test_ReddiLab <- function(true_ind_CT, spe_num_pos, spe_num_neg,
                         sigma, cut_off_CT, Num_Sample){
  
  # Place to store test results of array
  test_result_regular <- rep(NA, Num_Sample)
  
  # True CT of potential positive specimens
  true_ind_CT_retest <- true_ind_CT[spe_num_pos]
  
  # Number of retest with regular array testing
  num_retest_regular <- length(true_ind_CT_retest)
  
  # CT of potential positives upon retest
  ind_CT_retest <- ifelse(true_ind_CT_retest >= cut_off_CT, 34,
                          true_ind_CT_retest
                          + rnorm(n = num_retest_regular,
                                  mean = 0, sd = sigma))
  
  # Test results with retesting potential positives
  test_result_regular_WithRetest <- as.numeric(ind_CT_retest < cut_off_CT)
  
  # Result of retested specimens
  test_result_regular[spe_num_pos] <- test_result_regular_WithRetest
  
  # Assigning test results to sure negatives
  test_result_regular[spe_num_neg] <- 0
  
  # Listing outputs
  list(test_result_regular = test_result_regular,
       num_retest_regular = num_retest_regular,
       test_result_regular_WithRetest = test_result_regular_WithRetest)
}

# Function for running P-BEST 2021
PBest_test_ReddiLab <- function(measured_group_VL_PBest_Binary, M){
  
  # Direct translation of MATLAB code of P-Best authors to R code
  qMeasurement <- measured_group_VL_PBest_Binary
  
  maxNum <-  20
  
  dt = max(abs(t(M)%*%qMeasurement))
  
  tau = 0.005*dt;
  
  m <- length(qMeasurement)
  
  M_mat <- matlab_types$double(as.vector(as.double(M)),   
                               size = c(nrow(M), ncol(M)))
  
  qy <- matlab_types$double(as.list(as.numeric(qMeasurement)))
  
  qy <- eng$reshape(qy, as.integer(m), as.integer(1), nargout = as.integer(1))
  
  u = eng$opm(qy, M_mat, tau, as.numeric(maxNum))
  
  discreteOutput = eng$selectByError(u, M_mat, qy)
  
  detected_samples = eng$find(discreteOutput, nargout = as.integer(1))
  
  # Check if it is already a plain R numeric/double
  if(is.numeric(detected_samples) || is.double(detected_samples)){
    # Already an R number, just convert directly
    detected_idx <- as.integer(detected_samples)
  } else {
    # It is a Python object, extract _data as before
    ds_data <- reticulate::py_get_attr(detected_samples, "_data")
    detected_idx <- as.integer(unlist(reticulate::py_to_r(reticulate::iterate(ds_data))))
  }
  
  # Listing outputs
  list(detected_index = detected_idx, lambda = tau)
}

# Function to get results of all the algorithms
res_All_Algo_ReddiLab <- function(active_work_data, Num_Sample, Num_Group, Num_Group_Each_Sample,
                         M, M_Dorf, M_Dorf_opt_Se1, M_Dorf_opt_Se0.95,
                         sigma, cut_off_CT, cut_off_VL = 0.1590,
                         Spe_Amnt_microL = 8){
  # Simulating data
  sim_data <- sim_data_ReddiLab(data = active_work_data, M = M, M_Dorf = M_Dorf,
                                M_Dorf_opt_Se1 = M_Dorf_opt_Se1,
                                M_Dorf_opt_Se0.95 = M_Dorf_opt_Se0.95,
                                Spe_Amnt_microL = Spe_Amnt_microL,
                                sigma = sigma,
                                cut_off_CT = cut_off_CT)
  
  # Results of Dorfman testing
  res_Dorfman <- Dorfman_test_ReddiLab(measured_group_VL = sim_data$measured_group_VL_Dorf,
                                       cut_off_VL = cut_off_VL,
                                       cut_off_CT = cut_off_CT,
                                       Num_Sample = Num_Sample,
                                       true_ind_CT = sim_data$true_ind_CT,
                                       M_Dorf = M_Dorf, sigma = sigma)
  
  # Results of optimal Dorfman testing at Se = Sp = 1
  res_Dorfman_opt_Se1 <- Dorfman_test_ReddiLab(measured_group_VL = sim_data$measured_group_VL_Dorf_opt_Se1,
                                               cut_off_VL = cut_off_VL, cut_off_CT = cut_off_CT,
                                               Num_Sample = Num_Sample, true_ind_CT = sim_data$true_ind_CT,
                                               M_Dorf = M_Dorf_opt_Se1, sigma = sigma)
  
  # Results of optimal Dorfman testing at Se = 0.95, Sp = 1
  res_Dorfman_opt_Se0.95 <- Dorfman_test_ReddiLab(measured_group_VL = sim_data$measured_group_VL_Dorf_opt_Se0.95,
                                                  cut_off_VL = cut_off_VL, cut_off_CT = cut_off_CT,
                                                  Num_Sample = Num_Sample, true_ind_CT = sim_data$true_ind_CT,
                                                  M_Dorf = M_Dorf_opt_Se0.95, sigma = sigma)
  
  ################################## P-Best 2024 ###########################################
  # Cutoff to decide on whether the viral load of the group is large
  cut_off_VL_large <- 1
  # Result of PBest 2024
  res_PBest_2024 <- PBest_2024_ReddiLab(true_ind_CT = sim_data$true_ind_CT,
                               measured_Group_VL =  sim_data$measured_group_VL, M = M,
                               Num_Sample = Num_Sample, cut_off_VL_large = cut_off_VL_large,
                               Num_Group_Each_Sample = Num_Group_Each_Sample, sigma = sigma,
                               cut_off_CT = cut_off_CT)
  
  ################################## P-Best with Full pooling matrix ##########################
  # Place to store PBest's results
  test_result_PBest_full <- rep(NA, Num_Sample)
  
  # Binarizing measured group CT
  measured_group_VL_PBest_Binary <- as.numeric(sim_data$measured_group_VL >= cut_off_VL)
  
  # Running P-Best
  res_PBest <- PBest_test_ReddiLab(measured_group_VL_PBest_Binary = measured_group_VL_PBest_Binary, M = M)
  
  # Regularization parameter
  PBest_full_lambda <- res_PBest$lambda
  
  # Storing test results of detected specimens 
  test_result_PBest_full[res_PBest$detected_index] <- 1
  
  # Index of non-detected specimens 
  nondetected_index <- setdiff(c(1:Num_Sample), res_PBest$detected_index)
  
  # Storing test results of non-detected specimens
  test_result_PBest_full[nondetected_index] <- 0
  
  # Number of specimens PBest predicted to be positive using full pooling matrix
  num_specimens_PBest_full_predicted_positive <- sum(test_result_PBest_full)
  # Number of specimens PBest predicted to be negative using full pooling matrix
  num_specimens_PBest_full_predicted_negative <- sum(test_result_PBest_full == 0)
  
  #############################################################################################
  # Place to store test results of algorithms
  test_result <- test_result_regular_AR <- rep(NA, Num_Sample)
  
  # Group test results
  group_result <- as.numeric(sim_data$measured_group_VL >= cut_off_VL)
  
  # Number of groups tested positive
  group_result_sum <- sum(group_result)
  
  if(group_result_sum == 0){
    
    # Test results
    test_result <- test_result_regular_AR <- rep(0, Num_Sample)
    
    # Number of retest
    num_retest <- num_retest_regular <- num_retest_regular_AR <- 0
    
    # Indicates whether model was used
    Model_Used <- FALSE
    
    
  } else if(group_result_sum != 0 & group_result_sum < Num_Group_Each_Sample){
    
    # Index of positive tested rows
    pos_index <- which(group_result == 1)
    
    if(length(pos_index) == 1){
      
      # Index of potential positives
      spe_num <- which(M[pos_index, ] == 1)
      
    } else{
      
      # Rows of M corresponding to positives
      M_pos <- M[pos_index, ]
      
      # Index of potential positives
      spe_num <- which(colSums(M_pos) >= 1)
    }
    
    # Index of sure negative specimens
    index_negative <- setdiff(c(1:Num_Sample), spe_num)
    
    # Assigning to test result to sure negatives
    test_result[index_negative] <- 0
    
    # True CT of potential positive specimens
    true_ind_CT_retest <- sim_data$true_ind_CT[spe_num]
    
    # Number of retests
    num_retest <- num_retest_regular <- num_retest_regular_AR <- length(true_ind_CT_retest)
    
    # Error upon retest
    error_retest <- rnorm(n = num_retest, mean = 0, sd = sigma)
    
    # CT of potential positives upon retest
    ind_CT_retest <- ifelse(true_ind_CT_retest >= cut_off_CT, 34,
                            true_ind_CT_retest + error_retest)
    
    # Result of retested specimens
    test_result_WR <- as.numeric(ind_CT_retest < cut_off_CT)
    test_result[spe_num] <- test_result_WR
    
    # Test result of regular testing with always retest
    test_result_regular_AR <- test_result
    
    # Indicates whether model was used
    Model_Used <- FALSE
    
  } else {
    
    # Which groups are positive
    pos_index <- which(group_result == 1)
    
    # Pooling matrix where all groups are positive
    Mpos <- M[pos_index, ]
    
    # Column sum of the pooling matrix with all positive groups
    col_sum <- colSums(Mpos)
    
    # Index of possible positives
    spe_num_pos <- which(col_sum >= Num_Group_Each_Sample)
    
    # Index of sure positives
    spe_num_neg <- setdiff(c(1:Num_Sample), spe_num_pos)
    
    ############################################################################
    # Steps followed to decide whether the reduced pooling matrix has ambiguity
    if(length(spe_num_pos) == 1){
      
      # Indicator whether need to use model
      ambi_spe_indic <- 0
      
    } else{
      
      # Pooling matrix with all positive rows and all potential positives
      Mpos_update <- Mpos[, spe_num_pos]
      
      # Which row sums are equal to 1
      sum1_rows <- which(rowSums(Mpos_update) == 1)
      
      if(length(sum1_rows) < length(spe_num_pos)){
        # assigning a number greater than 0 to indicate model needed
        ambi_spe_indic <- 2
        
      } else{
        
        # Indicator of whether need to use model
        ambi_spe_indic <- sum(apply(Mpos_update[sum1_rows, ], 2, sum) == 0)
        
      }
    }
    
    ###########################################################################
    # If no ambiguity
    if(ambi_spe_indic == 0){
      
      # Test results
      test_result[spe_num_pos] <- 1
      
      test_result[spe_num_neg] <- 0
      
      # Number of retest
      num_retest <- num_retest_regular <- 0
      
      # Whether model used
      Model_Used <- FALSE
      
      ################## Array with always retest #######################
      test_result_regular_AR[spe_num_neg] <- 0
      
      # True CT of potential positive specimens
      true_ind_CT_retest <- sim_data$true_ind_CT[spe_num_pos]
      
      # Number of retests
      num_retest_regular_AR <- length(true_ind_CT_retest)
      
      # Error in CT upon retest
      error_retest <- rnorm(n = num_retest_regular_AR, mean = 0, sd = sigma)
      
      # CT of potential positives upon retest
      ind_CT_retest_WR <- ifelse(true_ind_CT_retest >= cut_off_CT, 34,
                                 true_ind_CT_retest + error_retest)
      
      test_result_WR <- as.numeric(ind_CT_retest_WR < cut_off_CT)
      
      # Result of retested specimens
      test_result_regular_AR[spe_num_pos] <- test_result_WR
      
    } else {
      
      # Model usage
      Model_Used <- TRUE
      
      # Number of retest
      num_retest <- 0
      
      # Removing row whose all the values are zero
      redundant_rows <- which(rowSums(Mpos_update) == 0)
      non_redunt_rows <- setdiff(c(1:nrow(Mpos_update)), redundant_rows)
      Mpos_update2 <- Mpos_update[non_redunt_rows, ]
      
      # Number of rows in trimmed M
      num_rows_Trim_M <- nrow(Mpos_update2)
      
      # NoF stands for number of folds
      NoF <- ifelse(num_rows_Trim_M < 5, num_rows_Trim_M, 5)
      
      # Measured group VL
      measured_group_VL <- sim_data$measured_group_VL
      
      # Updated positive index
      update_pos_idx <- pos_index[non_redunt_rows]
      
      # Measured group viral load of only positive groups
      measured_group_VL_onlyPos <- measured_group_VL[update_pos_idx]
      
      # Individual viral load cutoff
      cut_off_ind <- cut_off_VL
      
      ############################ Lasso VL #####################################
      # Place to store test results
      test_result_Lasso_VL <- rep(NA, Num_Sample)
      
      # Fitting model
      res_model_VL <- Model_VL(measured_group_VL = measured_group_VL_onlyPos,
                               Mpos_update = Mpos_update2, NoF = NoF,
                               cut_off_ind = cut_off_ind)
      # Storing test results
      test_result_Lasso_VL[spe_num_pos] <- res_model_VL$test_result_model_VL
      test_result_Lasso_VL[spe_num_neg] <- 0
      
      # Regularization parameter value
      lmd_Lasso_VL <- res_model_VL$lmd_VL
      
      ################################ SBL VL####################################
      test_result_SBL_VL <- rep(NA, Num_Sample)
      # Fitting SBL
      estimated_VL_SBL <- sbl$sbl(as.matrix(Mpos_update2), measured_group_VL_onlyPos)
      # Storing test results
      test_result_SBL <- as.numeric(estimated_VL_SBL >= cut_off_ind)
      test_result_SBL_VL[spe_num_pos] <- test_result_SBL
      test_result_SBL_VL[spe_num_neg] <- 0
      
      ################################ l1LS VL ###################################
      test_result_LS_VL <- rep(NA, Num_Sample)
      # Fitting Least square with L1 regularization
      GVL <- l1ls$np$array(measured_group_VL_onlyPos)
      estimated_VL_l1ls <- l1ls$l1ls(as.matrix(Mpos_update2), GVL, 0.10, 0)
      
      # Storing results
      test_result_LS <- as.numeric(estimated_VL_l1ls >= cut_off_ind)
      test_result_LS_VL[spe_num_pos] <- test_result_LS
      test_result_LS_VL[spe_num_neg] <- 0
      
      ################################# NNOMP VL #####################################
      test_result_NNOMP_VL <- rep(NA, Num_Sample)
      # Fitting NNOMP
      estimated_VL_NNOMP <- nnompcv$nnomp(Mpos_update2, as.integer(0), measured_group_VL_onlyPos,
                                          as.integer(0), as.integer(num_rows_Trim_M))
      
      # Storing test results
      test_result_NNOMP <- as.numeric(estimated_VL_NNOMP >= cut_off_ind)
      test_result_NNOMP_VL[spe_num_pos] <- test_result_NNOMP
      test_result_NNOMP_VL[spe_num_neg] <- 0
      
      ######################## P-Best with reduced pooling matrix #################
      test_result_PBest_reduced <- rep(NA, Num_Sample)
      
      # Direct translation of MATLAB code of P-Best authors
      measured_group_VL_PBest_Binary_onlyPos <- as.numeric(measured_group_VL_onlyPos >= cut_off_VL)
      
      # Fitting P-Best with reduced pooling matrix
      res_PBest_reduced <- PBest_test_ReddiLab(measured_group_VL_PBest_Binary = measured_group_VL_PBest_Binary_onlyPos, 
                                      M = Mpos_update2)
      
      # Regularization parameter value for PBest with reduced pooling matrix
      PBest_reduced_lambda <- res_PBest_reduced$lambda
      
      # Index of positive predicted specimens 
      pos_predicted <- spe_num_pos[res_PBest_reduced$detected_index]
      test_result_PBest_reduced[pos_predicted] <- 1
      
      # Index of negative predicted specimens
      neg_predicted <-  setdiff(spe_num_pos, pos_predicted)
      
      # Assigning test results
      test_result_PBest_reduced[neg_predicted] <-  0
      test_result_PBest_reduced[spe_num_neg] <- 0
      
      test_result_pos_model_PBest_reduced <- length(pos_predicted)
      test_result_neg_model_PBest_reduced <- length(neg_predicted)
      
      ##################### Regular array testing starts ################
      res_regular <- Regular_test_ReddiLab(true_ind_CT = sim_data$true_ind_CT,
                                           spe_num_pos = spe_num_pos,
                                           spe_num_neg = spe_num_neg, sigma = sigma,
                                           cut_off_CT = cut_off_CT,
                                           Num_Sample = Num_Sample)
      # number of retests
      num_retest_regular <- num_retest_regular_AR <- res_regular$num_retest_regular
      
    }
  }
  
  if(Model_Used){
    # Number of specimens needed model prediction
    num_specimens_needed_model_prediction <- length(spe_num_pos)
    
    # Lasso VL
    num_specimens_Lasso_predicted_positive_VL <- sum(res_model_VL$test_result_model_VL == 1)
    num_specimens_Lasso_predicted_negative_VL <- sum(res_model_VL$test_result_model_VL == 0)
    test_result_Lasso_VL <- test_result_Lasso_VL
    lmd_Lasso_VL <- lmd_Lasso_VL
    
    ## SBL VL
    num_specimens_SBL_predicted_positive_VL <- sum(test_result_SBL == 1)
    num_specimens_SBL_predicted_negative_VL <- sum(test_result_SBL == 0)
    test_result_SBL_VL <- test_result_SBL_VL
    
    ## Least Squares VL
    num_specimens_LS_predicted_positive_VL <- sum(test_result_LS == 1)
    num_specimens_LS_predicted_negative_VL <- sum(test_result_LS == 0)
    test_result_LS_VL <- test_result_LS_VL
    
    ## NNOMP VL
    num_specimens_NNOMP_predicted_positive_VL <- sum(test_result_NNOMP == 1)
    num_specimens_NNOMP_predicted_negative_VL <- sum(test_result_NNOMP == 0)
    test_result_NNOMP_VL <-test_result_NNOMP_VL
  
    ## P-Best reduced
    num_specimens_PBest_reduced_predicted_positive <- test_result_pos_model_PBest_reduced
    num_specimens_PBest_reduced_predicted_negative <- test_result_neg_model_PBest_reduced
    test_result_PBest_reduced <- test_result_PBest_reduced
    PBest_reduced_lambda <- PBest_reduced_lambda
    
    ## Regular
    num_specimens_tested_positive_regular <- sum(res_regular$test_result_regular_WithRetest == 1)
    num_specimens_tested_negative_regular <- sum(res_regular$test_result_regular_WithRetest == 0)
    
    test_result_regular <- test_result_regular_AR <- res_regular$test_result_regular
    
    
  } else{
    # Number of specimens needed model prediction
    num_specimens_needed_model_prediction <- NA
    
    # Lasso VL
    num_specimens_Lasso_predicted_positive_VL <- num_specimens_Lasso_predicted_negative_VL <- lmd_Lasso_VL <- NA
    
    ## SBL VL
    num_specimens_SBL_predicted_positive_VL <- num_specimens_SBL_predicted_negative_VL <- NA
    
    ## LS VL
    num_specimens_LS_predicted_positive_VL <- num_specimens_LS_predicted_negative_VL <- NA
    
    ## NNOMP VL
    num_specimens_NNOMP_predicted_positive_VL <- num_specimens_NNOMP_predicted_negative_VL <- NA
    
    # P-Best reduced
    num_specimens_PBest_reduced_predicted_positive <- num_specimens_PBest_reduced_predicted_negative <- NA
    PBest_reduced_lambda <- NA
    
    # Array
    num_specimens_tested_positive_regular <- num_specimens_tested_negative_regular <- NA
    
    # Test results from algorithms
    test_result_Lasso_VL <- test_result_SBL_VL <- test_result_LS_VL <- test_result_NNOMP_VL <-
      test_result_PBest_reduced <- test_result_regular <- test_result
    test_result_regular_AR <- test_result_regular_AR
  }
  
  # listing output
  list(true_status = sim_data$true_status,
       Model_Used = Model_Used,
       test_result_Lasso_VL = test_result_Lasso_VL,
       test_result_SBL_VL = test_result_SBL_VL,
       test_result_LS_VL = test_result_LS_VL,
       test_result_NNOMP_VL = test_result_NNOMP_VL,
       test_result_PBest_reduced = test_result_PBest_reduced,
       test_result_regular = test_result_regular,
       test_result_regular_AR = test_result_regular_AR,
       res_Dorfman = res_Dorfman,
       res_Dorfman_opt_Se1 = res_Dorfman_opt_Se1,
       res_Dorfman_opt_Se0.95 = res_Dorfman_opt_Se0.95,
       test_result_PBest_full = test_result_PBest_full,
       PBest_full_lambda = PBest_full_lambda,
       lmd_Lasso_VL = lmd_Lasso_VL,
       PBest_reduced_lambda = PBest_reduced_lambda,
       num_retest = num_retest,
       res_PBest_2024 = res_PBest_2024,
       num_retest_regular = num_retest_regular,
       num_retest_regular_AR = num_retest_regular_AR,
       num_specimens_needed_model_prediction = num_specimens_needed_model_prediction,
       num_specimens_Lasso_predicted_positive_VL = num_specimens_Lasso_predicted_positive_VL,
       num_specimens_Lasso_predicted_negative_VL = num_specimens_Lasso_predicted_negative_VL,
       num_specimens_SBL_predicted_positive_VL = num_specimens_SBL_predicted_positive_VL,
       num_specimens_SBL_predicted_negative_VL = num_specimens_SBL_predicted_negative_VL,
       num_specimens_LS_predicted_positive_VL = num_specimens_LS_predicted_positive_VL,
       num_specimens_LS_predicted_negative_VL = num_specimens_LS_predicted_negative_VL,
       num_specimens_NNOMP_predicted_positive_VL = num_specimens_NNOMP_predicted_positive_VL,
       num_specimens_NNOMP_predicted_negative_VL = num_specimens_NNOMP_predicted_negative_VL,
       num_specimens_PBest_reduced_predicted_positive = num_specimens_PBest_reduced_predicted_positive,
       num_specimens_PBest_reduced_predicted_negative = num_specimens_PBest_reduced_predicted_negative,
       num_specimens_PBest_full_predicted_positive = num_specimens_PBest_full_predicted_positive ,
       num_specimens_PBest_full_predicted_negative = num_specimens_PBest_full_predicted_negative ,
       num_specimens_tested_positive_regular = num_specimens_tested_positive_regular,
       num_specimens_tested_negative_regular = num_specimens_tested_negative_regular)
}


# Function to calculate accuracy of results
PSe <- function(status, result){
  # status is true status of each individual
  # result is the tested result of each individual
  
  # Number of tested positive among true positives
  tp <- sum(result[status == 1])
  fn <- sum(result[status == 1] == 0)
  
  # Number of tested negative among true negatives
  tn <- sum(result[status == 0] == 0)
  fp <- sum(result[status == 0] == 1)
  # ntn and ntp are number of true negative and positive
  list(TestNeg = tn, num_true_neg = sum(status == 0),
       TestPos = tp, num_true_pos = sum(status == 1),
       fp = fp, fn = fn) 
}

# Function for implementing P-BEST 2024
PBest_2024_ReddiLab <- function(true_ind_CT, measured_Group_VL, M,
                       Num_Sample, cut_off_VL_large = 1,
                       Num_Group_Each_Sample,
                       sigma = 1,
                       cut_off_CT = 33){
  
  # Making sure there is no negative viral loads
  measured_Group_VL <- ifelse(measured_Group_VL <= 0,
                              0, measured_Group_VL)
  
  # Direct translation of MATLAB code of P-Best authors to R code
  qMeasurement <- measured_Group_VL
  
  # Maximum group viral loads
  dt = max(abs(t(M)%*%qMeasurement))
  
  # Regularization parameter
  tau = 0.005*dt
  
  # Number of groups
  m <- length(qMeasurement)
  
  M_mat <- matlab_types$double(as.vector(as.double(M)),   # <-- Fix 1
                               size = c(nrow(M), ncol(M)))
  #qy    <- matlab_types$double(as.numeric(qMeasurement),  # <-- Fix 4
    #                           size = c(m, 1L))
  
  # Making object of group viral loads compatible for MATLAB
  qy <- matlab_types$double(as.list(as.numeric(qMeasurement)))
  qy <- eng$reshape(qy, as.integer(m), as.integer(1), nargout = as.integer(1))
  
  # Fitting GPSR to estimate individual viral loads
  estimate_VL <- eng$applyGPSR(qy, M_mat, tau)
  
  # Changing matlab.double object to R object
  estimate_VL_r <- as.numeric(np$array(estimate_VL))
  
  # Index of potential positive specimens 
  index_pot_pos <- which(estimate_VL_r > 0)
  
  # Index of negative specimens
  index_spe_neg <- setdiff(c(1:Num_Sample),
                           index_pot_pos)
  
  # If there is at least one potential positive
  if(!is.null(index_pot_pos)){
    
    # Reducing M with taking the columns of potential positive specimens
    M_new <- M[, index_pot_pos]
    
    if(length(index_pot_pos) >= 2){
      
      # Criterion 1: There must at least one group where it is the only potential positive
      # Taking row sum of the pooling matrix with positive estimated viral loads
      row_sum_M_new <- apply(M_new, 1, sum)
      # Finding the rows with only one specimen
      sum1_row <- which(row_sum_M_new == 1)
      
      # Place to store index of specimens that are in at least one group with being potential positive
      index_crit1 <- c()
      
      for(rs in sum1_row){
        # Index of the specimen  
        lone_spe <- which(M_new[rs, ] == 1)
        # Adding index of the specimen to the list
        index_crit1 <- c(index_crit1, lone_spe)
      }
      # Index of the specimens in M_new that satisfies criterion 1
      index_crit1 <- unique(index_crit1)
      # Index of the specimens in M which satisfies criterion 1
      true_index_crit1 <- index_pot_pos[index_crit1]
      # Index of specimens that does not satisfy criterion 1 
      suspected <- setdiff(index_pot_pos, true_index_crit1)
      
    } else{
      
      # Criterion 1: There must at least one group where it is the only potential positive
      # Index of specimens that are in at least one group with being potential positive
      index_crit1 <- true_index_crit1 <- index_pot_pos
      
      # Index of specimens that does not satisfy criterion 1 
      suspected <- c()
    }
    
    
    # If there is at least one potential positive satisfying criterion 1
    # Check criterion 2: All the groups the specimen is in has large viral loads
    if(!is.null(index_crit1)){
      # Place to store index of predicted positive specimens
      pred_pos <- c()
      # Checking the viral loads of all the groups each specimen satisfying criterion 1 is in; need all of them to be large to ensure criterion 2
      for(pc in true_index_crit1){
        # Finding the index of all the groups specimen is in
        spe_in_gr <- which(M[, pc] == 1)
        # Viral loads of all the groups specimen is in
        VL_spe_in_gr <- measured_Group_VL[spe_in_gr]
        # Number of groups with large viral loads
        num_large_VL_gr <- sum(VL_spe_in_gr >= cut_off_VL_large)
        # If number of groups with large viral loads is less than number of groups/sample
        if(num_large_VL_gr < Num_Group_Each_Sample){
          # Add the specimen to suspected list
          suspected <- c(suspected, pc)
          # Else criterion 2 is satisfied add the specimen to the list of predicted positives  
        } else pred_pos <- c(pred_pos, pc)
      }
    } else{
      # Index of predicted positives
      pred_pos <- c()
      # Index of suspected specimens
      suspected <- suspected
    } 
    
    # Number of specimens satisfying criterion 1 and 2
    len_pred_pos <- length(pred_pos)
    
    # If there is a specimen which satisfied criterion 1 and 2
    if(len_pred_pos >= 1){
      # True status of the specimens satisfying criterion 1 and 2
      true_status_pred_pos <- as.numeric(true_ind_CT[pred_pos] < cut_off_CT)
      # Test results
      test_result_pred_pos <- rep(1, len_pred_pos)
      # How many of the specimens satisfying criterion 1 and 2 are correctly tested positive 
      correct_pos <- PSe(status = true_status_pred_pos, result = test_result_pred_pos)$TestPos
      # else correct_pos is NA  
    } else correct_pos <- NA
    
    # Place to store test results
    test_result <- rep(NA, Num_Sample)
    # Assigning test results to specimens predicted positive
    test_result[pred_pos] <- 1
    # Assigning test results to specimens predicted negative
    test_result[index_spe_neg] <- 0
    # Number of retest
    Num_retest <- length(suspected)
    # True CT of specimens that did not satisfy criterion 1 or 2
    true_CT_suspected <- true_ind_CT[suspected]
    # Error in CT upon retest
    error_retest <- rnorm(Num_retest, 0, sigma)
    # True CT of suspected specimens upon retest
    true_CT_suspected_upon_retest <- ifelse(true_CT_suspected >= cut_off_CT,
                                            34,
                                            true_CT_suspected + error_retest)
    # Assigning test results of suspected specimens
    test_result[suspected] <- as.numeric(true_CT_suspected_upon_retest < cut_off_CT)
    
    # Else
  } else {
    
    test_result <- rep(0, Num_Sample)
    Num_retest <- 0
    len_pred_pos <- NA
    correct_pos <- NA
  }
  
  # Listing outputs
  list(test_result = test_result,
       Num_retest = Num_retest,
       len_pred_pos = len_pred_pos,
       correct_pos = correct_pos,
       lambda = tau)
}

################### Few names change in the paper ###########
# SBL is spare Bayesian Learning
# LS is least-squares
# P-Best_reduced is "P-BEST 2021*"
# P-Best_full is "P-BEST 2021"
# Regular testing or Array is "standard retesting"
#############################################################


