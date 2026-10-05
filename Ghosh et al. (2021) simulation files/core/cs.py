# vim: tabstop=2 expandtab shiftwidth=2 softtabstop=8
import numpy as np
import math
from sklearn.linear_model import Lasso, LassoLars, LassoCV, LassoLarsCV
import pylops
from joblib import Parallel, delayed
import pandas as pd
from core.comp import create_infection_array_with_num_cases, COMP
from inbuilt_algos import nnompcv
from inbuilt_algos import sbl
from inbuilt_algos import l1ls
import algos
## Importing matlab.engine to use matlab functions
import matlab.engine
import struct
# Starting matlab
eng = matlab.engine.start_matlab()
# Giving path of MATLAB functions
eng.addpath("/home/bilder/pranta/tapestry/PBestMatlabCode/", nargout = int(0))

from utils import output_validation_utils

from core import config

from core.matrices import *

# Numpy configuration
np.set_printoptions(precision=3)
# Numpy should raise Exception on division by zero so that we can catch programming errors
np.seterr(all='raise')

# Use compressed sensing to solve 0.5*||Mx - y||^2 + l * ||x||_1
class CS(COMP):
  def __init__(self, n, t, s, d, l, arr, M=None, mr=None):
    super().__init__(n, t, s, d, arr)
    if M is not None:
      assert n == M.shape[1]
      assert t == M.shape[0]
      self.M = M.T
      #print(self.M.shape)
    self.create_conc_matrix_from_infection_array(arr)
    self.l = l
    self.mr = mr
    # Pooling matrix to perform Dorfman testing
    self.DT = pd.read_csv("Dorf_Pool_matrix.csv").to_numpy()
    # This is the pooling matrix of optinal Dorfman testing at Se = Sp = 1
    # For different value of number positives per data this needs to be changed
    self.DT_opt_Se1 = pd.read_csv("DorfMatrix_p15by961_SeSp1.csv").to_numpy()
    
    # This is the pooling matrix of optinal Dorfman testing at Se = 0.95 and Sp = 1
    # For different value of number positives per data this needs to be changed
    self.DT_opt_Se95by100 = pd.read_csv("DorfMatrix_p15by961_SeSp1.csv").to_numpy()
  
  # Multiply actual conc matrix to M. This is the part done by mixing samples
  # and qpcr
  def get_quantitative_results(self, conc, add_noise=False,
      noise_magnitude=None, noise_model='exponential_gaussian'):
    conc = np.expand_dims(conc, axis=-1)
    #print(self.M.shape, conc.shape)
    y = np.matmul(self.M.T, conc).flatten()
    # Group viral loads when doing Dorfman testing with the same group size as with the pooling matrix
    y_dorf = np.matmul(self.DT, conc).flatten()
    
    # Group viral loads when doing Dorfman testing with optimal group size at Sensitivity = 1, Specificity = 1
    y_opt_dorf_Se1 = np.matmul(self.DT_opt_Se1, conc).flatten()
    # Group viral loads when doing Dorfman testing with optimal group size at Sensitivity = 0.95, Specificity = 1
    y_opt_dorf_Se95by100 = np.matmul(self.DT_opt_Se95by100, conc).flatten()
    
    # print('lol')
    sigval = 0.
    if add_noise:
      if noise_magnitude is not None:
        # This error is independent of the magnitude
        sigval = noise_magnitude
        error = np.random.normal(0., sigval)
        y = y + error
        raise ValueError("This noise model is incorrect and hence disabled. "
            "Enable this if you know what you're doing")
      elif config.noise_model == 'variable_gaussian':
        # This error is proportional to magnitude of y
        sigval = 0.01*np.absolute(y)
        error = np.random.normal(0., sigval)
        y = y + error
      elif config.noise_model == 'exponential_gaussian':
        # This noise accounts for cycle time variability
        # Cycle time is assumed to be Gaussian distributed, due to which log
        # of y is Gaussian. Hence 
        #p = 0.95
        error = np.random.normal(0., config.eps_std_dev, size=self.t)
        # Error Dorfman testing
        error_DT = np.random.normal(0., config.eps_std_dev, size=self.DT.shape[0])
        # Error with optimal Dorfman testing at Sensitivity = 1, Specificity = 1
        error_DT_opt_Se1 = np.random.normal(0., config.eps_std_dev, size=self.DT_opt_Se1.shape[0])
        # Error with optimal Dorfman testing at Sensitivity = 0.95, Specificity = 1
        error_DT_opt_Se95by100 = np.random.normal(0., config.eps_std_dev, size=self.DT_opt_Se95by100.shape[0])
        
        #print('Original y', y)
        #print('error exponents', error)
        y = y * ((1 + config.p) ** error)
        # Group viral loads with Dofman testing after adjusting for error
        y_dorf = y_dorf*((1 + config.p)**error_DT)
        # Group viral loads with optimal Dorfman testing at Sensitivity = 1, Specificity = 1 after adjusting for error
        y_opt_dorf_Se1 = y_opt_dorf_Se1*((1 + config.p)**error_DT_opt_Se1)
        # Group viral load with optimal Dorfman testing at Sensitivity = 0.95, Specificity = 1 after adjusting for error
        y_opt_dorf_Se95by100 = y_opt_dorf_Se95by100*((1 + config.p)**error_DT_opt_Se95by100)
        
      else:
        raise ValueError('Invalid noise model %s' % noise_model)


      if config.bit_flip_prob > 0 :
        raise ValueError('This is probably a mistake')
        print('before flip y = ', y)
        mask = (y > 0).astype(np.int32)
        flip = np.random.binomial(1, config.bit_flip_prob, self.t)
        flip = flip * mask
        y = y * (1 - flip)
        print('after flip y = ', y)
    # Returning outputs
    return y, y_dorf, y_opt_dorf_Se1, y_opt_dorf_Se95by100, sigval

  # Initial concentration of RNA in each sample
  def create_conc_matrix_from_infection_array(self, arr):
    # Fix tau to 0.01 * minimum value we expect in x
    # XXX: Actually tau should be defined by the algorithm.
    self.tau = 0.01 * 1 / config.scale
    #self.tau = 0.01 * 0.1
    #conc = 1 + np.random.poisson(lam=5, size=self.n)
    #conc = np.random.randint(low=config.x_low * config.scale, high=config.x_high *
    # config.scale + 1, size=self.n) / config.scale
    conc = np.random.uniform(config.x_low, config.x_high, size=self.n)
    #conc = 0.1 + 0.9 * np.random.rand(self.n)
    #conc = np.random.randint(low=1, high=11, size=self.n) / 10.
    #conc = np.ones(self.n)
    self.conc = conc * arr # Only keep those entries which are non-zero in arr

  # Solve the CS problem using Lasso
  #
  # y is results
  def decode_lasso(self, results, algo='lasso', prefer_recall=False,
      compute_stats=True):
    determined = 0
    overdetermined = 0
    # Add check if system is determined or overdetermined
    if self.t == self.n:
      determined = 1
    elif self.t > self.n:
      overdetermined = 1

    prob1 = None
    prob0 = None
    answer_high_precision = np.zeros(self.n)
    if algo == 'lasso':
      #lasso = LassoLars(alpha=self.l)
      #lasso = Lasso(alpha=self.l, positive = True, max_iter=10000)
      lasso = LassoCV(n_alphas=100, positive = True, max_iter=10000)
      lasso.fit(self.M.T, results)
      
      # Regularization parameter value
      lamda = lasso.alpha_
      
      # Estimated viral loads
      answer = lasso.coef_
      # Covering estimated viral loads to binary using positive cutoff 0
      infected = (answer != 0.).astype(np.int32)
      # Number of retest
      retest = 0
      
      # Checking ambiguity
      y = results
      Ambi = self.AmbiCheck(y)
      
      score = np.linalg.norm(answer - self.conc) / math.sqrt(self.t)

      if prob1 is None:
        assert prob0 is None
        prob1 = np.array(infected)
        prob0 = np.array(1 - infected)
    
      num_unconfident_negatives = 0
      if prefer_recall:
        # Report the unconfident -ves as +ve
        negatives = (infected == 0).astype(np.int32)
        unconfident_negatives = negatives * (prob0 < 0.6).astype(np.int32)
        num_unconfident_negatives = np.sum(unconfident_negatives)
        infected = infected + unconfident_negatives

      # Get definite defects
      bool_y = (y > 0).astype(np.int32)
      
      # Finding the sure positives, sure negatives, and potential positives with only the group test outcomes  
      _infected_comp, infected_dd, _score, _tp, _fp, _fn, surep, _unsurep, _ =\
        self.decode_comp_new(bool_y, compute_stats=compute_stats)

      #print(infected.shape)
      #print(infected_dd.shape)
      # Compare definite defects with ours to detect if our algorithm doesn't
      # detect something that should definitely have been detected
      wrongly_undetected = np.sum(infected_dd - infected_dd * infected)

      # Add infections from high precision algo
      infected_high_precision = (answer_high_precision > 0).astype(np.int32)
      # For ease of implementation we add above to infected_dd. This will become
      # sure_list later
      infected_dd = (infected_dd + infected_high_precision > 0).astype(np.int32)
      # Updating infected with the results obtained from only the group test outcomes
      infected = (infected + infected_dd > 0).astype(np.int32)
    
    # Performing P-Best with full pooling matrix so no COMP  
    elif algo == "PBest_full":
      
      # Coverting viral loads to binary
      bool_y = (results != 0.).astype(np.int32)
      # Since model is always fitted, it does not make sense to record ambiguity so NA.
      Ambi = pd.NA
      # binary y
      y = bool_y
      
      # Making inputs compatible for MATLAB functions
      A = np.ascontiguousarray(self.M.T, dtype = np.float64)
      A = matlab.double(A.tolist())
      
      maxNum =  20
      
      # Maximum viral loads 
      dt = np.max(np.abs(np.array(A).T @ np.array(y)))
      
      # Regularization parameter value
      tau = 0.005*dt;
      
      # Length of y
      m = len(y)
      
      # Making inputs compatible
      qy = matlab.double(y.astype(float).tolist())
      qy = eng.reshape(qy, int(m), int(1), nargout = 1)
      
      # Running PBest's MATLAB function opm()
      u = eng.opm(qy, A, tau, maxNum)
      
      # Runing PBest's MATLAB function selectByError()
      discreteOutput = eng.selectByError(u, A, qy)
      
      # Index of specimens detected as positive
      detected_samples = eng.find(discreteOutput, nargout = 1)
      
      # Adjusting the index of positive specimens as Python is 0-index
      if isinstance(detected_samples, (int, float)):
        # Scalar case
        detected_samples = np.array([int(detected_samples)]) - 1
        
      elif hasattr(detected_samples, '_data'):
        
        raw = bytes(detected_samples._data)
        
        n_elements = len(raw) // 8
        
        detected_samples = np.array([struct.unpack('d', raw[i*8:(i+1)*8])[0] 
                                  for i in range(n_elements)])
                                  
        detected_samples = detected_samples.astype(int) - 1
        
      else:
        detected_samples = np.array([int(i) for i in detected_samples]).flatten() - 1
      
      # Places to store infected results
      infected = np.zeros(self.n)
      answer = np.zeros(self.n)
      prob1_new = np.zeros(self.n)
      prob0_new = np.ones(self.n)
      determined = 1
      overdetermined = 0
      num_unconfident_negatives = 0
      
      # Regularization parameter 
      lamda = tau
      # Number of retest
      retest = 0
      # Updating infected based on the index of detected samples
      infected[detected_samples] = 1
      
      score = np.linalg.norm(answer - self.conc) / math.sqrt(self.t)

      # Finding the sure positives, sure negatives, and potential positives with only the group test outcomes  
      _infected_comp, infected_dd, _score, _tp, _fp, _fn, surep, _unsurep, _ =\
        self.decode_comp_new(bool_y, compute_stats=compute_stats)

      # Compare definite defects with ours to detect if our algorithm doesn't
      # detect something that should definitely have been detected
      wrongly_undetected = np.sum(infected_dd - infected_dd * infected)

      # Add infections from high precision algo
      infected_high_precision = (answer_high_precision > 0).astype(np.int32)
      # For ease of implementation we add above to infected_dd. This will become
      # sure_list later
      infected_dd = (infected_dd + infected_high_precision > 0).astype(np.int32)
    
    elif algo == "PBest-2024":
      
      # Direct translation of the P-BEST authors' MATLAB code
      qMeasurement = results
      
      # Making inputs compatible for MATLAB functions
      A = np.ascontiguousarray(self.M.T, dtype = np.float64)
      A = matlab.double(A.tolist())
      
      # Maximum viral loads 
      dt = np.max(np.abs(np.array(A).T @ np.array(qMeasurement)))
 
      # Regularization parameter
      tau = 0.005 * dt
      
      # Value of regularization parameter
      lamda = tau
      
      # Number of groups
      m = len(qMeasurement)
 
      # Make group viral loads compatible with the MATLAB engine
      qy = matlab.double(qMeasurement.tolist())
      qy = eng.reshape(qy, float(m), 1.0, nargout=1)   # (m x 1) column vector
 
      # Fit GPSR to estimate individual viral loads
      estimate_VL = eng.applyGPSR(qy, A, tau)
 
      # Convert the MATLAB double back to a NumPy array
      estimate_VL_r = np.asarray(estimate_VL).flatten()
 
      # Indices of potential positive specimens
      index_pot_pos = np.where(estimate_VL_r > 0)[0]
      index_pot_pos = np.asarray(index_pot_pos, dtype=int)
      
      # Indices of negative specimens
      index_spe_neg = np.setdiff1d(np.arange(self.n), index_pot_pos)
      index_spe_neg = np.asarray(index_spe_neg, dtype=int)
      
      # Cut off for large viral loads
      cut_off_VL_large = 0.2
      
      # 93 times 961 pooling matrix
      Num_Group_Each_Sample = 3
      
      # If there is at least one potential positive
      # (R's original `!is.null(index_pot_pos)` is always TRUE because which()
      #  returns integer(0), not NULL, when empty; using a length check matches
      #  the clear intent and the else-branch below.)
      if len(index_pot_pos) > 0:
 
        # Reduce M to the columns of potential positive specimens
        M_new = self.M.T[:, index_pot_pos]
        # --- Criterion 1: There is a group where this specimen is the only potential positive
        # Row sums of the reduced pooling matrix
        row_sum_M_new = M_new.sum(axis = 1)
        # Rows (groups) containing exactly one potential-positive specimen
        sum1_row = np.where(row_sum_M_new == 1)[0]
 
        # Place to store index of specimens satisfying criterion 1
        index_crit1 = []
        for rs in sum1_row:
          # index of specimen alone in the group
          lone_spe = np.where(M_new[rs, :] == 1)[0]
          # add the index to the list of index satisfying criterion 1
          index_crit1.extend(lone_spe.tolist())
        
        # Unique index of specimens
        index_crit1 = np.unique(index_crit1).astype(int)
 
        # True specimen indices (in M) satisfying criterion 1
        true_index_crit1 = index_pot_pos[index_crit1]
        
        # Specimens not satisfying criterion 1 are immediately added to the list of suspecteds
        suspected = list(np.setdiff1d(index_pot_pos, true_index_crit1))
 
        # Criterion 2: every group the specimen is in has a large viral loads
        if len(index_crit1) > 0:
          # place to store index of predicted positives
          pred_pos = []
          for pc in true_index_crit1:
            # Groups this specimen belongs to
            spe_in_gr = np.where(self.M.T[:, pc] == 1)[0]
            # Viral loads of those groups
            VL_spe_in_gr = qMeasurement[spe_in_gr]
            # How many of those groups have a large viral load
            num_large_VL_gr = np.sum(VL_spe_in_gr >= cut_off_VL_large)
            
            if num_large_VL_gr < Num_Group_Each_Sample:
                # add index to the list of suspecteds
                suspected.append(pc)          # fails criterion 2
                
            else:
                # add index to the list of predicted positives
                pred_pos.append(pc)            # predicted positive
        else:
          # Index of predicted positives
          pred_pos = []
          # suspected stays as is
 
        # Number of specimens satisfying criteria 1 and 2
        len_pred_pos = len(pred_pos)
 
        # Place to store test results
        infected = np.full(self.n, np.nan)
        
        if len_pred_pos > 0:
          # test results
          infected[np.array(pred_pos, dtype=int)] = 1
        # Assigning results of negative specimens
        infected[index_spe_neg] = 0
 
        # Retest the suspected specimens
        suspected = np.array(suspected, dtype=int)
        # number of retest
        retest = len(suspected)
        # true VL of suspected specimens
        true_VL_suspected = self.conc[suspected]
        # Assigning test results of binary specimens
        infected[suspected] = (true_VL_suspected > 0).astype(np.int32)
 
      else:
        # test results
        infected = np.zeros(Num_Sample)
        # Number of retest
        retest = 0
      
      # Ambiguity
      Ambi = pd.NA
      
      # All the quantities below this line are needed for completion
      answer = np.zeros(self.n)
      infected_internal = np.zeros(self.n)
      prob1_new = np.zeros(self.n)
      prob0_new = np.ones(self.n)
      determined = 1
      overdetermined = 0
      
      # All the quantities from this line are not relevant for us 
      prob1 = infected_internal.astype(np.float32)
      prob0 = (1 - infected_internal).astype(np.float32)

      score = np.linalg.norm(answer - self.conc) / math.sqrt(self.t)
      
      if prob1 is None:
        assert prob0 is None
        prob1 = np.array(infected)
        prob0 = np.array(1 - infected)
    
      num_unconfident_negatives = 0

      # Compare definite defects with ours to detect if our algorithm doesn't
      # detect something that should definitely have been detected
      wrongly_undetected = 0
      infected_dd = np.zeros(self.n)
      
    # PBest with reduced pooling matrix so there will be a COMP step  
    elif algo == "PBest_reduced":
      
      # Converting viral loads to binary
      bool_y = (results != 0.).astype(np.int32)
      # Finding the sure positives, sure negatives, and potential positives with only the group test outcomes  
      infected_comp, infected_dd, _score, _tp, _fp, _fn, surep, unsurep,\
        num_infected_in_test = self.decode_comp_new(bool_y, compute_stats)
      # Checking ambiguity
      Ambi = self.AmbiCheck(results)
      
      # Binary y
      y = bool_y
      
      # Getting the index of potential positive specimens that could not be determined with only group test outcomes
      non_zero_cols,  = np.nonzero(infected_comp)
      # Positive groups
      non_zero_rows,  = np.nonzero(y)
      
      # Subseting the pooling matrix with positive groups and index of potential positives
      A = self.M.T
      A = np.take(A, non_zero_cols, axis=1)
      A = np.take(A, non_zero_rows, axis=0)
      # Viral loads of positive groups
      y = y[non_zero_rows]
    
      maxNum =  20
      
      # Maximum viral loads
      dt = np.max(np.abs(np.array(A).T @ np.array(y)))
      
      # Regularization parameter value
      tau = 0.005*dt
      
      # Length of y
      m = len(y)
      
      # Making the inputs compatible for MATLAB function
      qy = matlab.double(y.astype(float).tolist())
      qy = eng.reshape(qy, int(m), int(1), nargout = 1)
      
      # Running PBest's MATLAB function opm()
      u = eng.opm(qy, A, tau, maxNum)
      
      # Running PBest's MATLAB function selectByError()
      discreteOutput = eng.selectByError(u, A, qy)
      
      # Index of specimens declared as positive
      detected_samples = eng.find(discreteOutput, nargout = 1)
      
      # Adjusting the index of specimens as Python is 0-index
      if isinstance(detected_samples, (int, float)):
        # Scalar case
        detected_samples = np.array([int(detected_samples)]) - 1
        
      elif hasattr(detected_samples, '_data'):
        raw = bytes(detected_samples._data)
        # Each double is 8 bytes
        n_elements = len(raw) // 8
        detected_samples = np.array([struct.unpack('d', raw[i*8:(i+1)*8])[0] 
                                  for i in range(n_elements)])
        detected_samples = detected_samples.astype(int) - 1  # 0-based indexing
      else:
        detected_samples = np.array([int(i) for i in detected_samples]).flatten() - 1
      
      # Place to store infected results
      infected = np.zeros(self.n)
      answer = np.zeros(self.n)
      prob1_new = np.zeros(self.n)
      prob0_new = np.ones(self.n)
      determined = 1
      overdetermined = 0
      # Number of retest
      retest = 0
      # Value of regularization parameter
      lamda = tau
      
      # Updating infected result based on the index of potential positives
      infected[non_zero_cols[detected_samples]] = 1
      
      score = np.linalg.norm(answer - self.conc) / math.sqrt(self.t)

      if prob1 is None:
        assert prob0 is None
        prob1 = np.array(infected)
        prob0 = np.array(1 - infected)
    
      num_unconfident_negatives = 0
      if prefer_recall:
        # Report the unconfident -ves as +ve
        negatives = (infected == 0).astype(np.int32)
        unconfident_negatives = negatives * (prob0 < 0.6).astype(np.int32)
        num_unconfident_negatives = np.sum(unconfident_negatives)
        infected = infected + unconfident_negatives

      # Compare definite defects with ours to detect if our algorithm doesn't
      # detect something that should definitely have been detected
      wrongly_undetected = np.sum(infected_dd - infected_dd * infected)

      # Add infections from high precision algo
      infected_high_precision = (answer_high_precision > 0).astype(np.int32)
      # For ease of implementation we add above to infected_dd. This will become
      # sure_list later
      infected_dd = (infected_dd + infected_high_precision > 0).astype(np.int32)
      
      # Updating infected with the results obtained from only the group test outcomes
      infected = (infected + infected_dd > 0).astype(np.int32)
    
    # Dorfman testing with optimal group size at Sensitivity = 1, Specificity = 1
    elif algo == "Opt_Dorfman_at_Se1":
      # Converting viral loads to binary
      bool_y = (results != 0.).astype(np.int32)
      # Finding the sure positives, sure negatives, and potential positives with only the group test outcomes  
      infected_comp, infected_dd, _score, _tp, _fp, _fn, surep, unsurep,\
        num_infected_in_test = self.decode_comp_new_Se1(bool_y, compute_stats)
      
      # Viral loads response
      y = results
      # Getting the index of potential positive specimens that could not be determined with only group test outcomes
      non_zero_cols,  = np.nonzero(infected_comp)
      # Positive groups
      non_zero_rows,  = np.nonzero(y)
      
      # Subseting the pooling matrix with positive groups and index of potential positives
      A = self.DT_opt_Se1
      A = np.take(A, non_zero_cols, axis=1)
      A = np.take(A, non_zero_rows, axis=0)
      
      # Viral loads of positive groups
      y = y[non_zero_rows]
      
      # Viral loads of specimens 
      x = self.conc
      # Viral loads of potential positive specimens
      x = x[non_zero_cols]
      # Converting viral loads of specimens to binary
      arr = (x>0).astype(np.int32)
      
      # Now solve using this new A and y
      n = A.shape[1]
      t = A.shape[0]
      # Parameters d and s do not matter. They are not used in the algorithm
      d = self.d
      s = self.s

      if t == 0:
        assert n == 0
      
      # Ambiguity is not applicable with Dorfman testing
      Ambi = pd.NA
      
      # Places to store results
      infected = np.zeros(self.n)
      answer = np.zeros(self.n)
      prob1_new = np.zeros(self.n)
      prob0_new = np.ones(self.n)
      determined = 1
      overdetermined = 0
      
      # Regularization parameter value
      lamda = None
      
      # Error in individual viral loads upon retest
      error_new = np.random.normal(0., config.eps_std_dev, size = n)
      
      # Viral loads upon retest  
      x_retest = x * ((1 + config.p) ** error_new)
      
      # Converting viral loads upon retest to binary  
      infected_internal = (x_retest > 0.).astype(np.int32)
      answer_internal = x_retest
      
      # Number of retest  
      retest = len(non_zero_cols)
      
      # prob1 and prob0 are not relevant   
      prob1 = infected_internal.astype(np.float32)
      prob0 = (1 - infected_internal).astype(np.float32)
      determined = 1 if t == n else 0
      overdetermined = 1 if t > n else 0
      
      # Updating infected based on retest
      for ans, val, idx in zip(answer_internal, infected_internal, non_zero_cols):
        infected[idx] = val
        answer[idx] = ans
      
      # prob1_new, prob0_new are not important for us but need to compute to run the code without error
      for p1, p0, idx in zip(prob1, prob0, non_zero_cols):
        prob1_new[idx] = p1
        prob0_new[idx] = p0
      score = np.linalg.norm(answer - self.conc) / math.sqrt(self.t)
      
      # prob1 and prob0 are also not relevant 
      if prob1 is None:
        assert prob0 is None
        prob1 = np.array(infected)
        prob0 = np.array(1 - infected)
    
      num_unconfident_negatives = 0

      # Compare definite defects with ours to detect if our algorithm doesn't
      # detect something that should definitely have been detected
      wrongly_undetected = np.sum(infected_dd - infected_dd * infected)
      # Add infections from high precision algo
      infected_high_precision = (answer_high_precision > 0).astype(np.int32)
    
    # Dorfman testing with optimal group size at Sensitivity = 0.95, Specicificity = 1
    elif algo == "Opt_Dorfman_at_Se95by100":
      
      # Converting viral loads to binary
      bool_y = (results != 0.).astype(np.int32)
      
      # Finding the sure positives, sure negatives, and potential positives with only the group test outcomes  
      infected_comp, infected_dd, _score, _tp, _fp, _fn, surep, unsurep,\
        num_infected_in_test = self.decode_comp_new_Se95by100(bool_y, compute_stats)
      # Viral loads response
      y = results
      # Getting the index of potential positive specimens 
      non_zero_cols,  = np.nonzero(infected_comp)
      # Index of positive groups
      non_zero_rows,  = np.nonzero(y)
      
      # Subseting A with only taking the positive groups and index of potential positive specimens 
      A = self.DT_opt_Se95by100
      A = np.take(A, non_zero_cols, axis=1)
      A = np.take(A, non_zero_rows, axis=0)
      
      # Viral loads of positive groups
      y = y[non_zero_rows]
      
      # Viral loads of all the specimens 
      x = self.conc
      # Viral loads of potential postives
      x = x[non_zero_cols]
      # Converting viral loads of potential positives to binary
      arr = (x>0).astype(np.int32)

      # Dimensions of the reduced pooling matrix 
      n = A.shape[1]
      t = A.shape[0]
      # Parameters d and s do not matter. They are not used in the algorithm
      d = self.d
      s = self.s

      if t == 0:
        assert n == 0
      
      # Ambiguity is NA with Dorfman testing
      Ambi = pd.NA
      
      # Places to store infected results
      infected = np.zeros(self.n)
      answer = np.zeros(self.n)
      prob1_new = np.zeros(self.n)
      prob0_new = np.ones(self.n)
      determined = 1
      overdetermined = 0
      # Value of regularization parameter
      lamda = None
      
      # Error in individual viral loads upon retest
      error_new = np.random.normal(0., config.eps_std_dev, size = n)
      # Viral loads of potential positives upon retest
      x_retest = x * ((1 + config.p) ** error_new)
      # Converting individual viral loads upon retest to binary  
      infected_internal = (x_retest > 0.).astype(np.int32)
      answer_internal = x_retest
      
      # Number of retest  
      retest = len(non_zero_cols)
      
      # Assigning results to infected
      for ans, val, idx in zip(answer_internal, infected_internal, non_zero_cols):
        infected[idx] = val
        answer[idx] = ans
      
      # All the quantities from this line are not relevant for us  
      prob1 = infected_internal.astype(np.float32)
      prob0 = (1 - infected_internal).astype(np.float32)
      determined = 1 if t == n else 0
      overdetermined = 1 if t > n else 0
      
      for p1, p0, idx in zip(prob1, prob0, non_zero_cols):
        prob1_new[idx] = p1
        prob0_new[idx] = p0
      
      score = np.linalg.norm(answer - self.conc) / math.sqrt(self.t)
      
      if prob1 is None:
        assert prob0 is None
        prob1 = np.array(infected)
        prob0 = np.array(1 - infected)
    
      num_unconfident_negatives = 0

      # Compare definite defects with ours to detect if our algorithm doesn't
      # detect something that should definitely have been detected
      wrongly_undetected = np.sum(infected_dd - infected_dd * infected)

      # Add infections from high precision algo
      infected_high_precision = (answer_high_precision > 0).astype(np.int32)
      # For ease of implementation we add above to infected_dd. This will become
      # sure_list later
      infected_dd = (infected_dd + infected_high_precision > 0).astype(np.int32)
      
    elif algo == 'Dorfman':
      
      # Binary response
      bool_y = (results != 0.).astype(np.int32)
      
      # Finding the sure positives, sure negatives, and potential positives with only the group test outcomes  
      infected_comp, infected_dd, _score, _tp, _fp, _fn, surep, unsurep,\
        num_infected_in_test = self.decode_comp_new3(bool_y, compute_stats)
      
      # Viral loads response
      y = results
      
      # Index of potential positive specimens
      non_zero_cols,  = np.nonzero(infected_comp)
      # Index of positive groups
      non_zero_rows,  = np.nonzero(y)
      
      # Pooling matrix
      A = self.DT
      # Subsetting pooling matrix taking only the rows of positive groups and columns of potential positives
      A = np.take(A, non_zero_cols, axis=1)
      A = np.take(A, non_zero_rows, axis=0)
      # Viral loads of positive groups
      y = y[non_zero_rows]
      
      # Individual viral loads
      x = self.conc
      # Viral loads of potential positives
      x = x[non_zero_cols]
      # Converting viral loads of potential positives to binary
      arr = (x>0).astype(np.int32)
      # Dimensions of reduced pooling matrix
      n = A.shape[1]
      t = A.shape[0]
      # parameters d and s do not matter. They are not used in the algorithm
      d = self.d
      s = self.s

      if t == 0:
        assert n == 0
      # Ambiguity is NA for Dorfman
      Ambi = pd.NA
      
      # Places to store results
      infected = np.zeros(self.n)
      answer = np.zeros(self.n)
      prob1_new = np.zeros(self.n)
      prob0_new = np.ones(self.n)
      determined = 1
      overdetermined = 0
      lamda = None
      
      # Error in individual viral loads  
      error_new = np.random.normal(0., config.eps_std_dev, size = n)
      # Individual viral loads upon retest   
      x_retest = x * ((1 + config.p) ** error_new)
      # Binarizing viral loads upon retest using cutoff 0  
      infected_internal = (x_retest > 0.).astype(np.int32)
      answer_internal = x_retest
      # Number of retest  
      retest = len(non_zero_cols)
      
      # Updating infected based on the results obtained with retest  
      for ans, val, idx in zip(answer_internal, infected_internal, non_zero_cols):
        infected[idx] = val
        answer[idx] = ans
      
      # All the quantities from this line are not relevant for us 
      prob1 = infected_internal.astype(np.float32)
      prob0 = (1 - infected_internal).astype(np.float32)
        
      determined = 1 if t == n else 0
      overdetermined = 1 if t > n else 0
      
      for p1, p0, idx in zip(prob1, prob0, non_zero_cols):
        prob1_new[idx] = p1
        prob0_new[idx] = p0
      
      score = np.linalg.norm(answer - self.conc) / math.sqrt(self.t)
      
      if prob1 is None:
        assert prob0 is None
        prob1 = np.array(infected)
        prob0 = np.array(1 - infected)
    
      num_unconfident_negatives = 0

      # Compare definite defects with ours to detect if our algorithm doesn't
      # detect something that should definitely have been detected
      wrongly_undetected = np.sum(infected_dd - infected_dd * infected)

      # Add infections from high precision algo
      infected_high_precision = (answer_high_precision > 0).astype(np.int32)
      # For ease of implementation we add above to infected_dd. This will become
      # sure_list later
      infected_dd = (infected_dd + infected_high_precision > 0).astype(np.int32)

    # "regular_ORUA": regular testing with only retest upon ambiguity     
    elif algo == "regular_ORUA":
      
      # Viral loads response
      y = results
      # Checking ambiguity
      Ambi = self.AmbiCheck(y)
      # Binary response
      bool_y = (y > 0).astype(np.int32)
      # Finding the sure positives, sure negatives, and potential with only the group test outcomes  
      infected_comp, infected_dd, _score, _tp, _fp, _fn, surep, unsurep, \
      num_infected_in_test = self.decode_comp_new(bool_y, compute_stats)
      
      # Index of potential positive specimens that could not be identified with COMP
      non_zero_cols,  = np.nonzero(infected_comp)
      # Index of positive groups
      non_zero_rows,  = np.nonzero(y)
      # Individual viral loads of all the specimens
      x = self.conc
      # Individual viral loads of specimens
      x = x[non_zero_cols]
      # Converting individual viral loads to binary
      arr = (x>0).astype(np.int32)
      
      # Places to store results
      infected = np.zeros(self.n)
      answer = np.zeros(self.n)
      prob1_new = np.zeros(self.n)
      prob0_new = np.ones(self.n)
      determined = 1
      overdetermined = 0
      
      # If there is ambiguity retest
      if Ambi == 1:
        # Number of retest
        retest = len(non_zero_cols)
        n = retest
        # Number of positive groups
        t = len(non_zero_rows)
        # Error in individual viral loads
        error_new = np.random.normal(0., config.eps_std_dev, size = n)
        # Individual viral loads upon retest
        x_retest = x * ((1 + config.p) ** error_new)
        # Converting viral loads upon retest to binary using cutoff of 0
        infected_internal = (x_retest > 0.).astype(np.int32)
        answer_internal = x_retest
        # Updating infected using the results obtained with retest
        for ans, val, idx in zip(answer_internal, infected_internal, non_zero_cols):
          infected[idx] = val
          answer[idx] = ans
        
        # Quantities needed to run code without any errors
        # These are not relevant to us
        prob1 = infected_internal.astype(np.float32)
        prob0 = (1 - infected_internal).astype(np.float32)
        score = 0.0
        
        determined = 1 if t == n else 0
        overdetermined = 1 if t > n else 0
        for p1, p0, idx in zip(prob1, prob0, non_zero_cols):
          prob1_new[idx] = p1
          prob0_new[idx] = p0
      # If there is no ambiguity no retest, automatically declare potential positives as positive and others as negative    
      else:
        # Number of retest
        retest = 0
        
        n = retest
        # Number of positive groups
        t = len(non_zero_rows)
        
        infected_internal = np.ones(len(non_zero_cols))
        answer_internal = x
        # Updating infected using the results obtained with only testing groups
        for ans, val, idx in zip(answer_internal, infected_internal, non_zero_cols):
          infected[idx] = val
          answer[idx] = ans
        # prob1, prob0, score, determined, overdetermined, prob1_new, prob0_new are not relevant for us
        prob1 = infected_internal.astype(np.float32)
        prob0 = (1 - infected_internal).astype(np.float32)
        score = 0.0
        determined = 1 if t == n else 0
        overdetermined = 1 if t > n else 0
        for p1, p0, idx in zip(prob1, prob0, non_zero_cols):
          prob1_new[idx] = p1
          prob0_new[idx] = p0

      if prob1 is None:
        assert prob0 is None
        prob1 = np.array(infected)
        prob0 = np.array(1 - infected)
      
      # num_unconfident_negatives is not relavant
      num_unconfident_negatives = 0
      
      # No regularization parameter neededed for regular testing
      lamda = None
      
      if prefer_recall:
        # Report the unconfident -ves as +ve
        negatives = (infected == 0).astype(np.int32)
        unconfident_negatives = negatives * (prob0 < 0.6).astype(np.int32)
        num_unconfident_negatives = np.sum(unconfident_negatives)
        infected = infected + unconfident_negatives

      # Compare definite defects with ours to detect if our algorithm doesn't
      # detect something that should definitely have been detected
      wrongly_undetected = np.sum(infected_dd - infected_dd * infected)

      # Add infections from high precision algo
      infected_high_precision = (answer_high_precision > 0).astype(np.int32)
      # For ease of implementation we add above to infected_dd. This will become sure_list later
      infected_dd = (infected_dd + infected_high_precision > 0).astype(np.int32)
      # Updating infected with the results obtained from only the group test outcomes
      infected = (infected + infected_dd > 0).astype(np.int32)
    
    # Regular testing with always retest
    elif algo == "regular_AR":
      # viral loads response
      y = results
      # Ambiguity is NA for regular testing
      Ambi = pd.NA
      # Binary response
      bool_y = (y > 0).astype(np.int32)
      # Finding the sure positives, sure negatives, and potential positives with only the group test outcomes  
      infected_comp, infected_dd, _score, _tp, _fp, _fn, surep, unsurep, \
      num_infected_in_test = self.decode_comp_new(bool_y, compute_stats)
      # index of potential positive specimens 
      non_zero_cols,  = np.nonzero(infected_comp)
      # Index of positive groups
      non_zero_rows,  = np.nonzero(y)
      # Viral loads of each specimen
      x = self.conc
      # Viral loads of potential positives
      x = x[non_zero_cols]
      # Converting viral loads to binary
      arr = (x>0).astype(np.int32)
      
      # Places to store results
      infected = np.zeros(self.n)
      answer = np.zeros(self.n)
      prob1_new = np.zeros(self.n)
      prob0_new = np.ones(self.n)
      determined = 1
      overdetermined = 0
      # No regularization parameter with regular testing
      lamda = None
      # Number of retest
      retest = len(non_zero_cols)
      n = retest
      # Number of postive groups
      t = len(non_zero_rows)
      # Error in viral loads upon retest  
      error_new = np.random.normal(0., config.eps_std_dev, size = n)
      # Viral loads upon retest  
      x_retest = x * ((1 + config.p) ** error_new)
      
      # Coverting viral loads upon retest to binary  
      infected_internal = (x_retest > 0.).astype(np.int32)
      answer_internal = x_retest
        
      # Updating infected based on the resuilts obtained upon retest
      for ans, val, idx in zip(answer_internal, infected_internal, non_zero_cols):
        infected[idx] = val
        answer[idx] = ans
      # prob1, prob0, prob1_new, prob0_new, score, determined, num_unconfident_negatives, overdetermined are not relevant in our settings; just need them to run the code without errors
      prob1 = infected_internal.astype(np.float32)
      prob0 = (1 - infected_internal).astype(np.float32)
        
      determined = 1 if t == n else 0
      overdetermined = 1 if t > n else 0
      
      for p1, p0, idx in zip(prob1, prob0, non_zero_cols):
        prob1_new[idx] = p1
        prob0_new[idx] = p0
      
      score = np.linalg.norm(answer - self.conc) / math.sqrt(self.t)

      if prob1 is None:
        assert prob0 is None
        prob1 = np.array(infected)
        prob0 = np.array(1 - infected)
    
      num_unconfident_negatives = 0
      if prefer_recall:
        # Report the unconfident -ves as +ve
        negatives = (infected == 0).astype(np.int32)
        unconfident_negatives = negatives * (prob0 < 0.6).astype(np.int32)
        num_unconfident_negatives = np.sum(unconfident_negatives)
        infected = infected + unconfident_negatives

      # Compare definite defects with ours to detect if our algorithm doesn't
      # detect something that should definitely have been detected
      wrongly_undetected = np.sum(infected_dd - infected_dd * infected)

      # Add infections from high precision algo
      infected_high_precision = (answer_high_precision > 0).astype(np.int32)
      # For ease of implementation we add above to infected_dd. This will become
      # sure_list later
      infected_dd = (infected_dd + infected_high_precision > 0).astype(np.int32)
      # Updating infected with the results obtained from only the group test outcomes
      infected = (infected + infected_dd > 0).astype(np.int32)
    
    # Did not make any changes in algo == 'OMP'   
    elif algo == 'OMP':
      temp_mat=(self.M.T).astype(float)
      temp_mat=pylops.MatrixMult(temp_mat)
      answer = pylops.optimization.sparsity.OMP(temp_mat, results, 10000,
          sigma=0.001)[0]
          
    elif algo== 'NNOMP':
      # Max d that can be detected by NNOMP is equal to number of rows
      answer=nnompcv.nnomp(self.M.T.astype('float'), 0, results, 0, self.t, cv=False)
      infected = (answer != 0.).astype(np.int32)
      # Number of retest
      retest = 0
      # Viral loads response
      y = results
      # Checking ambiguity
      Ambi = self.AmbiCheck(y)
      # No regularization parameter value 
      lamda = None
      
      score = np.linalg.norm(answer - self.conc) / math.sqrt(self.t)

      if prob1 is None:
        assert prob0 is None
        prob1 = np.array(infected)
        prob0 = np.array(1 - infected)
    
      num_unconfident_negatives = 0
      
      
      if prefer_recall:
        # Report the unconfident -ves as +ve
        negatives = (infected == 0).astype(np.int32)
        unconfident_negatives = negatives * (prob0 < 0.6).astype(np.int32)
        num_unconfident_negatives = np.sum(unconfident_negatives)
        infected = infected + unconfident_negatives
      # Get definite defects
      bool_y = (y > 0).astype(np.int32)
      # Finding the sure positives, sure negatives, and potential positives with only the group test outcomes   
      _infected_comp, infected_dd, _score, _tp, _fp, _fn, surep, _unsurep, _ =\
        self.decode_comp_new(bool_y, compute_stats=compute_stats)

      # Compare definite defects with ours to detect if our algorithm doesn't
      # detect something that should definitely have been detected
      wrongly_undetected = np.sum(infected_dd - infected_dd * infected)

      # Add infections from high precision algo
      infected_high_precision = (answer_high_precision > 0).astype(np.int32)
      # For ease of implementation we add above to infected_dd. This will become
      # sure_list later
      infected_dd = (infected_dd + infected_high_precision > 0).astype(np.int32)
      # Updating infected with the results obtained from only the group test outcomes 
      infected = (infected + infected_dd > 0).astype(np.int32)
    # Made no changes in NNOMPCV  
    elif algo=='NNOMPCV':
      temp_mat = (self.M.T).astype(float)
      mr = math.ceil(0.9*temp_mat.shape[1])
      m = temp_mat.shape[1]
      Ar = temp_mat[0:mr, :]
      Acv = temp_mat[mr+1:m, :]
      yr = results[0:mr]
      ycv = results[mr+1:m]
      #print('yo')
      # Max d that can be detected by NNOMP is equal to number of rows
      answer = nnompcv.nnomp(Ar, Acv, yr, ycv, self.t, cv=True)
    
    # No changes in NNOMP_loo_cv
    elif algo == 'NNOMP_loo_cv':
      answer, prob1, prob0 = self.decode_nnomp_multi_split_cv(results, 'loo_splits')
    
    # No changes in NNOMP_random_cv
    elif algo == 'NNOMP_random_cv':
      # Skip cross-validation for really small cases
      if np.sum(results) == 0:
        answer = np.zeros(self.n)
      elif self.t < 4:
        # Max d that can be detected by NNOMP is equal to number of rows
        answer = nnompcv.nnomp(self.M.T.astype('float'),0,results,0, self.t, cv=False)
      else:
        answer, prob1, prob0 = self.decode_nnomp_multi_split_cv(results, 'random_splits')
      
      infected = (answer != 0.).astype(np.int32)
      
    elif algo.startswith('combined_COMP_'):
    
      l = len('combined_COMP_')
      # Name of secondary algorithm
      secondary_algo = algo[l:]
      # Running decode_comp_combined() with the secondary algorithm
      answer, infected, prob1, prob0, determined, overdetermined, retest, Ambi, lamda =\
        self.decode_comp_combined(results, secondary_algo,
        compute_stats=compute_stats)
      # infected results
      infected = (answer != 0.).astype(np.int32)
      
      # prob1, prob0, prob1_new, prob0_new, score, determined, num_unconfident_negatives, overdetermined are not relevant in our settings; just need them to run the code without errors
      score = np.linalg.norm(answer - self.conc) / math.sqrt(self.t)
      
      if prob1 is None:
        assert prob0 is None
        prob1 = np.array(infected)
        prob0 = np.array(1 - infected)
    
      num_unconfident_negatives = 0
      
      if prefer_recall:
        # Report the unconfident -ves as +ve
        negatives = (infected == 0).astype(np.int32)
        unconfident_negatives = negatives * (prob0 < 0.6).astype(np.int32)
        num_unconfident_negatives = np.sum(unconfident_negatives)
        infected = infected + unconfident_negatives

      # Viral loads response
      y = results
      # Converting viral loads to binary
      bool_y = (y > 0).astype(np.int32)
      
      # Finding the sure positives, sure negatives, and potential positives with only the group test outcomes
      _infected_comp, infected_dd, _score, _tp, _fp, _fn, surep, _unsurep, _ =\
        self.decode_comp_new(bool_y, compute_stats=compute_stats)

      # Compare definite defects with ours to detect if our algorithm doesn't
      # detect something that should definitely have been detected
      wrongly_undetected = np.sum(infected_dd - infected_dd * infected)

      # Add infections from high precision algo
      infected_high_precision = (answer_high_precision > 0).astype(np.int32)
      # For ease of implementation we add above to infected_dd. This will become
      # sure_list later
      infected_dd = (infected_dd + infected_high_precision > 0).astype(np.int32)
      # Updating infected with the results obtained from only the group test outcomes
      infected = (infected + infected_dd > 0).astype(np.int32)
      
    # No changes in this algorithm
    elif algo.startswith('precise_SBL_'):
      # e.g. combined_SBL_clustered_combined_COMP_SBL
      # e.g. combined_SBL_clustered_COMP
      l = len('precise_SBL_')
      primary_algo = 'combined_COMP_SBL_clustered'
      secondary_algo = algo[l:]
      assert secondary_algo not in ['SBL_clustered', 'combined_COMP_SBL_clustered']

      y = results
      # First run SBL_clustered to get precise results
      # Then run secondary algorithm to get high recall results
      # "answer" is those from high recall ones. 
      # we'll create "answer_precise_SBL" and "infected_precise_SBL". 
      # infected_dd will become union of infected_dd and infected_precise_SBL
      # Hence surep will contain results from SBL_clustered as well.
      answer_high_precision = self.get_high_precision_algo_answer(primary_algo, y)
      answer = self.get_high_recall_algo_answer(secondary_algo, y)
    # Sparse bayesian learning  
    elif algo == 'SBL':
      # Pooling matrix
      A = self.M.T
      # Viral loads response
      y = results
      # Fitting SBL
      answer = sbl.sbl(A, y)
      # Converting estimates to binary
      infected = (answer != 0.).astype(np.int32)
      # Number of retest
      retest = 0
      # Ambiguity
      Ambi = self.AmbiCheck(y)
      # No regularization parameter
      lamda = None
      
      # prob1, prob0, prob1_new, prob0_new, score, determined, num_unconfident_negatives, overdetermined are not relevant in our settings; just need them to run the code without errors
      score = np.linalg.norm(answer - self.conc) / math.sqrt(self.t)
      if prob1 is None:
        assert prob0 is None
        prob1 = np.array(infected)
        prob0 = np.array(1 - infected)
    
      num_unconfident_negatives = 0
      if prefer_recall:
        # Report the unconfident -ves as +ve
        negatives = (infected == 0).astype(np.int32)
        unconfident_negatives = negatives * (prob0 < 0.6).astype(np.int32)
        num_unconfident_negatives = np.sum(unconfident_negatives)
        infected = infected + unconfident_negatives

      # Binary response
      bool_y = (y > 0).astype(np.int32)
      # Finding the sure positives, sure negatives, and potential positives with only the group test outcomes
      _infected_comp, infected_dd, _score, _tp, _fp, _fn, surep, _unsurep, _ =\
        self.decode_comp_new(bool_y, compute_stats=compute_stats)

      # Compare definite defects with ours to detect if our algorithm doesn't
      # detect something that should definitely have been detected
      wrongly_undetected = np.sum(infected_dd - infected_dd * infected)

      # Add infections from high precision algo
      infected_high_precision = (answer_high_precision > 0).astype(np.int32)
      # For ease of implementation we add above to infected_dd. This will become
      # sure_list later
      infected_dd = (infected_dd + infected_high_precision > 0).astype(np.int32)
      # Updating infected with the results obtained from only the group test outcomes
      infected = (infected + infected_dd > 0).astype(np.int32)
    # Least squares  
    elif algo == 'l1ls':
      # Pooling matrix
      A = self.M.T
      # Viral loads
      y = results
      # Fitting l1ls() depending on y
      if np.all(y == 0):
        answer = np.zeros(self.n)
      else:
        answer = l1ls.l1ls(A, y, self.l, self.tau)
      # Converting estimates to binary results 
      infected = (answer != 0.).astype(np.int32)
      # Number of retest
      retest = 0
      # Ambiguity
      Ambi = self.AmbiCheck(y)
      # Regularization parameter value
      lamda = self.l
      # prob1, prob0, prob1_new, prob0_new, score, determined, num_unconfident_negatives, overdetermined are not relevant in our settings; just need them to run the code without errors
      score = np.linalg.norm(answer - self.conc) / math.sqrt(self.t)

      if prob1 is None:
        assert prob0 is None
        prob1 = np.array(infected)
        prob0 = np.array(1 - infected)
    
      num_unconfident_negatives = 0
      
      if prefer_recall:
        # Report the unconfident -ves as +ve
        negatives = (infected == 0).astype(np.int32)
        unconfident_negatives = negatives * (prob0 < 0.6).astype(np.int32)
        num_unconfident_negatives = np.sum(unconfident_negatives)
        infected = infected + unconfident_negatives

      # Coverting viral loads to binary response
      bool_y = (y > 0).astype(np.int32)
      
      # Finding the sure positives, sure negatives, and potential positives with only the group test outcomes
      _infected_comp, infected_dd, _score, _tp, _fp, _fn, surep, _unsurep, _ =\
        self.decode_comp_new(bool_y, compute_stats=compute_stats)

      # Compare definite defects with ours to detect if our algorithm doesn't
      # detect something that should definitely have been detected
      wrongly_undetected = np.sum(infected_dd - infected_dd * infected)

      # Add infections from high precision algo
      infected_high_precision = (answer_high_precision > 0).astype(np.int32)
      # For ease of implementation we add above to infected_dd. This will become
      # sure_list later
      infected_dd = (infected_dd + infected_high_precision > 0).astype(np.int32)
      # Updating infected with the results obtained from only the group test outcomes
      infected = (infected + infected_dd > 0).astype(np.int32)
      
    elif algo == 'l1ls_cv':
      # Pooling matrix
      A = self.M.T
      # Viral loads response
      y = results
      
      # Fitting l1ls()
      sigval = 0.01 * np.mean(y)
      if np.all(y == 0):
        answer = np.zeros(self.n)
      else:
        answer = l1ls.l1ls_cv(A, y, sigval, self.tau)
      
      # Converting estimates to binary results
      infected = (answer != 0.).astype(np.int32)
      # Number of retest
      retest = 0
      # Regularization parameter value
      lamda = self.l
      # prob1, prob0, prob1_new, prob0_new, score, determined, num_unconfident_negatives, overdetermined are not relevant in our settings; just need them to run the code without errors
      score = np.linalg.norm(answer - self.conc) / math.sqrt(self.t)
      if prob1 is None:
        assert prob0 is None
        prob1 = np.array(infected)
        prob0 = np.array(1 - infected)
    
      num_unconfident_negatives = 0
      
      if prefer_recall:
        # Report the unconfident -ves as +ve
        negatives = (infected == 0).astype(np.int32)
        unconfident_negatives = negatives * (prob0 < 0.6).astype(np.int32)
        num_unconfident_negatives = np.sum(unconfident_negatives)
        infected = infected + unconfident_negatives

      # Converting viral loads to binary
      bool_y = (y > 0).astype(np.int32)
      # Finding the sure positives, sure negatives, and potential positives with only the group test outcomes
      _infected_comp, infected_dd, _score, _tp, _fp, _fn, surep, _unsurep, _ =\
        self.decode_comp_new(bool_y, compute_stats=compute_stats)

      # Compare definite defects with ours to detect if our algorithm doesn't
      # detect something that should definitely have been detected
      wrongly_undetected = np.sum(infected_dd - infected_dd * infected)

      # Add infections from high precision algo
      infected_high_precision = (answer_high_precision > 0).astype(np.int32)
      # For ease of implementation we add above to infected_dd. This will become
      # sure_list later
      infected_dd = (infected_dd + infected_high_precision > 0).astype(np.int32)
      # Updating infected with the results obtained from only the group test outcomes
      infected = (infected + infected_dd > 0).astype(np.int32)
    
    # No changes  
    elif algo in algos.algo_dict:
      params = { 'A' : self.M.T, 'y' : results }
      res = algos.algo_dict[algo](params)
      answer = res["x_est"]
      infected = (answer != 0.).astype(np.int32)
    else:
      raise ValueError('No such algorithm %s' % algo)


    if compute_stats:
      # re-compute surep from above infected_dd
      surep = np.sum(infected_dd)

      # Compute stats
      tpos = (infected * self.arr)
      fneg = (1 - infected) * self.arr
      fpos = infected * (1 - self.arr)
      
      tp = sum(tpos)
      fp = sum(fpos)
      fn = sum(fneg)
    
    # PBest with full pooling matrix can predict less correct positive than sure positives
    ## sometimes surrep could be greater than tp hence shoot error. To overcome this error 
    ## commenting it out
      #assert surep <= tp
      unsurep = tp + fp - surep
      
    else:
      tp = 0
      fp = 0
      fn = 0
      surep = 0
      unsurep = 0


    num_infected_in_test = np.zeros(self.t, dtype=np.int32)
    for test in range(self.t):
      for person in range(self.n):
        if infected[person] > 0 and self.M[person, test] == 1:
          num_infected_in_test[test] += 1
    # return outcomes from the algorithms
    return answer, infected, infected_dd, prob1, prob0, score, tp, fp, fn,\
        num_unconfident_negatives, determined, overdetermined, surep,\
        unsurep, wrongly_undetected, num_infected_in_test, retest, Ambi, lamda

  def get_high_precision_algo_answer(self, algo, y):
    assert algo == 'combined_COMP_SBL_clustered'
    x, infected, infected_dd, prob1, prob0, score, tp, fp, fn, uncon_negs, determined,\
        overdetermined, surep, unsurep, wrongly_undetected,\
        num_infected_in_test, lamda = self.decode_lasso(y, algo, prefer_recall=False,
            compute_stats=False)
    # ignore everything and just send x
    return x

  def get_high_recall_algo_answer(self, algo, y):
    if algo == 'COMP':
      bool_y = (y > 0).astype(np.int32)
      infected, infected_dd, score, tp, fp, fn, surep, unsurep,\
          num_infected_in_test = \
          self.decode_comp_new(bool_y, compute_stats=False)
      x = np.zeros(self.n)
      return infected
    else:
      x, infected, infected_dd, prob1, prob0, score, tp, fp, fn, uncon_negs, determined,\
          overdetermined, surep, unsurep, wrongly_undetected,\
          num_infected_in_test = self.decode_lasso(y, algo, prefer_recall=False,
              compute_stats=False)
      return x

  def decode_lasso_for_cv(self, train_Ms, train_ys, test_Ms, test_ys,
      algo='lasso', l=None, sigma=None):

    if algo == 'lasso' and l is None:
      raise ValueError('Need l for algo lasso')
    elif algo == 'OMP' and sigma is None:
      raise ValueError('Need sigma for algo OMP')

    scores = []
    for train_M, train_y, test_M, test_y in zip(train_Ms, train_ys, test_Ms,
        test_ys):
      #print('Doing lasso with')
      #print(train_M.shape, train_y.shape, test_M.shape, test_y.shape)
      if algo == 'lasso':
        #lasso = Lasso(alpha=l, max_iter=10000)
        lasso = LassoCV(n_alphas=100, positive = True, max_iter=10000)
        lasso.fit(train_M, train_y)
        pred_y = lasso.predict(test_M)
      elif algo == 'OMP':
        pass

      score = np.linalg.norm(test_y - pred_y) / len(test_y)
      scores.append(score)

    avg_score = np.average(scores)
    max_score = max(scores)
    median_score = np.median(scores)
    min_score = np.min(scores)
    return avg_score
    #return min_score

  # Get num random splits with given fraction. Sensing matrix will have at
  # most frac fraction of rows
  def return_random_splits(self, y, num, frac, mr=None):
    if mr is None:
      mr = math.floor(frac * self.t)
    else:
      assert mr < self.t
    r = self.t - mr
    # Following code only works for r > 1
    assert r > 1
    
    train_Ms = []
    test_Ms = []
    train_ys = []
    test_ys = []
    M = self.M.T # Uggh
    for i in range(num):
      perm = np.random.permutation(range(self.t))
      r_idx = perm[:r]
      m_idx = perm[r:]

      train_M = np.delete(M, r_idx, axis=0)
      train_y = np.delete(y, r_idx, axis=0)

      test_M = np.delete(M, m_idx, axis=0)
      test_y = np.delete(y, m_idx, axis=0)

      train_Ms.append(train_M)
      train_ys.append(train_y)
      test_Ms.append(test_M)
      test_ys.append(test_y)

    return train_Ms, train_ys, test_Ms, test_ys

  # Return splits for leave-one-out cross-validation
  def return_loo_cv_splits(self, y):
    train_Ms = []
    test_Ms = []
    train_ys = []
    test_ys = []

    # Unfortunately self.M is n x t so we need to transpose it
    M = self.M.T

    # Each row will be left out once as test_M
    for r in range(self.t):
      train_M = np.delete(M, r, axis=0)
      test_M = np.expand_dims(M[r], axis=0)
      train_y = np.delete(y, r, axis=0)
      test_y = np.array([y[r]])

      train_Ms.append(train_M)
      train_ys.append(train_y)
      test_Ms.append(test_M)
      test_ys.append(test_y)

    return train_Ms, train_ys, test_Ms, test_ys

  # Find best d by cross-validation using these splits
  #
  # Best d is the one found by majority of the splits
  def get_d_nnomp_cv(self, splits, max_d, resolve_method='voting', algo='NNOMP'):
    train_Ms, train_ys, test_Ms, test_ys = splits
    counts = np.zeros(max_d + 1)
    cum_error = np.zeros(max_d)
    # Keeps count of number of times each sample was declared as +ve
    x_ones = np.zeros(self.n)
    num_splits = len(train_Ms)
    for train_M, train_y, test_M, test_y in zip(train_Ms, train_ys, test_Ms,
        test_ys):
      x, error, d, errors = nnompcv.nnomp(train_M, test_M, train_y, test_y,
          max_d, cv=True)
      answer = (x > 0).astype(np.int32)
      x_ones += answer
      counts[d] += 1
      #print('Errors: ', np.array(errors))
      if errors:
        cum_error += errors

    best_d_maj = np.argmax(counts) + 1
    best_d_error = np.argmin(cum_error) + 1
    prob_of_one = x_ones / num_splits
    prob_of_zero = 1 - prob_of_one
    #print('prob of one:', prob_of_one)
    #print('prob of zero:', prob_of_zero)
    if resolve_method == 'voting':
      return best_d_maj, prob_of_one, prob_of_zero
    elif resolve_method == 'error':
      return best_d_error, prob_of_one, prob_of_zero
    else:
      raise ValueError('Invalid resolve method %s' % resolve_method)

  # Do leave one out splits
  # get best d from those splits
  # Run final nnomp algorithm using best d and entire matrix
  def decode_nnomp_multi_split_cv(self, y, method='random_splits'):
    if method == 'random_splits':
      splits = self.return_random_splits(y, 100, frac=0.7, mr=self.mr)
    elif method == 'loo_splits':
      splits = self.return_loo_cv_splits(y)

    best_d, prob1, prob0 = self.get_d_nnomp_cv(splits, max_d=self.t)
    if config.prefer_recall:
      best_d = 2 * best_d
    x = nnompcv.nnomp(self.M.T.astype('float'), 0, y, 0,
        best_d, cv=False)
    return x, prob1, prob0
    

  def do_cross_validation_get_lambda(self, y, sigval):
    lambda_min = max([sigval*math.sqrt(math.log(self.n))-5,0.01]);
    lambda_max = sigval*math.sqrt(math.log(self.n))+5;
    n_step = math.ceil((lambda_max - lambda_min) / 0.01)
    ll = np.linspace(lambda_min, lambda_max, n_step)
    #a = np.linspace(0.001, 0.01, num=10)
    #ll = np.concatenate([a, 10*a, 100*a, 1000*a, 10000*a, 100000*a])
    #ll = np.concatenate([a, 10*a, 100*a])
    #ll = np.linspace(0.001, 1., 1000)

    train_Ms = []
    test_Ms = []
    train_ys = []
    test_ys = []
    M = self.M.T
    # We'll do leave one out cross-validation
    for r in range(1):
      train_M = np.delete(M, r, axis=0)
      test_M = np.expand_dims(M[r], axis=0)
      train_y = np.delete(y, r, axis=0)
      test_y = np.array([y[r]])

      train_Ms.append(train_M)
      train_ys.append(train_y)
      test_Ms.append(test_M)
      test_ys.append(test_y)

    scores = []
    for l in ll:
      score = self.decode_lasso_for_cv(train_Ms, train_ys, test_Ms, test_ys,
          l=l)
      scores.append(score)
    scores = np.array(scores)
    idx = np.argmin(scores)
    #print(idx)
    self.l = ll[idx]
    #print('lambdas = ', ll)
    #print('scores = ', scores)
    print('Choosing lambda = %.4f' % self.l, 'score = %.4f' % score)
    return self.l
  
  # Function to check ambiguity
  def AmbiCheck(self, y, compute_stats = True):
    # Coverting viral loads to binary
    bool_y = (y > 0).astype(np.int32)
    # Finding the sure positives, sure negatives, and potential positives with only the group test outcomes
    infected_comp, infected_dd, _score, _tp, _fp, _fn, surep, unsurep, \
        num_infected_in_test = self.decode_comp_new(bool_y, compute_stats)
    # Index of potential positive specimens
    non_zero_cols,  = np.nonzero(infected_comp)
    # Index of positive groups
    non_zero_rows,  = np.nonzero(y)
    # Pooling matrix
    A = self.M.T
    # Subsetting the pooling matrix with taking columns potential positive specimens and index of positive groups
    A = np.take(A, non_zero_cols, axis=1)
    A = np.take(A, non_zero_rows, axis=0)
    # Viral loads of positive groups
    y = y[non_zero_rows]
    
    # Viral loads of all the specimens
    x = self.conc
    # Viral loads of potential positives
    x = x[non_zero_cols]
    # Converting viral loads to binary
    arr = (x>0).astype(np.int32)

    # Dimensions of the reduced pooling matrix
    n = A.shape[1]
    t = A.shape[0]
    
    if t == 0:
      assert n == 0
    # Checking whether there is ambiguity
    if t == 0:
      Ambi = 0
    elif t > 0 and n == 1:
      # Ambi = 0 implies no ambiguity
      Ambi = 0
    else:
      # Row sum
      row_sum = np.sum(A, axis = 1) == 1
      if np.sum(row_sum) < n:
        # Ambi = 1 implies ambiguity
        Ambi = 1
      else:
        # Column sum
        col_sum = np.sum(A, axis = 0)
        Ambi = (np.sum(col_sum == 0) > 0).astype(np.int32)
    # Return ambiguity    
    return Ambi 
  
  
  # Filter out those entries of x which are definitely 0 using COMP.
  # Remove corresponding columns from M.
  def decode_comp_combined(self, y, secondary_algo, test=False,
      compute_stats=True):
    # This assertion is needed because mr depends on number of rows.
    # Since number of rows will change for the internal CS, use frac instead
    assert self.mr == None

    bool_y = (y > 0).astype(np.int32)
    infected_comp, infected_dd, _score, _tp, _fp, _fn, surep, unsurep, \
        num_infected_in_test = self.decode_comp_new(bool_y, compute_stats)

    # Find the indices of 1's above. These will be retained. Rest will be
    # discarded
    #print('Comp output: ', infected_comp)
    non_zero_cols,  = np.nonzero(infected_comp)
    non_zero_rows,  = np.nonzero(y)
    #print('Indices of Non-zero columns:', non_zero_cols)
    #print('Indices of Non-zero rows:', non_zero_rows)

    A = self.M.T

    # Compute errors
    errors = output_validation_utils.detect_discrepancies_in_test(
        A.shape[0], bool_y, num_infected_in_test, log=False)
    
    total_errors = errors['err1'] + errors['err2']

    A = np.take(A, non_zero_cols, axis=1)
    #print(f'Remaining A : {A}, {A.shape}')
    A = np.take(A, non_zero_rows, axis=0)
    #print(f'Remaining A : {A}, {A.shape}')
    #print('y: ', y)
    #print('Non-zero rows:', non_zero_rows)
    #print('Non-zero rows len:', non_zero_rows.shape)
    #print('Shape of remaining A:', A.shape)
    #print('Remaining A: ', A)

    y = y[non_zero_rows]
    #print('Remaining y:', y)
    
    x = self.conc
    x = x[non_zero_cols]
    #print('Remaining x:', x)
    arr = (x>0).astype(np.int32)
    # Now solve using this new A and y

    # parameters d and s do not matter. They are not used in the algorithm
    n = A.shape[1]
    t = A.shape[0]
    d = self.d
    s = self.s
    l = 0.1
    #print(f'non_zero_cols: {non_zero_cols}')
    #print(f'non_zero_rows: {non_zero_rows}')
    #print(f'y: {y}')
    #print(f' t = {t}, n = {n}')
    #print(f'A = {A}, {A.shape}')

    # In case of invalid input, we may have the case that some rows are positive
    # but all columns are taken out. It should never happen that all rows are
    # negative but some columns remain from output of COMP. Hence second assertion
    # is invalid.
    if t == 0:
      assert n == 0
    #elif n == 0:
    #  assert t == 0
      
    infected = np.zeros(self.n)
    answer = np.zeros(self.n)
    prob1_new = np.zeros(self.n)
    prob0_new = np.ones(self.n)
    determined = 1
    overdetermined = 0
    
    #print("Retest:", len(non_zero_cols))
    
    # Calling internal algo is needed only when there is at least one infection
    #
    # It is also avoided when there are any discrepancies in the test. Only COMP
    # output is returned. Discrepancies usually happen due to spurious testing.
    # This behaviour may change later.
    if A.size != 0 and total_errors == 0:
      # Create another CS class to run the secondary algorithm
      # Better to set mr parameter to None since it depends on number of rows
      # and will change for this internal CS object. frac will be used instead
      # for deciding splits
      _cs = CS(n, t, s, d, l, arr, A, mr=None)
      _cs.conc = x
      
      answer_internal, infected_internal, infected_dd, prob1, prob0, score, tp, fp, fn, _, determined,\
        overdetermined, surep, unsurep, wrongly_undetected, num_infected_in_test, retest, Ambi, lamda =\
        _cs.decode_lasso(y, secondary_algo, compute_stats=compute_stats)
      
      for ans, val, idx in zip(answer_internal, infected_internal, non_zero_cols):
        infected[idx] = val
        answer[idx] = ans

      for p1, p0, idx in zip(prob1, prob0, non_zero_cols):
        prob1_new[idx] = p1
        prob0_new[idx] = p0
      
    if test:
      return infected, prob1_new, prob0_new, score, tp, fp, fn
    else:
      return answer, infected, prob1_new, prob0_new, determined, overdetermined, retest, Ambi, lamda

  def decode_qp(self, results):
    pass

  def print_matrix(self):
    pass

  def pickle_dump(self, filename):
    pass

if __name__ == '__main__':
  raise ValueError('Running experiments has been moved to cs_expts.py. '
      'Either use that or sel_matrix.py')

