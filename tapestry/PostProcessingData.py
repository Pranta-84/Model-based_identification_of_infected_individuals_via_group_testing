# This file is for calculating E(T)/I, PPA, NPA, PPV, and NPV using the results
## obtained with run_expts.py

# Importing required packages
import gzip
import pickle
import pandas as pd
import numpy as np
import statistics

# load the file
with gzip.open("/home/bilder/pranta/tapestry/expt_stats/optimized_M_93_961_kirkman/PBest-2024/15/expt_stats.p.gz", "rb") as f:
    data = pickle.load(f)

# Creating a place to store the quantities needed
rows = []
# Extracting required quantities and storing it to the object rows
for res in data:
    rows.append({
        "tp": res["tp"],
        "fp": res["fp"],
        "fn": res["fn"],
        "score": res["score"],
        "rmse": res["rmse"],
        "retest": res["retest"],
        "Ambi": res["Ambi"],
        "lambda": res["lamda"],
        "n": res["n"],
        "d": res["d"],
        "t": res["t"],
        "algo": res["algo"],
        "mlabel": res["mlabel"]
    })

# Converting the object rows to a a DataFrame
df = pd.DataFrame(rows)
# Number of rows or observations in the dataframe "df"
df.shape[0]

# Calculating Ambiguity proportion
np.sum(df["Ambi"])/10000

# Calculating the average value of the regularization
## Only applicable for algorithms "combined_COMP_l1ls" and "combined_COMP_lasso" 
statistics.mean(df["lambda"])

###### 5 positives/961 ####
# Mean of lambda with Least squares is 0.1
# Mean of lambda with lasso is 0.00126, maximum is 0.1547

###### 10 positives/961 ####
# Mean of lambda with Least squares is 0.1
# Mean of lambda with lasso is 0.00126, maximum is 0.0808

###### 15 positives/961 ####
# Mean of lambda with Least squares is 0.1
# Mean of lambda with lasso is 0.0069, maximum is 0.1068

###### 20 positives/961 ####
# Mean of lambda with Least squares is 0.1
# Mean of lambda with lasso is 0.0189, maximum is 0.1165


## Only for optimal Dorfman at Se = 1
### Need to update the pooling matrix based on number of positives per 961
DT_opt_Se1 = pd.read_csv("DorfMatrix_p15by961_SeSp1.csv").to_numpy()

### E(T)
ET = DT_opt_Se1.shape[0] + statistics.mean(df['retest'])
ET
### E(T)/I
ET/961

## Only for optimal Dorfman at Se = 0.95
### Need to update the pooling matrix based on number of positives per 961
DT_opt_Se95by100 = pd.read_csv("DorfMatrix_p10by961_Se0.95Sp1.csv").to_numpy()
### E(T)
ET = DT_opt_Se95by100.shape[0] + statistics.mean(df['retest'])
ET
### E(T)/I
ET/961



## Total number of retest
np.sum(df['retest'])

## Calculating avgerage number of retest 
avg_retest = statistics.mean(df['retest'])
avg_retest

## Only for Dorfman testing
### E(T)
ET_Dorf = 31 + avg_retest
ET_Dorf
### E(T)/I
ET_Dorf/961

# For all other algorithms
## E(T)
ET = 93 + statistics.mean(df['retest'])
ET
## E(T)/I
ET/961

# Calculating PPA
correct_pos = sum(df["tp"])
Total_true_pos = sum(df['d'])
PPA = correct_pos/Total_true_pos
PPA

# Calculating NPA
(961*10000 - sum(df['d']) - sum(df['fp']))/sum(961 - df['d'])

# Calculating PPV
Total_pred_pos = correct_pos + sum(df['fp'])
correct_pos/Total_pred_pos

# Calculating NPV
Total_pred_neg = (961*10000) - Total_pred_pos
(Total_pred_neg - sum(df['fn']))/Total_pred_neg


# Storing the results of all the algorithms for the 93 by 961 pooling matrix when
## there are 5 positive specimens per 961 specimens.
data_with5 = {
    "Algorithms": ["SBL", "NN-Least Square", "NNOMP", "Regular testing (retest upon ambiguity)",
    "Regular testing (always retest)", "Dorfman testing", "Opt_Dorfman_at_Se1", "Opt_Dorfman_at_Se0.95",
    "NN-LASSO", "PBest_reduced", "PBest_full", "P-BEST 2024"],
    "Cutoff": [0, 0, 0, 0, 0, 0, 0, 0, 0, pd.NA, pd.NA, pd.NA],
    "Number_of_positives_per_data":[5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    "Ambiguity_prop":[0.2143, 0.2143, 0.2143, 0.2143, pd.NA, pd.NA, pd.NA, pd.NA, 0.2143, 0.2143, pd.NA, pd.NA],
    "ET": [93, 93, 93, 94.77, 99.58, 176.82, 136.99, 137.84, 93, 93, 93, 121.51],
    "ET.per": [0.0967, 0.0967, 0.0967, 0.0986, 0.1036, 0.1840, 0.1426, 0.1434, 0.0967, 0.0967, 0.0967, 0.1264],
    "PPA": [0.9997, 1, 0.9993, 1, 1, 1, 1, 1, 0.9939, 1, 0.9990, 0.9838],
    "NPA": [0.9991, 0.9983, 0.9991, 0.9991, 1, 1, 1, 1, 0.9993, 0.9984, 0.9986, 1],
    "PPV": [0.8544, 0.7594, 0.8600, 0.8497, 1, 1, 1, 1, 0.8804, 0.7643, 0.7871, 1],
    "NPV": [0.9999, 1, 0.9999, 1, 1, 1, 1, 1, 0.9999, 1, 0.9999, 0.9999]
}

# Saving as dataframe 
res_data5 = pd.DataFrame(data_with5)
# Writing the file as .csv 
res_data5.to_csv("resultsTapestryandRegular_EachWith5.csv", index = False)

# Storing the false positives of all the algorithms for the 93 by 961 pooling matrix when
## there are 5 positive specimens per 961 specimens.
data_fp_5 = {
    "Algorithms": ["SBL", "NN-Least Square", "NNOMP", "Regular testing (retest upon ambiguity)", "Regular testing (always retest)", "Dorfman testing", "NN-LASSO"],
    "Number_of_positives_per_data":[5, 5, 5, 5, 5, 5, 5],
    "Min":[0, 0, 0, 0, 0, 0, 0],
    "Q1": [0, 1, 0, 0, 0, 0, 0],
    "Q2": [1, 1, 1, 1, 0, 0, 0],
    "Mean": [0.8607, 1.6010, 0.8340, 0.8720, 0, 0, 0.6783],
    "Q3": [1, 2, 1, 2, 0, 0, 1],
    "Max": [6, 8, 6, 3, 0, 0, 6],
    "Total": [8607, 16010, 8340, 8720, 0, 0, 6783]
}
# Saving as dataframe
res_data5_fp = pd.DataFrame(data_fp_5)
# Writing the file as .csv
res_data5_fp.to_csv("resultsTapestryandRegular_EachWith5_fp.csv", index=False)


# Storing the results of all the algorithms for the 93 by 961 pooling matrix when
## there are 10 positive specimens per 961 specimens.
data_with10 = {
    "Algorithms": ["SBL", "NN-Least Square", "NNOMP", "Regular testing (retest upon ambiguity)",
    "Regular testing (always retest)", "Dorfman testing", "Opt_Dorfman_at_Se1", "Opt_Dorfman_at_Se0.95",
    "NN-LASSO", "PBest_reduced", "PBest_full", "P-BEST 2024"],
    "Cutoff": [0, 0, 0, 0, 0, 0, 0, 0, 0, pd.NA, pd.NA, pd.NA],
    "Number_of_positives_per_data":[10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10],
    "Ambiguity_prop":[0.9998, 0.9998, 0.9998, 0.9998, pd.NA, pd.NA, pd.NA, pd.NA, 0.9998, 0.9998, pd.NA, pd.NA],
    "ET": [93, 93, 93, 118.24, 118.25, 300.51, 192.72, 192.79, 93, 93, 93, 137.95],
    "ET.per": [0.0967, 0.0967, 0.0967, 0.1230, 0.1230, 0.3127, 0.2005, 0.2006, 0.0967, 0.0967, 0.0967, 0.1435],
    "PPA": [0.9928, 1, 0.9772, 1, 1, 1, 1, 1, 0.9113, 0.9551, 0.9239, 0.9629],
    "NPA": [0.9914, 0.9840, 0.9918, 0.9999, 1, 1, 1, 1, 0.9935, 0.9899, 0.9920, 0.9999],
    "PPV": [0.5476, 0.3961, 0.5555, 0.9999, 1, 1, 1, 1, 0.5947, 0.4997, 0.5481, 0.9999],
    "NPV": [0.9999, 1, 0.9998, 1, 1, 1, 1, 1, 0.9991, 0.9995, 0.9992, 0.9996]
}
# Saving as dataframe
res_data10 = pd.DataFrame(data_with10)
# Writing the file as .csv
res_data10.to_csv("resultsTapestryandRegular_EachWith10.csv", index = False)



# Storing the results of all the algorithms for the 93 by 961 pooling matrix when
## there are 15 positive specimens per 961 specimens.
data_with15 = {
    "Algorithms": ["SBL", "NN-Least Square", "NNOMP", "Regular testing (retest upon ambiguity)",
    "Regular testing (always retest)", "Dorfman testing", "Opt_Dorfman_at_Se1", "Opt_Dorfman_at_Se0.95", 
    "NN-LASSO", "PBest_reduced", "PBest_full", "P-BEST 2024"],
    "Cutoff": [0, 0, 0, 0, 0, 0, 0, 0, 0, pd.NA, pd.NA, pd.NA],
    "Number_of_positives_per_data":[15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15],
    "Ambiguity_prop":[1, 1, 1, 1, pd.NA, pd.NA, pd.NA, pd.NA, 1, 1, pd.NA, pd.NA],
    "ET": [93, 93, 93, 153.75, 153.75, 406.21, 234.13, 234.13, 93, 93, 93, 146.58],
    "ET.per": [0.0967, 0.0967, 0.0967, 0.1599, 0.1599, 0.4227, 0.2436, 0.2436, 0.0967, 0.0967, 0.0967, 0.1525],
    "PPA": [0.9803, 1, 0.8968, 1, 1, 1, 1, 1, 0.6605, 0.5168, 0.5172, 0.9094],
    "NPA": [0.9735, 0.9516, 0.9799, 1, 1, 1, 1, 1, 0.9890, 0.9870, 0.9871, 0.9999],
    "PPV": [0.3701, 0.2469, 0.4139, 1, 1, 1, 1, 1, 0.4879, 0.3876, 0.3883, 0.9978],
    "NPV": [0.9997, 1, 0.9983, 1, 1, 1, 1, 1, 0.9946, 0.9923, 0.9923, 0.9986]
}
# Saving as dataframe
res_data15 = pd.DataFrame(data_with15)
# Writing the file as .csv
res_data15.to_csv("resultsTapestryandRegular_EachWith15.csv", index=False)


# Storing the false positives of all the algorithms for the 93 by 961 pooling matrix when
## there are 15 positive specimens per 961 specimens.
data_fp_15 = {
    "Algorithms": ["SBL", "NN-Least Square", "NNOMP", "Regular testing (retest upon ambiguity)", "Regular testing (always retest)", "Dorfman testing", "NN-LASSO"],
    "Number_of_positives_per_data":[15, 15, 15, 15, 15, 15, 15],
    "Min":[5, 12, 4, 0, 0, 0, 0],
    "Q1": [20, 39, 16, 0, 0, 0, 6],
    "Q2": [25, 46, 19, 0, 0, 0, 10],
    "Mean": [25.008, 45.93, 19.05, 0, 0, 0, 10.37],
    "Q3": [29, 53, 22, 0, 0, 0, 15],
    "Max": [55, 91, 33, 0, 0, 0, 31],
    "Total": [250080, 459255, 190498, 0, 0, 0, 103765]
}
# Saving as dataframe
res_data15_fp = pd.DataFrame(data_fp_15)
# Writing the file as .csv
res_data15_fp.to_csv("resultsTapestryandRegular_EachWith15_fp.csv", index=False)


# Storing the results of all the algorithms for the 93 by 961 pooling matrix when
## there are 20 positive specimens per 961 specimens.
data_with20 = {
    "Algorithms": ["SBL", "NN-Least Square", "NNOMP", "Regular testing (retest upon ambiguity)",
    "Regular testing (always retest)", "Dorfman testing", "Opt_Dorfman_at_Se1", "Opt_Dorfman_at_Se0.95",
    "NN-LASSO", "PBest_reduced", "PBest_full", "PBest-2024"],
    "Cutoff": [0, 0, 0, 0, 0, 0, 0, 0, 0, pd.NA, pd.NA, pd.NA],
    "Number_of_positives_per_data":[20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20],
    "Ambiguity_prop":[1, 1, 1, 1, pd.NA, pd.NA, pd.NA, pd.NA, 1, 1, pd.NA, pd.NA],
    "ET": [93, 93, 93, 205.38, 205.38, 496.67, 269.77, 270.19, 93, 93, 93, 190.73],
    "ET.per": [0.0967, 0.0967, 0.0967, 0.2137, 0.2137, 0.5168, 0.2807, 0.2811, 0.0967, 0.0967, 0.0967, 0.1985],
    "PPA": [0.9378, 1, 0.7112, 1, 1, 1, 1, 1, 0.3609, 0.2784, 0.2784, 0.9482],
    "NPA": [0.9414, 0.9018, 0.9696, 1, 1, 1, 1, 1, 0.9921, 0.9846, 0.9846, 0.9999],
    "PPV": [0.2538, 0.1779, 0.3319, 1, 1, 1, 1, 1, 0.4932, 0.2784, 0.2784, 0.9988],
    "NPV": [0.9986, 1, 0.9937, 1, 1, 1, 1, 1, 0.9865, 0.9846, 0.9846, 0.9989]
}
# Saving as dataframe
res_data_20 = pd.DataFrame(data_with20)
# Writing as .csv
res_data_20.to_csv("resultsTapestryandRegular_EachWith20.csv", index=False)

# Storing the false positives of all the algorithms for the 93 by 961 pooling matrix when
## there are 20 positive specimens per 961 specimens.
data_fp_20 = {
    "Algorithms": ["SBL", "NN-Least Square", "NNOMP", "Regular testing (retest upon ambiguity)", "Regular testing (always retest)", "Dorfman testing", "NN-LASSO"],
    "Number_of_positives_per_data":[20, 20, 20, 20, 20, 20, 20],
    "Min":[18, 37, 9, 0, 0, 0, 0],
    "Q1": [45, 79, 25, 0, 0, 0, 1],
    "Q2": [54, 91, 29, 0, 0, 0, 6],
    "Mean": [54.844, 91.88, 28.54, 0, 0, 0, 7.24],
    "Q3": [64, 104, 32, 0, 0, 0, 12],
    "Max": [111, 173, 47, 0, 0, 0, 34],
    "Total": [548440, 918801, 285406, 0, 0, 0, 72363]
}
# Saving as dataframe
res_fp_20 = pd.DataFrame(data_fp_20)
# Saving as .csv
res_fp_20.to_csv("resultsTapestryandRegular_EachWith20_fp.csv", index=False)





