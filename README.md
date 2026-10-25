# Are we asking too much of the data? Model-based identification of infected individuals via group testing

This repository contains all the software programs used for the paper "Are we asking too much of the data? Model-based identification of infected individuals via group testing" submitted to *Statistics in Medicine*.

### Setup
To run these programs, you need R, python, and MATLAB (with a valid license) installed on your device. Because the programs use all three simultaneously, the three must be able to communicate with one another. The steps below are one way to start the setup.

First create a conda environment named "r-matlab311" with Python 3.11 by running the following in the terminal:
```{}
conda create -n r-matlab311 python=3.11 numpy
```
Next, activate the environment:
```{}
conda activate r-matlab311
```
Finally, install the matlab engine for Python within the environment:
```{}
pip install matlabengine==24.1.*
```

Note that these steps are a starting point for setting up a R, python, and MATLAB connection on your device. Some commands may need to be adjusted for your own environment. Once the setup is complete, all the programs are expected to run without any errors. 

### Ghosh et al. (2021) simulation

The "Ghosh et al. (2021) simulation files" folder contains all the files used to obtain the results presented in Section 3.1 of the paper. Most of files were written by the authors of [Ghosh et al. (2021)](https://ieeexplore.ieee.org/document/9416868). This folder contains their code along with our additions or rearrangements. The original version of their code is available in their GitHub [repository](https://github.com/atoms-to-intelligence/tapestry/tree/master), and we follow the instructions given there to run their programs. We do not make any changes that alter the implementation of their algorithms. Rather, we add or re-arrange some of their code to implement standard group testing algorithms (explained in section 2.1 of the paper) and P-BEST alongside their algorithms. All of our additions are documented in the file `/Ghosh et al. (2021) simulation files/AdditionsToTapestryCode.md` file. We implement all the algorithms using the program `/Ghosh et al. (2021) simulation files/tools/run_expts.py`. 

### Abbott SARS-CoV-2 assay simulation
The "Abbott SARS-CoV-2 assay R files" folder contains all the R files used to obtain results presented in Section 3.2 of the paper. The program "NecessaryFunctions.R" contains functions implementing all the algorithms. The program "PaperCodeAbbottAssay.R" implements all the algorithms and can be used to reproduce the results. All the lines are commented. Please read these comments for better understanding.

### Clemson assay
The "Clemson assay R files" folder contains all the R files used to obtain results presented in Section 4 of the paper. The file "ReddiLabLognormalData.csv" contains the approximated REDDI data. The program "FunctionsForApproximatedREDDILabDataAnalysis.R" contains functions implementing each algorithm. The program "ApproximatedREDDILabDataAnalysis.R" implements all the algorithms and can be used to reproduce the results. All the lines are commented. Please read these comments for better understanding.

### General information
The folder "/PBestMATLABCode/" within each folder contains MATLAB programs written by Shental et al. (2021). We do not make any changes to these programs for our investigation. To implement P-BEST 2021, we use all the programs, whereas to implement [P-BEST 2024](https://doi.org/10.1038/s43856-024-00531-w), we only use "applyGPSR.m". These programs can also be accessed from the "/mfiles/" folder of the [P-BEST 2021 repository](https://github.com/NoamShental/PBEST). This is the repository corresponding to [Shental et al. (2021)](https://www.science.org/doi/10.1126/sciadv.abc5961). The "PoolingMatrices" folder contains all the pooling matrices promised in the paper and its web appendix.

### Contact
Please contact pranta.du.stat@gmail.com if you have any questions.
