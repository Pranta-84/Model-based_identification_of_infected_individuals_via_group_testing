# Are we asking too much of the data? Model-based identification of infected individuals via group testing

This repository contains all the software programs used for the paper "Are we asking too much of the data? Model-based identification of infected individuals via group testing" submitted to *Statistics in Medicine*.

The "Abbott SARS-CoV-2 assay R files" folder contains all the R files used to obtain results presented in Section 3.2 of the paper. The program "NeccessaryFunctions.R" contains functions implementing all the algorithms. The program "PaperCodeAbbottAssay.R" implements all the algorithms and can be used to reproduce the results. All the lines are properly commented- please read the comments for better understanding.

The "Clemson assay R files" folder contains all the R files used to obtain results presented in Section 4 of the paper. 

The folder "/PBestMATLABCode/" within each folder contains MATALB programs written by Shental et al. (2021). We do not make any changes to these programs for our investigation. To implement P-BEST 2021, we use all the programs, whereas to implement [P-BEST 2024](https://doi.org/10.1038/s43856-024-00531-w), we only use "applyGPSR.m". These programs can also be accessed from the "/mfiles/" folder of the [P-BEST 2021 repo](https://github.com/NoamShental/PBEST). This is the repository corresponding to [Shental et al. (2021)](https://www.science.org/doi/10.1126/sciadv.abc5961).

The "PoolingMatrices" folder contains all the pooling matrices promised in the paper and its web appendix.
