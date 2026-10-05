# Are we asking too much of the data? Model-based identification of infected individuals via group testing

This repository contains all the software programs used for the paper "Are we asking too much of the data? Model-based identification of infected individuals via group testing" submitted to *Statistics in Medicine*.

The "Ghosh et al. (2021) simulation files" folder contains all the files used to obtain the results presented in Section 3.1 of the paper. Most of files were written by the authors of [Ghosh et al. (2021)](https://ieeexplore.ieee.org/document/9416868). This folder contains their code along with our additions or rearrangements. The original version of their code is available in their GitHub [repository](https://github.com/atoms-to-intelligence/tapestry/tree/master) and we follow the instructions mentioned in that to run their programs. We do not make any changes that alter the implementation of their algorithms. Rather, we add or re-arrange some of their code to implement standard group testing algorithms (explained in section 2.1 of the paper) and P-BEST alongside their algorithms. All of our additions are documented in the file `/Ghosh et al. (2021) simulation files/AdditionsToTapestryCode.md` file. We implement all the algorithms using the program `/Ghosh et al. (2021) simulation files/tools/run_expts.py`. 

The "Abbott SARS-CoV-2 assay R files" folder contains all the R files used to obtain results presented in Section 3.2 of the paper. The program "NeccessaryFunctions.R" contains functions implementing all the algorithms. The program "PaperCodeAbbottAssay.R" implements all the algorithms and can be used to reproduce the results. All the lines are properly commented- please read the comments for better understanding.

The "Clemson assay R files" folder contains all the R files used to obtain results presented in Section 4 of the paper. 

The folder "/PBestMATLABCode/" within each folder contains MATALB programs written by Shental et al. (2021). We do not make any changes to these programs for our investigation. To implement P-BEST 2021, we use all the programs, whereas to implement [P-BEST 2024](https://doi.org/10.1038/s43856-024-00531-w), we only use "applyGPSR.m". These programs can also be accessed from the "/mfiles/" folder of the [P-BEST 2021 repository](https://github.com/NoamShental/PBEST). This is the repository corresponding to [Shental et al. (2021)](https://www.science.org/doi/10.1126/sciadv.abc5961).

The "PoolingMatrices" folder contains all the pooling matrices promised in the paper and its web appendix.
