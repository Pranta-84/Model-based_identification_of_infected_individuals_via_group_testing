# Additions on tapestry code

Most of our additions are to the function `decode_lasso()` in `cs.py`, where we add non-tapestry algorithms (Dorfman, standard retesting, and P-BEST). We also re-arrange some of the existing lines. None of the additions or re-arrangements alter the flow of operations for the tapestry algorithms; they ensure that the non-tapestry algorithms run correctly under the Ghosh et al. (2021) settings without errors.

A function `Ambicheck()` is added to `cs.py` (lines 1620 - 1668), which checks whether there is ambiguity in the simplified pooling matrix (explained in Section 2.2 of the paper). A few lines (62 - 67, 91 - 95, 101 - 105) are added to the function `get_quantitative_results()` in `cs.py` to generate group viral loads using the Dorfman pooling matrix. A few other lines are added in `cs.py` as well; all added lines are commented.

A few lines (156 - 172) are added to the function `do_single_expt()` in `cs_expts.py` to implement Dorfman testing. We also add the functions `decode_comp_new3()`, `decode_comp_new_Se1()`, `decode_comp_new_Se95by100(), get_infected_dd2(), get_infected_dd_Se1(),` and `get_infected_dd_Se95by100()` to `comp.py` in the `/core/` folder to enable "regular Dorfman" and "optimal Dorfman" testing.

We add `PostProcessingData.py`, which analyze the results of each algorithm presented in `/expt_stats/optimized_M_93_961_kirkman/` folder.
