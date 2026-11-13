# dinkelbachbounds
This repository contains empirical application code and estimator construction for the paper "Partial Identification and Estimation in a Class of Potential Outcomes Models: A Convex-Fractional Approach". Applications include JobCorps (Lee, 2009, Section 8), and Thinking Healthy Program (THP, Baranov et al., 2020) and the Primary School Deworming Project (Miguel and Kremer, 2004, Section 9). Some implementation choices are recorded in the file `ASSUMPTIONS.md`. More detailed instructions are available in the file `dinkelbach_codeguide.pdf`.

## Overview
We separate the scripts into three categories:

1. Data preparation (`data_prep.py`, `data_prep\_baranov.py`, `data_prep_deworming.py`): each corresponding file loads the data set, and maps it to treatment $D$, selection $S$, outcome $Y$, covariates $X$, weights and clusters.
2. Estimation (`estimator.py`, `lp_bounds.py`): the first file includes the first stage estimation, the cross-fitted localized one-step estimator, and inference; the second file computes the classical Lee bounds, exact marginal (no-covariate) bounds, and the discrete LP heuristic of Section 2.3.
3. Runner scripts denoted as `run_*.py` create and store the results from the paper. The following table is a short summary.

| Script | Output | Output files | Paper references
|---|---|---|---|
| `run_section8.py` | JobCorps share bounds with bands; bounds by identification regime; complier rescaling; Masten-Pourier frontier | `table1_validation.csv`, `table5_bounds.csv`, `table7_compliers.csv`, `mp_frontier.csv` | Tables 1, 5, 7 (Sections 3.1, 8) |
|`run_table2.py` | JobCorps discrete-heuristic LP bounds by number of wage bins | `table2_lp_granularity.csv` | Table 3 (Sections 2.3, 3.2) |
|`run_baranov.py` | THP: ITT, Lee, exact marginal and covariate bounds for every outcome | `baranov_bounds.csv` | Tables 8A, 9 (Section 9) |
|`run_deworming.py` | Deworming: ITT, Lee, exact marginal and covariate bounds; share bounds by specification |`deworming_bounds.csv`, `deworming_table1.csv` | Tables 8C, 11 (Section 9) |
|`run\_app\_tables.py` | THP analogues of Tables 1 and 3 | `app\_table1\_analogues.csv`, `app_table3_analogues.csv` | not in the paper 

## How to run the code
Each runner script starts from the folder that contains the scripts, e.g. `python run_section8.py`. All runner scripts create and write to `../results/`, a path relative to the *working* directory. Note that the data paths are computer-dependent, stored in the data-preparation modules (`DATA_PATH`, `CSV`), and therefore must be edited once the files are run on another computer. Approximate running times on a M1 MacBook Pro: `run_section8.py` about 30 minutes (of which Table 1 bands with Lee's 28 covariates take about 23), `run_table2.py` under one minute, `run_deworming.py` about 40 seconds.

## References
Victoria Baranov, Sonia Bhalotra, Pietro Biroli, and Joanna Maselko. *Maternal depression, women’s empowerment, and parental investment: Evidence from a randomized controlled trial.* American Economic Review, 110 (3):824–59, March 2020.

David S. Lee. *Training, wages, and sample selection: Estimating sharp bounds on treatment effects.* Review of Economic Studies, 76(3):1071–1102, 2009.

Edward Miguel and Michael Kremer. *Worms: Identifying impacts on education and health in the presence of treatment externalities.* Econometrica, 72(1):159–217, 2004.