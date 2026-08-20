# TailCausalityTests

Code for the paper **"Identification and Inference for Causal Effects in Extremes under General Conditions"** by Lisa Leimenstoll and Melanie Schienle.

It includes:
- Simulation studies for causal tail tests in R in `Siumlation_Study/`,
- Applications to train delays, river flows, and financial data in `Application/`

## Abstract 

Understanding the propagation of extreme events is important in many economic and environmental applications, yet most econometric methods for causal inference focus on average effects rather than tail behavior. This paper studies the identification of causal relations in extremes within a linear structural causal model associated with a directed acyclic graph. We analyze the asymptotic behavior of the Causal Tail Coefficient (CTC), a measure of causal dependence between extreme realizations of variables, when the innovations of the structural model follow regularly varying distributions. In contrast to the existing literature, we allow the variables in the system to exhibit heterogeneous tail indices and consider the presence of potentially heavy-tailed confounders. We derive theoretical results characterizing the limiting behavior of the CTC under these conditions and show how differences in tail behavior can provide identifying information about the causal structure. Based on these results, we propose a testing framework that distinguishes direct causal effects from confounding in heavy-tailed environments. The finite-sample performance of the proposed methodology is investigated through simulation studies and illustrated using applications to climate and financial data. The results demonstrate that causal relations among extreme events can be identified under the proposed assumptions, even when variables exhibit substantially different tail behavior.

---

## R code structure (scripts and simulations)

The R code is organized as follows:

### Core methods

- `functions.R`  
  Helper functions for tail indices, causal tail coefficients, and tests.

- `testing_strategy.R`  
  Functions to test for causal tail direction and confounding.

### Simulation studies

The following scripts in `Application/Simulation_Study/` generate the results reported in Section 3 of the paper, ordered according to the section numbering:

- `convergence_simulation.R`  
  Produces the results in **Section 5.2 (CTC Estimates: Finite-Sample Rate of Convergence)**.

- `simulation_study_configurations.R`  
  Generates the simulation data for **Section 5.3 (Finite-Sample Performance of the Proposed Tests)**.

- `simulation_study_configurations_evaluation.R`  
  Produces the results and plots for **Section 5.3 (Finite-Sample Performance of the Proposed Tests)**.

- `confounder_simulation.R`  
  Produces the results in **5.3.3 (Confounder-Test)**.

- `k_simulation.R`  
  Produces the results in **5.1 (Choice of the Tuning Parameter k: Independent Pre-Study)**.

- `lingam_pretest_comp.R`  
  Produces the results in **Section 5.2 (CTC Estimates: Finite-Sample Rate of Convergence)**.

---

## Data sources and preprocessing

This repository contains three main applications of the causal tail tests, each based on external data sources. The code in this repository does **not** redistribute the raw data. Instead, we provide scripts to download (where possible) or describe how to obtain the data, and to reproduce the preprocessing steps.

### 1. Precipitation and Train Delays in Switzerland

- **Train departure times**  
  Open public transport data from:
  - opentransportdata.swiss: https://opentransportdata.swiss  
    (real-time / historical timetable and operational data for Swiss trains), downloaded on 13.08.2024.
  - Download historical operational data for the relevant period (May 2021–July 2024).
  - In the code, this data is read and filtered in:
    - `Application/train_processing/01_filter_train_data_switzerland.R`
    

- **Precipitation**  
  Hourly precipitation data from:
  - MeteoSchweiz (MeteoSwiss): https://www.meteoswiss.admin.ch  
    (Swiss Federal Office of Meteorology and Climatology), downloaded on 14.08.2024.
  - Obtain hourly precipitation measurements for the Zurich station (e.g. station code `REH`) from the MeteoSwiss data portal.  
  - Store these data files under `Application/data/weather/`.  
  - Matching of train delays to hourly precipitation is performed in:
   - `Application/train_processing/02_prep_train_data_switzerland.R`

The corresponding R code is located in `Application/train_processing/` and `Application/`.  
The final merged dataset is written as `Application/data/data_combined_train_weather.csv` (not included in the repository).

---

### 2. Precipitation and River Flows (Bavaria, Germany)

- **River discharge**  
  - Bayerisches Landesamt für Umwelt (LfU): https://www.lfu.bayern.de  
    (river gauge data, e.g. Danube at Passau, Main at Würzburg/Schweinfurt)
    - Donau/Passau, downloaded on 25.11.2024  
    - Main/Würzburg, downloaded on 26.11.2024  
    - Main/Schweinfurt, downloaded on 03.09.2025  
  - Export as CSV and place under `Application/data/`  
    (e.g. `passau.csv`, `wurzburg.csv`, `schweinfurt.csv`).

- **Precipitation**  
  - Deutscher Wetterdienst (DWD): https://www.dwd.de  
    (German National Meteorological Service – Climate Data Center, daily precipitation).
  - Download daily precipitation series for station IDs `03878` and `05705`, and store them under `Application/data/`  
    (e.g. `passau_precipitation.txt`, `wurzburg_precipitation.txt`).

The corresponding R code can be found in:
- `Application/riverflow_application.R`

---

### 3. Financial Stock Markets and Cryptocurrencies

- **S&P 500 index prices (open, close)**  
  - Yahoo Finance: https://finance.yahoo.com  
    (symbol `^GSPC`, accessed in R via `quantmod::getSymbols()`).

- **Bitcoin prices**  
  - Coin Metrics: https://coinmetrics.io  
    (e.g. daily close prices from their market data CSV/API).  
  - Stored locally (`Application/data/coin-metrics.csv`), downloaded on 06.12.2024.

- **Volatility and equity indices as confounders**
  - **CBOE Volatility Index (VIX)**  
    - Yahoo Finance: symbol `^VIX`, via R `quantmod::getSymbols()`.
  - **MSCI Europe Index**  
    - MSCI (2023) – daily index levels (MSCI Europe).  
    - Must be obtained from MSCI’s data services and stored locally  
      (e.g. `Application/data/MSCI_eu.csv`), downloaded on 08.05.2025.

The main R code for this application is located in:
- `Application/finance_processing/finance_causal_tail_analysis.R`

---


## External code and license

This project includes adapted versions of code from the following sources:

1. **Gnecco et al. (2019) – causalXtreme**  
   - Nicola Gnecco, Nicolai Meinshausen, Jonas Peters, and Sebastian Engelke, 2019.  
   - Original repository: https://github.com/nicolagnecco/causalXtreme  
   - License: GPL-3.0  

   **Adapted code in this repository:**
   - `Application/Simulation_Study/k_simulation.R`  
     The code is distributed under the same GPL-3.0 license, with minor modifications:
       - Addition of a second tail index and Pareto distribution in `simulate_data`.
       - Addition of Pareto distribution in `simulate_noise`.
       - In `simulation_0`: change of arguments, addition of a second tail index.
       - Addition of a second tail index and Pareto distribution in `my_args`.
       - Changes of plot labels and adjustment for two tail indices.
       - Additional simulation study for the percentage of wrong causal inference between two variables.
   - An adapted version of the `causal_tail_coefficient` function is located in `functions.R`, with the following modification:
       - Addition of a `min` argument to analyse the lower tail.

2. **Pasche et al. (2020) – ExtremalCausalModelling**  
   - O. C. Pasche, V. Chavez-Demoulin, and A. C. Davison, 2020.  
   - Original repository: https://github.com/opasche/ExtremalCausalModelling  

   **Adapted code in this repository:**
   - Adapted functions are located in `functions.R`, with minor modifications:
       - Change of the test statistic in `causal_tail_permutation_test` to test two-sided.
