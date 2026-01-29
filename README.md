# TailCausalityTests
Code for the paper "General tests for pairwise causality in extremes" by Lisa Leimenstoll and Melanie Schienle

## External code and license

This project includes a modified version of code from:

Nicola Gnecco, Nicolai Meinshausen, Jonas Peters, and Sebastian Engelke, 2019.
Original repository: https://github.com/nicolagnecco/causalXtreme
License: GPL-3.0

The adapted code is located in `Simulation_Study/k_simulation.R` and is distributed under
the same GPL-3.0 license, with minor modifications:
  - add second tail index and pareto distribution in simulate_data function
  - add pareto distribution in simulate_noise
  - in simulation_0 change arguments, add second tail index
  - add second tail index and pareto distribution in my_args
  - change of plot labels and adjust for two tail indices
  - add additional simualtion study for percentage of wrong causal inference 
    between two variables

An adapeted version of the causal_tail_coefficient function is located in functions.R`, with minor modifications:

- add of min argument to analyse the lower tail

And a modified version from:
O. C. Pasche, V. Chavez-Demoulin and A. C. Davison. 2020.
https://github.com/opasche/ExtremalCausalModelling 
License: GPL-3.0

The adapeted functions are located in functions.R`, with minor modifications:
- Change of test statistic in caustal_tail_perumation_test to test two sided

## Data sources and preprocessing

This repository contains three main applications of the causal tail tests, each based on external data sources. The code in this repository does **not** redistribute the raw data. Instead, we provide scripts to download (where possible) or describe how to obtain the data, and to reproduce the preprocessing steps.

### 1. Precipitation and train delays in Switzerland

- **Train departure times**  
  Open public transport data from:
  - opentransportdata.swiss: https://opentransportdata.swiss  
    (real-time / historical timetable and operational data for Swiss trains), downloaded on 13.08.2024
  - Download historical operational data for the relevant period (May 2021–July 2024)
  - In the code, this data is read and filtered in scripts such as:
    - `Application/train_processing/01_filter_train_data_switzerland.R`
    - `Application/train_processing/02_prep_train_data_switzerland.R`


- **Precipitation**  
  Hourly precipitation data from:
  - MeteoSchweiz (MeteoSwiss): https://www.meteoswiss.admin.ch  
    (Swiss Federal Office of Meteorology and Climatology), downloaded on 14.08.2024
  - Obtain hourly precipitation measurements for the Zurich station (e.g. station code `REH`) from the MeteoSwiss data portal.  
  - The repository assumes you store these data files under `Application/data/weather/`.  
  - The matching of train delays to hourly precipitation is performed in:
    - `Application/train_processing/03_match_train_weather.R`


The corresponding R code is located in `Application/train_processing/` and `Application/`, and the final merged dataset is written as `Application/data/data_combined_train_weather.csv` (not included in the repository).

---

### 2. Precipitation and river discharges (Bavaria, Germany)


- **River discharge**  
  - Bayerisches Landesamt für Umwelt (LfU): https://www.lfu.bayern.de  
    (river gauge data, e.g. Danube at Passau, Main at Würzburg/Schweinfurt)
    - Donau/Passau, download on 25.11.2024
    - Main/Würzburg, download on 26.11.2024
    - Main/Schweinfurt, download on 03.09.2025
  - Export as CSV and place under `Application/data/` (e.g. `passau.csv`, `wurzburg.csv`, `schweinfurth.csv`)

- **Precipitation**  
  - Deutscher Wetterdienst (DWD): https://www.dwd.de  
    (German National Meteorological Service – Climate Data Center, daily precipitation)
  - Download daily precipitation series for station IDs `03878` and `05705`,  and store them under `Application/data/` (e.g. `passau_precipitation.txt`, `wurzburg_precipitation.txt`).


The corresponding R code can be found in:

- `Application/riverflow_application.R`

---

### 3. Financial stock markets and cryptocurrencies


- **S&P 500 index prices (open, close)**  
  - Yahoo Finance: https://finance.yahoo.com  
    (symbol `^GSPC`, accessed in R via `quantmod::getSymbols()`)

- **Bitcoin prices**  
  - Coin Metrics: https://coinmetrics.io  
    (e.g. daily close prices from their market data CSV/API)
  -  stored locally ('Application/data/coin-metrics.csv`)
  -  downloaded on 06.12.2024

- **Volatility and equity indices as confounders**

  - CBOE Volatility Index (VIX):
    - Yahoo Finance: symbol `^VIX`, via R `quantmod::getSymbols()`.
  - MSCI Europe Index:
    - MSCI (2023) – daily index levels (MSCI Europe).
    - Must be obtained from MSCI’s data services, stored locally (e.g. `Application/data/MSCI_eu.csv`).
    - Downloaded on 08.05.2025
      

The main R code for this application is located in:

- `Application/finance_processing/finance_causal_tail_analysis.R`


---

For all three applications, the raw data are **not** redistributed in this repository due to licensing and size constraints. Instead, this README and the corresponding scripts indicate exactly which external sources to use and how to prepare the data in order to reproduce the analyses.

## R files

