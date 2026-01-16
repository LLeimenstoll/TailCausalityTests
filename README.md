# TailCausalityTests
Code for the paper "General tests for pairwise causality in extremes" by Lisa Leimenstoll and Melanie Schienle

## External code and license

This project includes a modified version of code from:

Gnecco, Nicola, Nicolai Meinshausen, Jonas Peters, and Sebastian Engelke, 2019.
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

**Data sources**

- **Train operations (expected and actual departure/arrival times)**  
  Open public transport data from:
  - opentransportdata.swiss: https://opentransportdata.swiss  
    (real-time / historical timetable and operational data for Swiss trains), downloaded on 13.08.2024

- **Precipitation (Zurich)**  
  Hourly precipitation data from:
  - MeteoSchweiz (MeteoSwiss): https://www.meteoswiss.admin.ch  
    (Swiss Federal Office of Meteorology and Climatology), downloaded on 14.08.2024


**How to obtain the data**

1. **Train data (opentransportdata.swiss)**  
   - Register at https://opentransportdata.swiss and request an API key if needed.  
   - Download historical operational data for the relevant period (May 2021–July 2024) for Swiss long-distance trains.  
   - In the code, this data is read and filtered in scripts such as:
     - `Application/train_processing/01_download_filter_IC1.R`
     - `Application/train_processing/02_build_train_delays.R`

2. **Weather data (MeteoSchweiz)**  
   - Obtain hourly precipitation measurements for the Zurich station (e.g. station code `REH`) from the MeteoSwiss data portal.  
   - The repository assumes you store these data files under `Application/data/weather/`.  
   - The matching of train delays to hourly precipitation is performed in:
     - `Application/train_processing/03_match_train_weather.R`


The corresponding R code is located in `Application/train_processing/` and `Application/`, and the final merged dataset is written as `Application/data/data_combined_train_weather.csv` (not included in the repository).

---

### 2. Precipitation and river discharges (Bavaria, Germany)

**Data sources**

- **River discharge (daily/instantaneous flows)**  
  - Bayerisches Landesamt für Umwelt (LfU): https://www.lfu.bayern.de  
    (river gauge data, e.g. Danube at Passau, Main at Würzburg/Schweinfurt)
    - Donau/Passau, download on 25.11.2024
    - Main/Würzburg, download on 26.11.2024
    - Main/Schweinfurt, download on 03.09.2025

- **Precipitation (daily totals)**  
  - Deutscher Wetterdienst (DWD): https://www.dwd.de  
    (German National Meteorological Service – Climate Data Center, daily precipitation)

**How to obtain the data**

1. **River discharge data**

   - From the LfU portal, download discharge time series
   - Export them as CSV and place under `Application/data/` (e.g. `Passau.csv`, `Wurzburg.csv`, `Schweinfurth.csv`).

2. **Precipitation data**

   - From the DWD Climate Data Center, download daily precipitation series for stations covering the corresponding catchments (e.g. Passau, Würzburg).
   - Store them under `Application/data/` (e.g. `passau_precipitation.txt`, `wurzburg_precipitation.txt`).

The corresponding R code can be found in:

- `Application/river_processing/riverflows_precipitation_analysis.R`

---

### 3. Financial stock markets and cryptocurrencies

**Data sources**

- **S&P 500 index prices (open, close)**  
  - Yahoo Finance: https://finance.yahoo.com  
    (symbol `^GSPC`, accessed via `quantmod::getSymbols()`)

- **Bitcoin prices**  
  - Coin Metrics: https://coinmetrics.io  
    (e.g. daily close prices from their market data CSV/API)
    donwloaded on 06.12.2024

- **Volatility and equity indices as confounders**

  - CBOE Volatility Index (VIX):
    - Yahoo Finance: symbol `^VIX`, via `quantmod::getSymbols()`.
  - MSCI Europe Index:
    - MSCI (2023) – daily index levels (MSCI Europe).
      - Must be obtained from MSCI’s data services, stored locally (e.g. `Application/data/MSCI_eu.csv`).
      - Downloaded on 08.05.2025

**How to obtain the data**

1. **S&P 500 and VIX via Yahoo Finance**

   - The scripts use the R package **quantmod**:

     ```r
     getSymbols("^GSPC", from = startdate, to = enddate, auto.assign = TRUE)
     getSymbols("^VIX",  from = startdate, to = enddate, auto.assign = TRUE)
     ```

   - You can adjust `startdate` and `enddate` in the scripts under `Application/finance_processing/`.

2. **Bitcoin prices (Coin Metrics)**

   - Download historical BTC price data (e.g. daily close) from Coin Metrics:
     - https://coinmetrics.io
   - Save as CSV (e.g. `coin-metrics.csv`) under `Application/data/`.

3. **MSCI Europe index**

   - Obtain daily MSCI Europe index levels for the corresponding period from MSCI’s data services:
     - https://www.msci.com  
   - Store as `Application/data/MSCI_eu.csv` with appropriate date and index columns.

The main R code for this application is located in:

- `Application/finance_processing/finance_causal_tail_analysis.R`

This script shows the complete pipeline from data retrieval via `quantmod` and CSV input to the tail-causal inference results.

---

For all three applications, the raw data are **not** redistributed in this repository due to licensing and size constraints. Instead, this README and the corresponding scripts indicate exactly which external sources to use and how to prepare the data in order to reproduce the analyses.
