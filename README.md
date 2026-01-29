# TailCausalityTests

Code for the paper **"General tests for pairwise causality in extremes"** by Lisa Leimenstoll and Melanie Schienle.

It includes:
- Simulation studies for causal tail tests in R in `Siumlation_Study/`,
- Applications to train delays, river flows, and financial data in `Application/`,

---

## Data sources and preprocessing

This repository contains three main applications of the causal tail tests, each based on external data sources. The code in this repository does **not** redistribute the raw data. Instead, we provide scripts to download (where possible) or describe how to obtain the data, and to reproduce the preprocessing steps.

### 1. Precipitation and train delays in Switzerland

- **Train departure times**  
  Open public transport data from:
  - opentransportdata.swiss: https://opentransportdata.swiss  
    (real-time / historical timetable and operational data for Swiss trains), downloaded on 13.08.2024.
  - Download historical operational data for the relevant period (May 2021–July 2024).
  - In the code, this data is read and filtered in:
    - `Application/train_processing/01_filter_train_data_switzerland.R`
    - `Application/train_processing/02_prep_train_data_switzerland.R`

- **Precipitation**  
  Hourly precipitation data from:
  - MeteoSchweiz (MeteoSwiss): https://www.meteoswiss.admin.ch  
    (Swiss Federal Office of Meteorology and Climatology), downloaded on 14.08.2024.
  - Obtain hourly precipitation measurements for the Zurich station (e.g. station code `REH`) from the MeteoSwiss data portal.  
  - Store these data files under `Application/data/weather/`.  
  - Matching of train delays to hourly precipitation is performed in:
    - `Application/train_processing/03_match_train_weather.R`

The corresponding R code is located in `Application/train_processing/` and `Application/`.  
The final merged dataset is written as `Application/data/data_combined_train_weather.csv` (not included in the repository).

---

### 2. Precipitation and river discharges (Bavaria, Germany)

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
       - Change of the test statistic in `causal_tail_permutation_test` to test two-sided.

