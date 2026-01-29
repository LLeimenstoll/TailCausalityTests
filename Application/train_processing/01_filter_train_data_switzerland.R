# Script to load raw train data and keep only line IC1 ----------------------

library(dplyr)   # data manipulation (not heavily used here, but available)
library(tidyr)   # tidying tools (not used in this script)
library(readr)   # read_csv() for faster CSV I/O

# Set working directory to GitHub repo root ---------------------------------
setwd()

# NOTE: data_path must be defined somewhere before this script, e.g.:
data_train_path <- "Application/data/train/"

# Get list of subdirectories containing "ist-daten" -------------------------
file_list <- list.files(
  data_train_path,
  pattern = "ist-daten"
)

# Restrict to specific subset of folders for each year (here elements 25 to 36) ----------
file_list <- file_list[25:36]

# ---------------------------------------------------------------------------
# Part 1: Filter each raw file to retain only IC1 line
# ---------------------------------------------------------------------------

for (file in file_list) {
  print(file)
  
  # List all data files in this subfolder whose names contain "istdaten" ----
  data_list <- list.files(
    paste0(data_train_path, file, "/"),
    pattern = "istdaten"
  )
  
  for (data in data_list) {
    print(data)
    
    # Read raw CSV with ';' separator (original SBB format) -----------------
    train_data <- read.csv(paste0(data_train_path, file, "/", data), sep = ";")
    
    # Keep only rows corresponding to line "IC1" ----------------------------
    train_data_new <- train_data[train_data$LINIEN_TEXT == "IC1", ]
    
    # Overwrite original file with filtered IC1 data ------------------------
    write.csv(train_data_new,
              file = paste0(data_train_path, file, "/", data),
              row.names = FALSE)
  }
}

# ---------------------------------------------------------------------------
# Part 2: Create merged dataset for one year (across all folders)
# ---------------------------------------------------------------------------

data_2022 <- NULL  # initialize combined data object

for (file in file_list) {
  print(file)
  
  # List all (already filtered) files in this subfolder ---------------------
  data_list <- list.files(
    paste0(data_train_path, file, "/"),
    pattern = "istdaten"
  )
  
  for (data in data_list) {
    if (exists("data_2022")) {
      # Append to existing combined dataset --------------------------------
      temp_dataset <- read_csv(paste0(data_train_path, file, "/", data))
      data_2022 <- rbind(data_2022, temp_dataset)
      rm(temp_dataset)
    }
    
    # If merged dataset does not exist yet, initialize it -------------------
    if (!exists("data_2022")) {
      data_2022 <- read_csv(paste0(data_train_path, file, "/", data))
    }
  }
}

summary(data_2022)

# Save combined yearly dataset (here: e.g. for 2022) ------------------------
write.csv(data_2022,
          file = paste0(data_train_path, "data_2022.csv"),
          row.names = FALSE)

tail(data_2022)

# Check proportion of a specific status (AN_PROGNOSE_STATUS) ----------------
s <- summary(as.factor(data_2022$AN_PROGNOSE_STATUS))
s / nrow(data_2022)

# ---------------------------------------------------------------------------
# Part 3: Bind all years together (2021–2024)
# ---------------------------------------------------------------------------

data_2021 <- read_csv(paste0(data_train_path, "data_2021.csv"))
data_2022 <- read_csv(paste0(data_train_path, "data_2022.csv"))
data_2023 <- read_csv(paste0(data_train_path, "data_2023.csv"))
data_2024 <- read_csv(paste0(data_train_path, "data_2024.csv"))

# Stack rows from all years into one big dataset ----------------------------
data_gesamt <- rbind(data_2021, data_2022, data_2023, data_2024)

# Save final full-period dataset --------------------------------------------
write.csv(data_gesamt,
          file = paste0(data_train_path, "data_gesamt.csv"),
          row.names = FALSE)