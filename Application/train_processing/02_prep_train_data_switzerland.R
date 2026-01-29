# Match train data (Switzerland) with weather data --------------------------

library(dplyr)
library(tidyr)
library(readr)
library(lubridate)  # for year(), round_date()
library(ggplot2)

# Set working directory to GitHub repo root ---------------------------------
setwd()

# Define paths --------------------------------------------------------------
data_train_path   <- "Application/data/train/"
data_path         <- "Application/data/"
data_weather_path <- "Application/data/weather/"

# Load combined train data (all years, line already filtered) ---------------
data_gesamt <- read_csv(paste0(data_train_path, "data_gesamt.csv"))

# Drop columns that are not needed ------------------------------------------
drop_cols <- c("BETREIBER_ID", "...1", "...2", "...3",
               "FAHRT_BEZEICHNER", "BETREIBER_ABK", "BETREIBER_NAME",
               "PRODUKT_ID", "LINIEN_TEXT")
data_gesamt[, drop_cols] <- NULL

# Create trip identifier per operating day × line ID ------------------------
data_gesamt$LINIEN_ID <- as.factor(data_gesamt$LINIEN_ID)
data_gesamt$key       <- as.factor(paste0(data_gesamt$BETRIEBSTAG,
                                          data_gesamt$LINIEN_ID))

# Remove duplicate rows -----------------------------------------------------
df <- data_gesamt[!duplicated(data_gesamt), ]
nrow(df)
nrow(data_gesamt)

# List of stations on the IC1 corridor between St. Gallen and Genève-Airport
station_list <- c(
  "Bern", "Zürich HB", "Zürich Oerlikon", "Zürich Flughafen",
  "Winterthur", "Wil SG", "Uzwil", "Flawil",
  "Gossau SG", "St. Gallen", "Fribourg/Freiburg",
  "Lausanne", "Genève", "Genève-Aéroport"
)

unique(df$HALTESTELLEN_NAME)

# Group data by key and collect visited stations per trip -------------------
grouped_stations <- df %>%
  group_by(key) %>%
  summarise(names_list = list(HALTESTELLEN_NAME), .groups = "drop")
grouped_stations$names_list

# Helper: check if all stations in v1 are contained in v2[[1]] --------------
is.contained <- function(v1, v2) {
  x <- TRUE
  for (i in v1) {
    if (!(i %in% v2[[1]])) {
      x <- FALSE
    }
  }
  return(x)
}

# Get keys (trips) that include all stations in station_list ----------------
keys_complete <- c()
for (r in 1:nrow(grouped_stations)) {
  print(grouped_stations$key[r])
  if (is.contained(station_list, grouped_stations$names_list[r])) {
    keys_complete <- append(keys_complete, grouped_stations$key[r])
  }
}

length(keys_complete)
write.csv(keys_complete, paste0(data_path, "keys_complete.csv"),
          row.names = FALSE)

nrow(grouped_stations)

# Keep only complete trips (all stations present) ---------------------------
data_comp <- df[df$key %in% keys_complete, ]
nrow(data_comp)

# Inspect trips with non-real-time arrivals/departures at any station -------
data_comp[data_comp$AB_PROGNOSE_STATUS != "REAL" &
            data_comp$AN_PROGNOSE_STATUS != "REAL", ]

# Split trips into two directions:
# r1: St. Gallen -> Genève-Aéroport
# r2: Genève-Aéroport -> St. Gallen
keys_r1   <- c()
keys_r2   <- c()
keys_rest <- c()

for (key in unique(data_comp$key)) {
  tmp <- data_comp[data_comp$key == key, ]
  tmp <- tmp[order(tmp$ABFAHRTSZEIT), ]
  if (tmp$HALTESTELLEN_NAME[1] == "St. Gallen") {
    keys_r1 <- append(keys_r1, key)
  } else if (tmp$HALTESTELLEN_NAME[1] == "Genève-Aéroport") {
    keys_r2 <- append(keys_r2, key)
  } else {
    keys_rest <- append(keys_rest, key)
  }
}

length(keys_r1)
length(keys_r2)
length(keys_rest)

data_r1 <- data_comp[data_comp$key %in% keys_r1, ]
data_r2 <- data_comp[data_comp$key %in% keys_r2, ]

# Check if first station in r1 has valid real-time departure ----------------
keys_unreal_r1 <- c()
for (key in unique(data_r1$key)) {
  tmp <- data_r1[data_r1$key == key, ]
  # Train is not cancelled at St. Gallen?
  if (tmp$FAELLT_AUS_TF[tmp$HALTESTELLEN_NAME == "St. Gallen"] == FALSE) {
    # Departure prognosis missing or not REAL?
    if (is.na(tmp$AB_PROGNOSE_STATUS[tmp$HALTESTELLEN_NAME == "St. Gallen"])) {
      keys_unreal_r1 <- append(keys_unreal_r1, key)
      print(key)
    } else {
      if (tmp$AB_PROGNOSE_STATUS[tmp$HALTESTELLEN_NAME == "St. Gallen"] != "REAL") {
        keys_unreal_r1 <- append(keys_unreal_r1, key)
        print(key)
      }
    }
  }
}
length(keys_unreal_r1)

# Check if first station in r2 has valid real-time departure ----------------
keys_unreal_r2 <- c()
for (key in unique(data_r2$key)) {
  tmp <- data_r2[data_r2$key == key, ]
  # Train is not cancelled at Genève-Aéroport
  if (tmp$FAELLT_AUS_TF[tmp$HALTESTELLEN_NAME == "Genève-Aéroport"][1] == FALSE) {
    # Departure prognosis missing
    if (is.na(tmp$AB_PROGNOSE_STATUS[tmp$HALTESTELLEN_NAME == "Genève-Aéroport"])) {
      keys_unreal_r2 <- append(keys_unreal_r2, key)
      print(key)
    } else {
      if (tmp$AB_PROGNOSE_STATUS[tmp$HALTESTELLEN_NAME == "Genève-Aéroport"] != "REAL") {
        keys_unreal_r2 <- append(keys_unreal_r2, key)
        print(key)
      }
    }
  }
}
length(keys_unreal_r2)
data_r2[data_r2$key == keys_unreal_r2[2], ]

# Filter out trips with problematic first station prognosis -----------------
data_r1_1 <- subset(data_r1, !(key %in% keys_unreal_r1))
nrow(data_r1_1); nrow(data_r1)

data_r2_1 <- subset(data_r2, !(key %in% keys_unreal_r2))
nrow(data_r2_1); nrow(data_r2)

# Check if final station data is available  ---------------------------------
# For r1: final station is Genève-Aéroport
keys_unreal_r1_1 <- c()
for (key in unique(data_r1_1$key)) {
  tmp <- data_r1_1[data_r1_1$key == key, ]
  if (tmp$FAELLT_AUS_TF[tmp$HALTESTELLEN_NAME == "Genève-Aéroport"][1] == FALSE) {
    if (is.na(tmp$AN_PROGNOSE_STATUS[tmp$HALTESTELLEN_NAME == "Genève-Aéroport"][1])) {
      keys_unreal_r1_1 <- append(keys_unreal_r1_1, key)
      print(key)
    } else {
      if (tmp$AN_PROGNOSE_STATUS[tmp$HALTESTELLEN_NAME == "Genève-Aéroport"][1] != "REAL") {
        keys_unreal_r1_1 <- append(keys_unreal_r1_1, key)
        print(key)
      }
    }
  }
}
length(keys_unreal_r1_1)

# For r2: final station is St. Gallen
keys_unreal_r2_1 <- c()
for (key in unique(data_r2_1$key)) {
  tmp <- data_r2_1[data_r2_1$key == key, ]
  if (tmp$FAELLT_AUS_TF[tmp$HALTESTELLEN_NAME == "St. Gallen"] == FALSE) {
    if (is.na(tmp$AN_PROGNOSE_STATUS[tmp$HALTESTELLEN_NAME == "St. Gallen"][1])) {
      keys_unreal_r2_1 <- append(keys_unreal_r2_1, key)
      print(key)
    } else {
      if (tmp$AN_PROGNOSE_STATUS[tmp$HALTESTELLEN_NAME == "St. Gallen"][1] != "REAL") {
        keys_unreal_r2_1 <- append(keys_unreal_r2_1, key)
        print(key)
      }
    }
  }
}
length(keys_unreal_r2_1)

# Filter again, now also on final station data ------------------------------
data_r1_2 <- subset(data_r1_1, !(key %in% keys_unreal_r1_1))
nrow(data_r1_1); nrow(data_r1_2)

data_r2_2 <- subset(data_r2_1, !(key %in% keys_unreal_r2_1))
nrow(data_r2_1); nrow(data_r2_2)

######################## Calculate Delays ###################################

# Define start and end station for delay computation ------------------------
start_st <- "Zürich HB"      # can also be "St. Gallen" in other variants
end_st   <- "Bern"           # can also be "Genève-Aéroport"

# r1: direction ZB (Zürich -> Bern in your notation) ------------------------
delay_table_r1 <- data.frame(
  key      = unique(data_r1_2$key),
  time_ab  = "01.01.2000 00:00:00",
  time_an  = "01.01.2000 00:00:00",
  direction = "ZB"
)
delay_table_r1$delay <- 0

for (r in 1:nrow(delay_table_r1)) {
  key <- delay_table_r1$key[r]
  tmp <- data_r1_2[data_r1_2$key == key, ]
  
  # If train cancelled at start or end station, code as large delay (1000) --
  if (tmp$FAELLT_AUS_TF[tmp$HALTESTELLEN_NAME == end_st][1]   == "TRUE" ||
      tmp$FAELLT_AUS_TF[tmp$HALTESTELLEN_NAME == start_st][1] == "TRUE") {
    delay_table_r1$delay[r] <- 1000
  } else {
    # Planned departure/arrival times --------------------------------------
    t_an <- tmp$ANKUNFTSZEIT[tmp$HALTESTELLEN_NAME == end_st][1]
    t_ab <- tmp$ABFAHRTSZEIT[tmp$HALTESTELLEN_NAME == start_st][1]
    t_an <- strptime(paste0(t_an, ":00"), "%d.%m.%Y %H:%M:%S")
    t_ab <- strptime(paste0(t_ab, ":00"), "%d.%m.%Y %H:%M:%S")
    diff_exp <- as.numeric(difftime(t_an, t_ab, units = "mins"))
    
    # Real-time departure/arrival from prognosis ---------------------------
    r_an <- tmp$AN_PROGNOSE[tmp$HALTESTELLEN_NAME == end_st][1]
    r_ab <- tmp$AB_PROGNOSE[tmp$HALTESTELLEN_NAME == start_st][1]
    delay_table_r1$time_an[r] <- r_an
    delay_table_r1$time_ab[r] <- r_ab
    r_an <- strptime(paste0(r_an, ":00"), "%d.%m.%Y %H:%M:%S")
    r_ab <- strptime(paste0(r_ab, ":00"), "%d.%m.%Y %H:%M:%S")
    diff_real <- as.numeric(difftime(r_an, r_ab, units = "mins"))
    
    # Delay = actual travel time - planned travel time ---------------------
    delay_table_r1$delay[r] <- diff_real - diff_exp
  }
}

hist(delay_table_r1$delay, breaks = 100)
anteil_ausfall_1 <- nrow(delay_table_r1[delay_table_r1$delay == 1000, ]) /
  nrow(delay_table_r1); anteil_ausfall_1

quantile(delay_table_r1$delay[delay_table_r1$delay != 1000], 0.999)
hist(delay_table_r1$delay[delay_table_r1$delay != 1000], breaks = 100)
max(delay_table_r1$delay[delay_table_r1$delay != 1000])

# r2: opposite direction (BZ) -----------------------------------------------
delay_table_r2 <- data.frame(
  key       = unique(data_r2_2$key),
  time_an   = "01.01.2000 00:00:00",
  time_ab   = "01.01.2000 00:00:00",
  direction = "BZ"
)
delay_table_r2$delay <- 0

for (r in 1:nrow(delay_table_r2)) {
  key <- delay_table_r2$key[r]
  tmp <- data_r2_2[data_r2_2$key == key, ]
  
  # If train cancelled at start or end station -----------------------------
  if (tmp$FAELLT_AUS_TF[tmp$HALTESTELLEN_NAME == start_st][1] == "TRUE" ||
      tmp$FAELLT_AUS_TF[tmp$HALTESTELLEN_NAME == end_st][1]   == "TRUE") {
    delay_table_r2$delay[r] <- 1000
  } else {
    # Planned times: now start_st is arrival, end_st is departure ----------
    t_an <- tmp$ANKUNFTSZEIT[tmp$HALTESTELLEN_NAME == start_st][1]
    t_ab <- tmp$ABFAHRTSZEIT[tmp$HALTESTELLEN_NAME == end_st][1]
    t_an <- strptime(paste0(t_an, ":00"), "%d.%m.%Y %H:%M:%S")
    t_ab <- strptime(paste0(t_ab, ":00"), "%d.%m.%Y %H:%M:%S")
    diff_exp <- as.numeric(difftime(t_an, t_ab, units = "mins"))
    
    r_an <- tmp$AN_PROGNOSE[tmp$HALTESTELLEN_NAME == start_st][1]
    r_ab <- tmp$AB_PROGNOSE[tmp$HALTESTELLEN_NAME == end_st][1]
    delay_table_r2$time_an[r] <- r_an[1]
    delay_table_r2$time_ab[r] <- r_ab[1]
    r_an <- strptime(paste0(r_an, ":00"), "%d.%m.%Y %H:%M:%S")
    r_ab <- strptime(paste0(r_ab, ":00"), "%d.%m.%Y %H:%M:%S")
    diff_real <- as.numeric(difftime(r_an, r_ab, units = "mins"))
    
    delay_table_r2$delay[r] <- diff_real - diff_exp
  }
}

hist(delay_table_r2$delay, breaks = 100)
anteil_ausfall_2 <- nrow(delay_table_r2[delay_table_r2$delay == 1000, ]) /
  nrow(delay_table_r2); anteil_ausfall_2

quantile(delay_table_r2$delay[delay_table_r2$delay != 1000], 0.999)
hist(delay_table_r2$delay[delay_table_r2$delay != 1000], breaks = 500, freq = FALSE)
max(delay_table_r2$delay[delay_table_r2$delay != 1000])

########## Combine both directions ##########################################

train_delays      <- rbind(delay_table_r1, delay_table_r2)
train_delays_save <- apply(train_delays, 2, as.character)  # ensure no factors

write.csv(train_delays_save,
          file = paste0(data_train_path, "train_delays.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# Load and combine weather station data (Swiss stations)
# ---------------------------------------------------------------------------

file_list    <- list.files(data_weather_path, pattern = "swiss_stations")
data_weather <- NULL

for (file in file_list) {
  print(file)
  
  if (exists("data_weather")) {
    temp_dataset <- read.csv(paste0(data_weather_path, file))
    data_weather <- rbind(data_weather, temp_dataset)
    rm(temp_dataset)
  }
  if (!exists("data_weather")) {
    data_weather <- read.csv(paste0(data_weather_path, file))
  }
}

summary(data_weather)
head(data_weather)

data_weather$time    <- strptime(data_weather$time, "%Y-%m-%d %H:%M:%S")
data_weather$station <- as.factor(data_weather$nat_abbr)

# Extract Bern (BER) and Zürich (REH) stations ------------------------------
data_weather_bern <- subset(data_weather, station == "BER")
data_weather_zh   <- subset(data_weather, station == "REH")

# ---------------------------------------------------------------------------
# Match train delays with hourly precipitation (Bern & Zurich)
# ---------------------------------------------------------------------------

data_trains <- read.csv(paste0(data_train_path, "train_delays.csv"))

data_trains$time_ab <- strptime(data_trains$time_ab, "%d.%m.%Y %H:%M:%S")
data_trains$time_an <- strptime(data_trains$time_an, "%d.%m.%Y %H:%M:%S")

# Round times to full hours for merging with hourly weather data -----------
data_trains$time_ab_r <- round_date(data_trains$time_ab, unit = "hours")
data_trains$time_an_r <- round_date(data_trains$time_an, unit = "hours")

head(data_trains)

# Initialize merged data frame with precipitation columns -------------------
data_all <- data_trains

data_all$rain       <- NA
data_all$rain_3     <- NA
data_all$rain_6     <- NA
data_all$rain_12    <- NA
data_all$rain_24    <- NA
data_all$rain_zh    <- NA
data_all$rain_3_zh  <- NA
data_all$rain_6_zh  <- NA
data_all$rain_12_zh <- NA
data_all$rain_24_zh <- NA

# Remove cancelled trips (delay == 1000) ------------------------------------
data_all <- data_all[data_all$delay < 1000, ]

# Keep weather data from 2020 onwards ---------------------------------------
data_weather_bern <- data_weather_bern[year(data_weather_bern$time) > 2019, ]
data_weather_zh   <- data_weather_zh[year(data_weather_zh$time)   > 2019, ]

# Loop over all train trips and attach precipitation summaries --------------
for (i in 1:nrow(data_all)) {
  print(i)
  
  # Determine which timestamp to use depending on direction -----------------
  if (data_all$direction[i] == "ZB") {
    time   <- data_all$time_an_r[i]
    time_1 <- data_all$time_ab_r[i]  # currently unused
  } else {
    time   <- data_all$time_ab_r[i]
    time_1 <- data_all$time_an_r[i]
  }
  
  # Find matching weather row for Bern (arrival/departure time) ------------
  row <- which(data_weather_bern$time == time)
  
  # 1h, 3h, 6h, 12h, 24h precipitation sums in Bern ------------------------
  data_all$rain[i]    <- data_weather_bern$rre150h0[row]
  data_all$rain_3[i]  <- sum(data_weather_bern$rre150h0[(row - 2):row])
  data_all$rain_6[i]  <- sum(data_weather_bern$rre150h0[(row - 5):row])
  data_all$rain_12[i] <- sum(data_weather_bern$rre150h0[(row - 11):row])
  data_all$rain_24[i] <- sum(data_weather_bern$rre150h0[(row - 23):row])
  
  # Same for Zurich (ZH) ----------------------------------------------------
  data_all$rain_zh[i]    <- data_weather_zh$rre150h0[row]
  data_all$rain_3_zh[i]  <- sum(data_weather_zh$rre150h0[(row - 2):row])
  data_all$rain_6_zh[i]  <- sum(data_weather_zh$rre150h0[(row - 5):row])
  data_all$rain_12_zh[i] <- sum(data_weather_zh$rre150h0[(row - 11):row])
  data_all$rain_24_zh[i] <- sum(data_weather_zh$rre150h0[(row - 23):row])
}

head(data_all)
summary(data_all)

# Save final combined train–weather dataset ---------------------------------
write.csv(data_all, paste0(data_path, "data_combined_train_weather.csv"),
          row.names = FALSE)