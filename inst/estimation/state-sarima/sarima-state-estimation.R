## state-level SARIMA fits
## 5 Oct 2016: started for regional level
## 3 Nov 2017: adapted for state-level
## 9 Oct 2018: updated for 2018/2019 season, including sarimaTD migration
## Nicholas Reich, Casey Gibson

library(plyr);library(dplyr)
library(tidyr)
library(lubridate)
library(cdcFlu20182019)
library(forecast)
library(sarimaTD)

## do all regional fits in parallel
library(doMC)
registerDoMC(25)

data(state_flu_data)
state_names <- unique(state_flu$region)
state_names <- state_names[-which(state_names %in% c("Florida"))]

#state_names <- c("Alabama","Alaska")

## Florida and Louisiana: drop completely
## PR: starts in end of 2013
## VI: starts in end of 2015

region_seasons <- expand.grid(
  region = state_names,
  #first_test_season = paste0(2011:2018, "/", 2012:2019), 
  first_test_season = "2010/2011", # uncomment this line if you only want one season as test
  stringsAsFactors = FALSE
)

## get fits with seasonal differencing before call to auto.arima
## where to store SARIMA fits
path <- paste0("inst/estimation/state-sarima/fits-seasonal-differencing/")

#foreach(i = seq_len(nrow(region_seasons))) %dopar% {
for(i in 1:nrow(region_seasons)){
    message(paste(Sys.time(), "starting", region_seasons$region[i]))
  fit_region_sarima(data = state_flu,
    region = region_seasons$region[i],
    first_test_season = region_seasons$first_test_season[i],
    seasonal_difference = TRUE,
    transformation = "box-cox", 
    prediction_target_var = "unweighted_ili",
    path = path)
}
