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
registerDoMC(4)

#data(flu_data)
flu_data <- as.data.frame(readRDS("data/fully_observed_data_formatted.rds"))
flu_data <- flu_data[order(flu_data$epiweek),]
# for (row in 1:nrow(flu_data)){
#   if (flu_data[row,]$week >20){
#       current_year <- flu_data[row,]$year
#       current_year_1 <- current_year + 1
#       flu_data[row,"season"] <- paste0(current_year,"/",current_year_1)
#   }else{
#     current_year <- flu_data[row,]$year
#     current_year_1 <- current_year - 1
#     flu_data[row,"season"] <- paste0(current_year_1,"/",current_year)
#   }
# }
for (analysis_time_season in c( "2011/2012","2012/2013","2013/2014")){
  

  region_names <- as.character(unique(flu_data$region))
  
  ## Florida and Louisiana: drop completely
  ## PR: starts in end of 2013
  ## VI: starts in end of 2015
  
  region_seasons <- expand.grid(
    region = region_names,
    #first_test_season = paste0(2011:2018, "/", 2012:2019), 
    first_test_season =analysis_time_season , # uncomment this line if you only want one season as test
    stringsAsFactors = FALSE
  )
  
  ## get fits with seasonal differencing before call to auto.arima
  ## where to store SARIMA fits
  path <- paste0("inst/estimation/region-sarima/fits-seasonal-differencing/")
  
  #foreach(i = seq_len(nrow(region_seasons))) %dopar% {
  for(i in 1:nrow(region_seasons)){
      message(paste(Sys.time(), "starting", region_seasons$region[i]))
    fit_region_sarima(data = flu_data,
      region = region_seasons$region[i],
      first_test_season = region_seasons$first_test_season[i],
      seasonal_difference = TRUE,
      transformation = "box-cox", 
      prediction_target_var = "weighted_ili",
      path = path)
  }
  
  
  
  path <- paste0("inst/estimation/region-sarima/fits-no-seasonal-differencing/")
  
  #foreach(i = seq_len(nrow(region_seasons))) %dopar% {
  for(i in 1:nrow(region_seasons)){
    message(paste(Sys.time(), "starting", region_seasons$region[i]))
    fit_region_sarima(data = flu_data,
                      region = region_seasons$region[i],
                      first_test_season = region_seasons$first_test_season[i],
                      seasonal_difference = FALSE,
                      transformation = "box-cox", 
                      prediction_target_var = "weighted_ili",
                      path = path)
  }
}
