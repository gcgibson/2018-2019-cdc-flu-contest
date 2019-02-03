## This code is based on inst/data-processing/create-clean-flu-data.R

library(plyr)
library(dplyr)
library(tidyr)
library(lubridate)
library(MMWRweek)
library(cdcfluview)
library(cdcFlu20182019)
library(forecast)
#library(kcde)
#library(copula)
library(FluSight)
library(gridExtra)
library(foreach)
library(doMC)
registerDoMC(cores=4)
seasonal_difference <- TRUE
delay_adjustment_list <- c("M1")

region_str_array_eval <- c("National",paste0(1:10))
region_str_true <- c("nat",paste0("hhs",1:10))



method <- paste0("sarima_seasonal_difference_", seasonal_difference)
submissions_save_path <- paste0("inst/submissions/region-", method)
data <-readRDS("./data/flu_data_with_backfill_edit.rds")
lag_df <- read.csv("./data/lag_df")

for (analysis_time_season in c("2015/2016")){
  for (delay_adjustment in delay_adjustment_list){
    for (test_week_formatted  in c(seq(40,52),seq(1,20))) {
      if (test_week_formatted < 40){
        test_season_formatted <- substr(analysis_time_season,6,9)
      } else{
        test_season_formatted <- substr(analysis_time_season,1,4)
        
      }
      if (test_week_formatted <=9){
        test_week_formatted <- paste0("0",test_week_formatted)
      }
      
      
      current_observed_data_with_lags <- data[ data$issue <= as.numeric(paste0(test_season_formatted,test_week_formatted)),]
      current_observed_data <- current_observed_data_with_lags %>% group_by(region,epiweek) %>%
        filter(lag == max(lag))
      current_observed_data <- current_observed_data[order(current_observed_data$epiweek),]
      if (delay_adjustment == "NONE") {
        current_observed_data <- as.data.frame(current_observed_data)
        simulate_trajectories_sarima_params <- list(
          fits_filepath = paste0("inst/estimation/region-sarima/",
                                 ifelse(seasonal_difference,
                                        "fits-seasonal-differencing",
                                        "fits-no-seasonal-differencing")),
          prediction_target_var = "weighted_ili",
          seasonal_difference = seasonal_difference,
          transformation = "box-cox",
          first_test_season = analysis_time_season,
          do_sampling_lag = FALSE
        )
        
        weeks_in_first_season_year <-
          get_num_MMWR_weeks_in_first_season_year(analysis_time_season)
        
        res <- get_submission_via_trajectory_simulation(
          data = current_observed_data,
          analysis_time_season = analysis_time_season,
          first_analysis_time_season_week = 10, # == week 40 of year
          last_analysis_time_season_week = weeks_in_first_season_year - 11, # at week 41, we do prediction for a horizon of one week ahead
          prediction_target_var = "weighted_ili",
          incidence_bins = data.frame(
            lower = c(0, seq(from = 0.05, to = 12.95, by = 0.1)),
            upper = c(seq(from = 0.05, to = 12.95, by = 0.1), Inf)),
          incidence_bin_names = as.character(seq(from = 0, to = 13, by = 0.1)),
          n_trajectory_sims = 10000,
          #  n_trajectory_sims = 100,
          simulate_trajectories_function = sample_predictive_trajectories_arima_wrapper,
          simulate_trajectories_params = simulate_trajectories_sarima_params,
          regional_switch="Country")
  
      }else if (delay_adjustment=="M1"){
        if (test_week_formatted >=40){
          for (lag_itr in seq(40,test_week_formatted)){
            current_lag <- as.numeric(test_week_formatted) -lag_itr
            for (test_region in unique(data$region)){
              prop_estimate <- mean(lag_df[lag_df$Region == test_region & lag_df$week < 201540 , paste0("X",current_lag)],na.rm = T)
              current_observed_data[current_observed_data$region == test_region & current_observed_data$epiweek == paste0(test_season_formatted,lag_itr),]$wili <-
                current_observed_data[current_observed_data$region == test_region &current_observed_data$epiweek == paste0(test_season_formatted,lag_itr),]$wili/prop_estimate
            }
          }
        } else{
          for (lag_itr in seq(40,52)){
            current_lag <- 52 -lag_itr
            for (test_region in unique(data$region)){
              prop_estimate <- mean(lag_df[lag_df$Region == test_region & lag_df$week < 201540 , paste0("X",current_lag)],na.rm = T)
              current_observed_data[current_observed_data$region == test_region & current_observed_data$epiweek == paste0(test_season_formatted,lag_itr),]$wili <-
                current_observed_data[current_observed_data$region == test_region &current_observed_data$epiweek == paste0(test_season_formatted,lag_itr),]$wili/prop_estimate
            }
          }
          for (lag_itr in seq(1,as.numeric(test_week_formatted))){
            current_lag <- as.numeric(test_week_formatted) -lag_itr
            for (test_region in unique(data$region)){
              prop_estimate <- mean(lag_df[lag_df$Region == test_region & lag_df$week < 201540 , paste0("X",current_lag)],na.rm = T)
              current_observed_data[current_observed_data$region == test_region & current_observed_data$epiweek == paste0(test_season_formatted,lag_itr),]$wili <-
                current_observed_data[current_observed_data$region == test_region &current_observed_data$epiweek == paste0(test_season_formatted,lag_itr),]$wili/prop_estimate
            }
          }
        }
        current_observed_data <- as.data.frame(current_observed_data)
        simulate_trajectories_sarima_params <- list(
          fits_filepath = paste0("inst/estimation/region-sarima/",
                                 ifelse(seasonal_difference,
                                        "fits-seasonal-differencing",
                                        "fits-no-seasonal-differencing")),
          prediction_target_var = "weighted_ili",
          seasonal_difference = seasonal_difference,
          transformation = "box-cox",
          first_test_season = analysis_time_season,
          do_sampling_lag=FALSE
        )
        
        weeks_in_first_season_year <-
          get_num_MMWR_weeks_in_first_season_year(analysis_time_season)
        
        res <- get_submission_via_trajectory_simulation(
          data = current_observed_data,
          analysis_time_season = analysis_time_season,
          first_analysis_time_season_week = 10, # == week 40 of year
          last_analysis_time_season_week = weeks_in_first_season_year - 11, # at week 41, we do prediction for a horizon of one week ahead
          prediction_target_var = "weighted_ili",
          incidence_bins = data.frame(
            lower = c(0, seq(from = 0.05, to = 12.95, by = 0.1)),
            upper = c(seq(from = 0.05, to = 12.95, by = 0.1), Inf)),
          incidence_bin_names = as.character(seq(from = 0, to = 13, by = 0.1)),
          n_trajectory_sims = 10000,
          #  n_trajectory_sims = 100,
          simulate_trajectories_function = sample_predictive_trajectories_arima_wrapper,
          simulate_trajectories_params = simulate_trajectories_sarima_params,
          regional_switch="Country")
        
      } else if (delay_adjustment == "M2"){
        simulate_trajectories_sarima_params <- list(
          fits_filepath = paste0("inst/estimation/region-sarima/",
                                 ifelse(seasonal_difference,
                                        "fits-seasonal-differencing",
                                        "fits-no-seasonal-differencing")),
          prediction_target_var = "weighted_ili",
          seasonal_difference = seasonal_difference,
          transformation = "box-cox",
          first_test_season = analysis_time_season,
          do_sampling_lag=TRUE
        )
        
        weeks_in_first_season_year <-
          get_num_MMWR_weeks_in_first_season_year(analysis_time_season)
        current_observed_data <- as.data.frame(current_observed_data)
        
        res <- get_submission_via_trajectory_simulation(
          data = current_observed_data,
          analysis_time_season = analysis_time_season,
          first_analysis_time_season_week = 10, # == week 40 of year
          last_analysis_time_season_week = weeks_in_first_season_year - 11, # at week 41, we do prediction for a horizon of one week ahead
          prediction_target_var = "weighted_ili",
          incidence_bins = data.frame(
            lower = c(0, seq(from = 0.05, to = 12.95, by = 0.1)),
            upper = c(seq(from = 0.05, to = 12.95, by = 0.1), Inf)),
          incidence_bin_names = as.character(seq(from = 0, to = 13, by = 0.1)),
          n_trajectory_sims = 10000,
          #  n_trajectory_sims = 100,
          simulate_trajectories_function = sample_predictive_trajectories_arima_wrapper,
          simulate_trajectories_params = simulate_trajectories_sarima_params,
          regional_switch="Country")
        
      }
      
      
      res_file <- file.path(submissions_save_path,
        paste0(
          "EW", sprintf("%02d", tail(current_observed_data$week, 1)),
          "-", tail(current_observed_data$year, 1),
          "-ReichLab_", method,"-delay-",delay_adjustment,
          ".csv"))
      write.csv(res,
        file = res_file,
        row.names = FALSE)
    }
  }
}
#(FluSight::verify_entry_file(res_file))

### Plots for sanity

# make_predictions_plots(
#   preds_save_file = res_file,
#   plots_save_file = paste0(
#     submissions_save_path,
#     "/plots/",
#     tail(data$year, 1),
#     "-",
#     sprintf("%02d", tail(data$week, 1)),
#     "-KOT-", method, "-plots.pdf"),
#   data = data
# )
