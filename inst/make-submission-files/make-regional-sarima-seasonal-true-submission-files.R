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
library(nnet)
library(lme4)

registerDoMC(cores=2)
seasonal_difference <- TRUE
delay_adjustment_list <- c("M5")

region_str_array_eval <- c("National",paste0(1:10))
region_str_true <- c("nat",paste0("hhs",1:10))



method <- paste0("sarima_seasonal_difference_", seasonal_difference)
submissions_save_path <- paste0("inst/submissions/region-", method)
data <-readRDS("./data/flu_data_with_backfill_edit.rds")
lag_df <- read.csv("./data/lag_df")
lag_df0 <- read.csv("./data/lag_df0")
subset_lag_df <- lag_df0[lag_df0$week < 201540,]

loess_fit <- nnet(X0~Incidence +season_week + Region, subset_lag_df, size=12, maxit=500, linout=T, decay=0.01)
lm_fit <- lm(X0~Incidence +season_week + Region, data=subset_lag_df)
lm_fit_hierarchical <- lmer(X0~Incidence +season_week + (1|Region), data=subset_lag_df)

fully_observed_data <- as.data.frame(readRDS("./data/fully_observed_data_formatted.rds"))


for (analysis_time_season in c("2017/2018")){
  for (delay_adjustment in delay_adjustment_list){
    if(analysis_time_season == "2017/2018"){
      end_week <- 20
    }else{
      end_week <- 20
    }
    for (test_week_formatted in c(seq(40,52),seq(1,end_week)))   {
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
  
      } else if (delay_adjustment == "M3"){
        
          for (test_region_m3 in unique(current_observed_data$region)){
            if (test_week_formatted == 40){
              
            } else if (test_week_formatted == "01"){
              lag_before_this_one <- data[data$region == test_region_m3 &data$epiweek == paste0(as.numeric(paste0(as.numeric(test_season_formatted)-1,52))) & data$lag==0,]$weighted_ili/data[data$region == test_region_m3 &data$epiweek == paste0(as.numeric(paste0(as.numeric(test_season_formatted)-1,52))) &
                                                                                                                                                                                                data$lag==1,]$weighted_ili
            }else if (test_week_formatted  <=10){
              lag_before_this_one <- data[data$region == test_region_m3 & data$epiweek == paste0(as.numeric(paste0(test_season_formatted,paste0("0",as.numeric(test_week_formatted)-1)))) &
                                            data$lag==0,]$weighted_ili/data[data$region == test_region_m3 & data$epiweek == paste0(as.numeric(paste0(test_season_formatted,paste0("0",as.numeric(test_week_formatted)-1)))) &
                                                                              data$lag==1,]$weighted_ili
            }else{
              
            
              lag_before_this_one <- data[data$region == test_region_m3 &data$epiweek == paste0(as.numeric(paste0(test_season_formatted,as.numeric(test_week_formatted)-1))) &
                                            data$lag==0,]$weighted_ili/data[data$region == test_region_m3 &data$epiweek == paste0(as.numeric(paste0(test_season_formatted,as.numeric(test_week_formatted)-1))) &
                                                                              data$lag==1,]$weighted_ili  
            }
            if (test_week_formatted != 40){
               tmp <- current_observed_data[current_observed_data$region == test_region_m3,]$weighted_ili
               tmp[length(tmp)] <-tmp[length(tmp)]/lag_before_this_one
              current_observed_data[current_observed_data$region == test_region_m3,]$weighted_ili <-tmp
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
        
      }else if (delay_adjustment == "TRUE"){
        current_observed_data <- fully_observed_data[fully_observed_data$epiweek <= paste0(test_season_formatted,test_week_formatted),]
        
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
              prop_estimate <- median(lag_df[ lag_df$week < 201540 , paste0("X",current_lag)],na.rm = T)
              current_observed_data[current_observed_data$region == test_region & current_observed_data$epiweek == paste0(test_season_formatted,lag_itr),]$wili <-
                current_observed_data[current_observed_data$region == test_region &current_observed_data$epiweek == paste0(test_season_formatted,lag_itr),]$wili/prop_estimate
            }
          }
        } else{
          for (lag_itr in seq(40,52)){
            current_lag <- 52 -lag_itr
            for (test_region in unique(data$region)){
              prop_estimate <- median(lag_df[lag_df$week < 201540 , paste0("X",current_lag)],na.rm = T)
              current_observed_data[current_observed_data$region == test_region & current_observed_data$epiweek == paste0(test_season_formatted,lag_itr),]$wili <-
                current_observed_data[current_observed_data$region == test_region &current_observed_data$epiweek == paste0(test_season_formatted,lag_itr),]$wili/prop_estimate
            }
          }
          for (lag_itr in seq(1,as.numeric(test_week_formatted))){
            current_lag <- as.numeric(test_week_formatted) -lag_itr
            for (test_region in unique(data$region)){
              prop_estimate <- median(lag_df[ lag_df$week < 201540 , paste0("X",current_lag)],na.rm = T)
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
        
      } else if (delay_adjustment=="M4"){
            for (test_region in unique(data$region)){
              test_dat <- data.frame(Region=test_region,season_week = as.numeric(test_week_formatted),Incidence=current_observed_data[current_observed_data$region == test_region & current_observed_data$epiweek == paste0(test_season_formatted,test_week_formatted),]$weighted_ili)
              prop_estimate <-predict(loess_fit,newdata=test_dat)
              current_observed_data[current_observed_data$region == test_region & current_observed_data$epiweek == paste0(test_season_formatted,test_week_formatted),]$wili <-                current_observed_data[current_observed_data$region == test_region &current_observed_data$epiweek == paste0(test_season_formatted,test_week_formatted),]$wili/prop_estimate
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
        
      }else if (delay_adjustment=="M5"){
        for (test_region in unique(data$region)){
          test_dat <- data.frame(Region=test_region,season_week = as.numeric(test_week_formatted),Incidence=current_observed_data[current_observed_data$region == test_region & current_observed_data$epiweek == paste0(test_season_formatted,test_week_formatted),]$weighted_ili)
          prop_estimate <-predict(lm_fit,newdata=test_dat)
          current_observed_data[current_observed_data$region == test_region & current_observed_data$epiweek == paste0(test_season_formatted,test_week_formatted),]$wili <-                current_observed_data[current_observed_data$region == test_region &current_observed_data$epiweek == paste0(test_season_formatted,test_week_formatted),]$wili/prop_estimate
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
        
      }else if (delay_adjustment=="M6"){
        for (test_region in unique(data$region)){
          test_dat <- data.frame(Region=test_region,season_week = as.numeric(test_week_formatted),Incidence=current_observed_data[current_observed_data$region == test_region & current_observed_data$epiweek == paste0(test_season_formatted,test_week_formatted),]$weighted_ili)
          prop_estimate <-predict(lm_fit_hierarchical,newdata=test_dat)
          current_observed_data[current_observed_data$region == test_region & current_observed_data$epiweek == paste0(test_season_formatted,test_week_formatted),]$wili <-                current_observed_data[current_observed_data$region == test_region &current_observed_data$epiweek == paste0(test_season_formatted,test_week_formatted),]$wili/prop_estimate
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
