library(EnvStats)
library(ggplot2)

get_model_prob <- function(season,week,target,model,region) {
  model_csv <- read.csv(paste0("./inst/submissions/region-sarima_seasonal_difference_TRUE/EW",week,"-",top_level_season_formatted,"-ReichLab_sarima_seasonal_difference_TRUE-delay-",model,".csv"))
  if(target == "1 wk ahead" | target == "2 wk ahead" | target == "3 wk ahead" | target == "4 wk ahead"){
    prob_and_truth <- get_model_prob_k_week_ahead(season,week,model,region,model_csv,target)
  } else if (target == "Season onset"){
    prob_and_truth <- get_model_prob_season_onset(season,week,model,region,model_csv,target)
  }
  return (prob_and_truth)
}


get_model_prob_k_week_ahead <- function(season,week,model,region,model_csv,target){
  formatted_region <- region_str_array_eval[match(region,region_str_data_set)]
  truth <- round(fully_observed_data[fully_observed_data$region == formatted_region & fully_observed_data$epiweek ==paste0(season,week),]$weighted_ili,1)
  truth_l <- max(0,truth-.1)
  truth_r <- min(13, truth+.1)
  
  prob <- model_csv[model_csv$Location ==region & model_csv$Target == target & model_csv$Type == "Bin" & model_csv$Bin_start_incl ==  truth , ]$Value
  prob_l <- model_csv[model_csv$Location ==region & model_csv$Target == target & model_csv$Type == "Bin" & model_csv$Bin_start_incl ==  truth_l , ]$Value
  prob_r <- model_csv[model_csv$Location ==region & model_csv$Target == target & model_csv$Type == "Bin" & model_csv$Bin_start_incl ==  truth_r , ]$Value
  
  return (c(prob,prob_l,prob_r,truth))
}


get_model_prob_season_onset <- function(season,week,model,region,model_csv,target){
  if (week >= 40){
    analysis_time_season <- paste0(season,"/",as.numeric(season) + 1)
  } else{
    analysis_time_season <- paste0(as.numeric(season-1),"/",as.numeric(season))
  }
  formatted_region <- region_str_array_eval[match(region,region_str_data_set)]
  trajectory <- fully_observed_data[fully_observed_data$region == formatted_region & fully_observed_data$season == analysis_time_season ,]$weighted_ili
  baseline <- get_onset_baseline(region = formatted_region, season = analysis_time_season)
  
  truth_unformatted <- min(which (trajectory > baseline))
  truth <- fully_observed_data[fully_observed_data$region == formatted_region & fully_observed_data$season == analysis_time_season & fully_observed_data$season_week == truth_unformatted,]$week
  truth_l <- max(1,truth-1)
  truth_r <- min(52, truth+1)
  prob <- model_csv[model_csv$Location ==region & model_csv$Target == target & model_csv$Type == "Bin" & model_csv$Bin_start_incl ==  truth , ]$Value
  prob_l <- model_csv[model_csv$Location ==region & model_csv$Target == target & model_csv$Type == "Bin" & model_csv$Bin_start_incl ==  truth_l , ]$Value
  prob_r <- model_csv[model_csv$Location ==region & model_csv$Target == target & model_csv$Type == "Bin" & model_csv$Bin_start_incl ==  truth_r , ]$Value
  return (c(prob,prob_l,prob_r,truth))
  
}


gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}


fully_observed_data <- readRDS("./data/fully_observed_data_formatted.rds")

lag_df0 <- read.csv("lag_df0")

region_str_array_eval <- c("National",paste0("Region ",1:10))
region_str_data_set <- c("US National","HHS Region 1",  "HHS Region 2",  "HHS Region 3",  "HHS Region 4" ,"HHS Region 5",  "HHS Region 6",  "HHS Region 7",  "HHS Region 8",  "HHS Region 9","HHS Region 10")

region_str_true <- c("nat",paste0("hhs",1:10))


targets <- c("Season onset")#,"2 wk ahead","3 wk ahead","4 wk ahead")
test_seasons <-c("2015","2016","2017")
test_models <- c("TRUE","NONE","M1")
test_regions <- region_str_data_set

result_df <- matrix(NA,ncol= 9)

for (region in test_regions){
  for (model in test_models){ 
    for (target in targets){
      for (top_level_season in test_seasons){
        for (week in c(seq(41,52),paste0("0",1:9),seq(10,20))){
          if (week >= 41){
            top_level_season_formatted <- top_level_season
          } else{
            top_level_season_formatted <- as.numeric(top_level_season) + 1
          }
           model_prob_and_truth <- get_model_prob(top_level_season_formatted,week,target,model,region) 
           result_df <- rbind(result_df,c(top_level_season,week,target,model,region,model_prob_and_truth))
        }
      }   
    }
  }
}

result_df <- result_df[2:nrow(result_df),]

colnames(result_df) <- c("Season","week","target","model","region","p_c","p_l","p_r","truth")

result_df <- as.data.frame(result_df)

result_df$p_c <- as.numeric(as.character(result_df$p_c))
result_df$p_r <- as.numeric(as.character(result_df$p_r))
result_df$p_l <- as.numeric(as.character(result_df$p_l))
result_df$week <- as.numeric(as.character(result_df$week))
result_df$total_prob <- result_df$p_c +result_df$p_r +result_df$p_l

result_df$total_log_prob <- pmax(log(result_df$total_prob),-10)
result_df$ll <- rep(NA, nrow(result_df))

mean(result_df[result_df$model=="TRUE" & result_df$region=="HHS Region 2" ,]$total_log_prob,na.rm = TRUE)
mean(result_df[result_df$model=="NONE",]$total_log_prob,na.rm = TRUE)
mean(result_df[result_df$model=="M1",]$total_log_prob,na.rm = TRUE)




for (region in unique(result_df$region)){
  for (season in unique(result_df$Season)){
    for (week in unique(result_df$week)){
      for (model in unique(result_df$model)){
        for (target in unique(result_df$target)){
          
          formatted_season <- paste0(season,"/",as.numeric(season)+1)
          l_0 <- data[data$region == region_str_array_eval[match(region,region_str_data_set)] & data$season == formatted_season & data$week == week & data$lag == 0,]$weighted_ili
          l_infty <- data[data$region == region_str_array_eval[match(region,region_str_data_set)] & data$season == formatted_season & data$week == week & data$lag == max(data[data$region == region_str_array_eval[match(region,region_str_data_set)] & data$season == formatted_season & data$week == week,]$lag),]$weighted_ili
          result_df[result_df$Season ==season &result_df$target == target & result_df$week == week & result_df$region == region & result_df$model==model,]$ll <-l_0-l_infty 
        }
      }
    }
  }
}
