library(EnvStats)
library(ggplot2)

get_model_prob <- function(season,week,target,model,region) {
  model_csv <- read.csv(paste0("./inst/submissions/region-sarima_seasonal_difference_TRUE/EW",week,"-",season,"-ReichLab_sarima_seasonal_difference_TRUE-delay-",model,".csv"))
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
  point <- model_csv[model_csv$Location ==region & model_csv$Target == target & model_csv$Type == "Point"  , ]$Value
  
  return (c(prob,prob_l,prob_r,truth,point))
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
  
  onset_week <- get_onset_week(round(trajectory,1),baseline = baseline,3,first_season_week = 31, get_num_MMWR_weeks_in_first_season_year(analysis_time_season))
  truth <- season_week_to_year_week(onset_week,first_season_week = 31,weeks_in_first_season_year =get_num_MMWR_weeks_in_first_season_year(analysis_time_season) )
  
  if (truth == 52){
    truth_l <- 51
    truth_r <- 1
  } else if(truth == 1){
    truth_l <- 52
    truth_r <- 2
  } else if (truth == 40){
    truth_l <- 40
    truth_r <- 41
  } else if (truth == 20){
    truth_l <- 19
    truth_r <- 20
  }else{
    truth_l <- truth -1
    truth_r <- truth + 1
  }
  
  
  
  prob <- model_csv[model_csv$Location ==region & model_csv$Target == target & model_csv$Type == "Bin" & model_csv$Bin_start_incl ==  truth , ]$Value
  prob_l <- model_csv[model_csv$Location ==region & model_csv$Target == target & model_csv$Type == "Bin" & model_csv$Bin_start_incl ==  truth_l , ]$Value
  prob_r <- model_csv[model_csv$Location ==region & model_csv$Target == target & model_csv$Type == "Bin" & model_csv$Bin_start_incl ==  truth_r , ]$Value
  point <- model_csv[model_csv$Location ==region & model_csv$Target == target & model_csv$Type == "Point"  , ]$Value
  
  return (c(prob,prob_l,prob_r,truth,point))
  
}


gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}


fully_observed_data <- data.frame(readRDS("./data/fully_observed_data_formatted.rds"))
fully_observed_data <- fully_observed_data[order(fully_observed_data$epiweek),]
data <- readRDS("./data/flu_data_with_backfill_edit.rds")
lag_df0 <- read.csv("lag_df0")
lag_df <- read.csv("data/lag_df")
region_str_array_eval <- c("National",paste0("Region ",1:10))
region_str_data_set <- c("US National","HHS Region 1",  "HHS Region 2",  "HHS Region 3",  "HHS Region 4" ,"HHS Region 5",  "HHS Region 6",  "HHS Region 7",  "HHS Region 8",  "HHS Region 9","HHS Region 10")

region_str_true <- c("nat",paste0("hhs",1:10))


targets <- c("Season onset","1 wk ahead","2 wk ahead","3 wk ahead","4 wk ahead")
test_seasons <-c("2015")
test_models <- c("FSMOOTHED","TRUE","NONE")
test_regions <- region_str_data_set

result_df <- matrix(NA,ncol= 10)

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

colnames(result_df) <- c("Season","week","target","model","region","p_c","p_l","p_r","truth","point")

result_df <- as.data.frame(result_df)

result_df$p_c <- as.numeric(as.character(result_df$p_c))
result_df$p_r <- as.numeric(as.character(result_df$p_r))
result_df$p_l <- as.numeric(as.character(result_df$p_l))
result_df$week <- as.numeric(as.character(result_df$week))
result_df$total_prob <- result_df$p_c +result_df$p_r +result_df$p_l

result_df$total_log_prob <- pmax(log(result_df$total_prob),-10)
result_df$ll <- rep(NA, nrow(result_df))

mean(result_df[result_df$model=="TRUE"  ,]$total_log_prob,na.rm = TRUE)
mean(result_df[result_df$model=="NONE",]$total_log_prob,na.rm = TRUE)
#mean(result_df[result_df$model=="M1",]$total_log_prob,na.rm = TRUE)

library(ggplot2)
library(dplyr)
ggplot(result_df[result_df$Season == "2013"  & result_df$target == "1 wk ahead",],aes(x=week,y=total_log_prob,col=region)) + geom_point() + facet_grid(model ~.)

ggplot(result_df[result_df$Season == "2013"  & result_df$target == "1 wk ahead",],aes(x=week,y=total_log_prob,col=region)) + geom_point() + facet_grid(model ~.)

current_observed_data <- data %>% group_by(region,epiweek) %>%
  filter(lag == 0)
current_observed_data <- data.frame(current_observed_data[order(current_observed_data$epiweek),])


data_for_plot <- rbind(current_observed_data,fully_observed_data)
  
data_for_plot$f <- c(rep("Currently Observed",nrow(current_observed_data)),rep("Fully Observed",nrow(fully_observed_data)))
lag_0_by_week <- lag_df %>% group_by(season_week,Region) %>% summarize(X0=var(X0,na.rm = T))
lag_0_by_week <- data.frame(lag_0_by_week)
lag_0_by_week$f <- rep("Ratio",nrow(lag_0_by_week))

ggplot(data_for_plot[data_for_plot$region == "Region 7" & data_for_plot$season == "2010/2011",],aes(x=season_week,y=weighted_ili,col=f)) + geom_line() +
  geom_line(data=lag_0_by_week[lag_0_by_week$Region == "Region 7",],aes(x=season_week,y=X0))



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
