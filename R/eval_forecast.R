library(EnvStats)
library(ggplot2)

do_eval <- function(){
  get_model_prob_k_week_ahead <- function(season,week,model,region,model_csv,target){
  formatted_region <- region_str_array_eval[match(region,region_str_data_set)]
  if (as.numeric(week) == 52){
    formatted_week_plus_1 <- 1
  } else{
    formatted_week_plus_1 <- as.numeric(week) + 1 
  }
  if (as.numeric(week) < 9 | as.numeric(week) == 52){
    if (as.numeric(week) != 52) {
      truth <- round(fully_observed_data[fully_observed_data$region == formatted_region & fully_observed_data$epiweek ==paste0(season,"0",formatted_week_plus_1),]$weighted_ili,1)
    } else{
      truth <- round(fully_observed_data[fully_observed_data$region == formatted_region & fully_observed_data$epiweek ==paste0(as.numeric(season)+1,"0",formatted_week_plus_1),]$weighted_ili,1)
      
    }
  } else{
    truth <- round(fully_observed_data[fully_observed_data$region == formatted_region & fully_observed_data$epiweek ==paste0(season,formatted_week_plus_1),]$weighted_ili,1)
    
  }
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


get_model_prob_season_peak_week <- function(season,week,model,region,model_csv,target){
  if (week >= 40){
    analysis_time_season <- paste0(season,"/",as.numeric(season) + 1)
  } else{
    analysis_time_season <- paste0(as.numeric(season-1),"/",as.numeric(season))
  }
  formatted_region <- region_str_array_eval[match(region,region_str_data_set)]
  trajectory <- fully_observed_data[fully_observed_data$region == formatted_region & fully_observed_data$season == analysis_time_season ,]$weighted_ili
  
  truth <- season_week_to_year_week(which.max(trajectory),first_season_week = 31,weeks_in_first_season_year =get_num_MMWR_weeks_in_first_season_year(analysis_time_season) )
  
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


get_model_prob_season_peak_percentage <- function(season,week,model,region,model_csv,target){
  if (week >= 40){
    analysis_time_season <- paste0(season,"/",as.numeric(season) + 1)
  } else{
    analysis_time_season <- paste0(as.numeric(season-1),"/",as.numeric(season))
  }
  formatted_region <- region_str_array_eval[match(region,region_str_data_set)]
  trajectory <- fully_observed_data[fully_observed_data$region == formatted_region & fully_observed_data$season == analysis_time_season ,]$weighted_ili
  truth <- max(trajectory)
  truth <- round(truth,1)
  truth_l <- max(0,truth-.1)
  truth_r <- min(13, truth+.1)
  
  
  
  
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
region_str_array_eval <- c("National",paste0("Region ",1:10))
region_str_data_set <- c("US National","HHS Region 1",  "HHS Region 2",  "HHS Region 3",  "HHS Region 4" ,"HHS Region 5",  "HHS Region 6",  "HHS Region 7",  "HHS Region 8",  "HHS Region 9","HHS Region 10")

region_str_true <- c("nat",paste0("hhs",1:10))


targets <- c("Season peak percentage","Season peak week","Season onset","1 wk ahead","2 wk ahead","3 wk ahead","4 wk ahead")
test_seasons <-c("2015","2016","2017")
test_models <- c("TRUE","NONE","M1","M2","M3","M4","M5","M6")
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
levels(result_df$model) <- c("Naive","Sampling","Naive (Week)","Naive (Region)","Hierarchical Regression","Non-linear","Unrevised","Revised")

mean(result_df[result_df$model=="Revised" &  result_df$target == "Season onset"   ,]$total_log_prob,na.rm = TRUE)
mean(result_df[result_df$model=="Unrevised"  &    result_df$target == "Season onset" ,]$total_log_prob,na.rm = TRUE)
mean(result_df[result_df$model=="Naive (Week)" &  result_df$target == "Season onset" ,]$total_log_prob,na.rm = TRUE)
mean(result_df[result_df$model=="Sampling" &   result_df$target == "Season onset" ,]$total_log_prob,na.rm = TRUE)

library(ggplot2)
library(dplyr)

season_week <- c()
for (i in 1:nrow(result_df)){
  season_week <- c(season_week,year_week_to_season_week(result_df[i,]$week,result_df[i,]$Season))
}

result_df$season_week <- season_week
ggplot(result_df[(result_df$model == "Unrevised" | result_df$model == "Revised"  )  & result_df$target =="Season peak percentage"  & result_df$Season == 2016 ,],aes(x=jitter(season_week),y=(total_log_prob),col=region)) + geom_point() + facet_grid(model ~.) + theme_bw() + ylab("Log Prob")  + xlab("Season week")

ggplot(result_df[result_df$Season == "2013"  & result_df$target == "1 wk ahead",],aes(x=week,y=total_log_prob,col=region)) + geom_point() + facet_grid(model ~.) 

result_df$Model <- result_df$model
result_df[result_df$model == "Sampling" & (result_df$target == "1 wk ahead" |result_df$target == "2 wk ahead" | result_df$target == "3 wk ahead" | result_df$target == "4 wk ahead"  ),]$total_log_prob <- NA

ggplot(result_df,aes(x=target,y=total_log_prob,col=model)) +  stat_summary(fun.y="mean", geom="point") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab("Average Log Probability") +  xlab("Target") + facet_grid(~Season) 
ggplot(result_df[(result_df$target!= "1 wk ahead" &result_df$target!= "2 wk ahead"& result_df$target!= "3 wk ahead"& result_df$target!= "4 wk ahead" ) & (result_df$model == "Unrevised" | result_df$model == "Revised" | result_df$model == "Sampling"),],aes(x=target,y=total_log_prob,col=model)) +  stat_summary(fun.y="mean", geom="point") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab("Average Log Probability") +  xlab("Target") + facet_wrap(~region,ncol=4) 

#levels(result_df$region) <- c("Region 1","Region 10",paste0("Region ",2:9),"National")
lag_df$region <- lag_df$Region


ggplot(result_df[(result_df$model == "M1" | result_df$model == "NONE" | result_df$model ==  "TRUE")& (result_df$target == "1 wk ahead" | result_df$target == "2 wk ahead" | result_df$target == "3 wk ahead" | result_df$target == "4 wk ahead"),],aes(x=target,y=total_log_prob,col=model)) +  stat_summary(fun.y="mean", geom="point") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab("Average Log Probability") +  xlab("Target") + facet_grid(~region) 


result_df_means <- data.frame(result_df %>% group_by(target,model,region,Season) %>% summarize(mt=mean(total_log_prob,na.rm = T)))


ggplot(result_df[(result_df$model == "Sampling" | result_df$model == "Revised" | result_df$model ==  "Unrevised")& (result_df$target == "Season onset" | result_df$target == "Season peak percentage"),],aes(x=target,y=total_log_prob,col=model)) +  stat_summary(fun.y="mean", geom="point") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab("Average Log Probability") +  xlab("Target")+ facet_grid(Season~region)  
ggplot(result_df_means[result_df_means$target == "1 wk ahead",],aes(x=model,y=mt,group=1)) + geom_line()+ facet_grid(Season~region) + theme_bw()

#print ("hello")
# for (region in unique(result_df$region)){
#   for (season in unique(result_df$Season)){
#     for (week in unique(result_df$week)){
#       for (model in unique(result_df$model)){
#         for (target in unique(result_df$target)){
# 
#           formatted_season <- paste0(season,"/",as.numeric(season)+1)
#           l_0 <- data[data$region == region_str_array_eval[match(region,region_str_data_set)] & data$season == formatted_season & data$week == week & data$lag == 0,]$weighted_ili
#           l_infty <- data[data$region == region_str_array_eval[match(region,region_str_data_set)] & data$season == formatted_season & data$week == week & data$lag == max(data[data$region == region_str_array_eval[match(region,region_str_data_set)] & data$season == formatted_season & data$week == week,]$lag),]$weighted_ili
#           result_df[result_df$Season ==season &result_df$target == target & result_df$week == week & result_df$region == region & result_df$model==model,]$ll <-l_0/l_infty
#         }
#       }
#     }
#   }
# }


ggplot(result_df[result_df$target == "1 wk ahead",],aes(x=ll,y=total_log_prob)) + geom_point() + geom_smooth() + theme_bw() + ylab("Log Prob") + xlab("Reporting Ratio") + ggtitle(label = "1 wk ahead for 2012/2013 & 2013/2014")




## Exploration plots
model_csv_1 <- read.csv(paste0("./inst/submissions/region-sarima_seasonal_difference_TRUE/EW10-",2016,"-ReichLab_sarima_seasonal_difference_TRUE-delay-","NONE",".csv"))
model_csv_2 <- read.csv(paste0("./inst/submissions/region-sarima_seasonal_difference_TRUE/EW04-",2016,"-ReichLab_sarima_seasonal_difference_TRUE-delay-","M2",".csv"))

model_csv <- rbind(model_csv_1,model_csv_2)
model_csv$m <-c( rep("1",nrow(model_csv_1)),rep("2",nrow(model_csv_2)))
library(forcats)
library(dplyr)
model_csv$Bin_start_incl <- as.factor(model_csv$Bin_start_incl)

model_csv_plot <- model_csv[model_csv$Location == "HHS Region 2" &model_csv$Target == "Season onset" & model_csv$Type == "Bin",] %>% 
  dplyr::mutate(Bin_start_incl = factor(Bin_start_incl, 
                                        levels = c(seq(40,52),seq(20)))) 

p1 <- ggplot(data_for_plot[data_for_plot$region == "Region 2" & data_for_plot$season == "2015/2016",],aes(x=season_week,y=weighted_ili,col=f)) + geom_point() + theme_bw() + geom_abline(intercept=get_onset_baseline("National",season="2015/2016"),slope=0) 
p2 <- ggplot(model_csv_plot,aes(x=Bin_start_incl,y=Value,color=m)) + geom_point() + theme_bw() 



##

for (week in unique(result_df$week)){
  print (result_df[result_df$week == week ,]$point)
}

result_df$point <- as.numeric(as.character(result_df$point))
result_df$truth <- as.numeric(as.character(result_df$truth))
ggplot(result_df[result_df$target == "1 wk ahead" & result_df$Season == 2016 ,],aes(x=season_week,y=point,col=model)) + geom_line() +
  geom_line(data=result_df[result_df$target == "1 wk ahead" & result_df$Season == 2016,],aes(x=season_week,y=truth),col='blue') + facet_grid(~region) + theme_bw()
}
