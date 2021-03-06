library(EnvStats)
library(ggplot2)
do_eval_retro <- function(){
  get_model_prob <- function(season,week,target,model,region) {
    model_csv <- read.csv(paste0("./inst/submissions/region-sarima_seasonal_difference_TRUE/EW",week,"-",season,"-ReichLab_sarima_seasonal_difference_TRUE-delay-",model,".csv"))
    if(target == "1 wk ahead" | target == "2 wk ahead" | target == "3 wk ahead" | target == "4 wk ahead"){
      prob_and_truth <- get_model_prob_k_week_ahead(season,week,model,region,model_csv,target)
    } else if (target == "Season onset"){
      prob_and_truth <- get_model_prob_season_onset(season,week,model,region,model_csv,target)
    } else if (target == "Season peak week"){
      prob_and_truth <- get_model_prob_season_peak_week(season,week,model,region,model_csv,target)
    } else if (target == "Season peak percentage"){
      prob_and_truth <- get_model_prob_season_peak_percentage(season,week,model,region,model_csv,target)
    }
    return (prob_and_truth)
  }
  
  
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

    if (onset_week != "none"){
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
    }else{
      truth <- onset_week
      truth_l <-onset_week
      truth_r <- onset_week
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
  test_seasons <-c("2011","2012","2013", "2014")
  test_models <- c("NONE","TRUE")#c("FSMOOTHED","TRUE","NONE","M1","M2","M3","M4","M5","M6")
  test_regions <- region_str_data_set
  
  train_result_df <- matrix(NA,ncol= 10)
  
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
           # print (c(top_level_season_formatted,week,target,model,region))
            model_prob_and_truth <- get_model_prob(top_level_season_formatted,week,target,model,region) 
            
            train_result_df <- rbind(train_result_df,c(top_level_season,week,target,model,region,model_prob_and_truth))
          }
        }   
      }
    }
  }
  
  train_result_df <- train_result_df[2:nrow(train_result_df),]
  
  colnames(train_result_df) <- c("Season","week","target","model","region","p_c","p_l","p_r","truth","point")
  
  train_result_df <- as.data.frame(train_result_df)
  
  train_result_df$p_c <- as.numeric(as.character(train_result_df$p_c))
  train_result_df$p_r <- as.numeric(as.character(train_result_df$p_r))
  train_result_df$p_l <- as.numeric(as.character(train_result_df$p_l))
  train_result_df$week <- as.numeric(as.character(train_result_df$week))
  train_result_df$total_prob <- train_result_df$p_c +train_result_df$p_r +train_result_df$p_l
  
  train_result_df$total_log_prob <- pmax(log(train_result_df$total_prob),-10)
  train_result_df$ll <- rep(NA, nrow(train_result_df))
  
  mean(train_result_df[train_result_df$model=="TRUE" &  train_result_df$target == "Season peak week"   ,]$total_log_prob,na.rm = TRUE)
  mean(train_result_df[train_result_df$model=="NONE"  &    train_result_df$target == "Season peak week" ,]$total_log_prob,na.rm = TRUE)
  mean(train_result_df[train_result_df$model=="M1" &  train_result_df$target == "Season peak percentage" ,]$total_log_prob,na.rm = TRUE)
  mean(train_result_df[train_result_df$model=="M2" &   train_result_df$target == "Season peak week" ,]$total_log_prob,na.rm = TRUE)
  
  train_result_df <- readRDS("train_result_df")
  library(ggplot2)
  library(dplyr)
  
  season_week <- c()
  for (i in 1:nrow(train_result_df)){
    season_week <- c(season_week,year_week_to_season_week(train_result_df[i,]$week,train_result_df[i,]$Season))
  }
  
  train_result_df$season_week <- season_week
  ggplot(train_result_df[(train_result_df$model == "Unrevised" | train_result_df$model == "Revised"  )  & train_result_df$target =="Season peak percentage"  & train_result_df$Season == 2016 ,],aes(x=jitter(season_week),y=(total_log_prob),col=region)) + geom_point() + facet_grid(model ~.) + theme_bw() + ylab("Log Prob")  + xlab("Season week")
  
  ggplot(train_result_df[train_result_df$Season == "2013"  & train_result_df$target == "1 wk ahead",],aes(x=week,y=total_log_prob,col=region)) + geom_point() + facet_grid(model ~.) 
  
  
  #levels(train_result_df$model) <-c("FSMOOTHED","Mean","Sampling","M3","Mean/Season Week","Hierarchical","Nonlinear","None","True") 
  
  train_result_df_means <- data.frame(train_result_df %>% group_by(target,model,region,Season) %>% summarize(mt=mean(total_log_prob,na.rm = T)))
  long <- melt(train_result_df_means, id.vars = c("region","Season","mt"))
  
  levels(train_result_df_means$model) <- c("Unrevised","Revised")
  train_result_df$diff <- train_result_df[train_result_df$model == "NONE",]$total_log_prob-train_result_df[train_result_df$model == "TRUE",]$total_log_prob
  ggplot(train_result_df,aes(x=region,y=diff)) + geom_violin() + theme_bw()  +facet_wrap(~target) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab("Difference in Log Score") +ggtitle("Difference Log-scores for Revised vs Unrevised")+ theme(plot.title = element_text(hjust = 0.5)) 
    #
  #levels(train_result_df$region) <- c("Region 1","Region 10",paste0("Region ",2:9),"National")
  lag_df$region <- lag_df$Region
  
  
  ggplot(train_result_df[(train_result_df$model == "M1" | train_result_df$model == "NONE" | train_result_df$model ==  "TRUE")& (train_result_df$target == "1 wk ahead" | train_result_df$target == "2 wk ahead" | train_result_df$target == "3 wk ahead" | train_result_df$target == "4 wk ahead"),],aes(x=target,y=total_log_prob,col=model)) +  stat_summary(fun.y="mean", geom="point") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab("Average Log Probability") +  xlab("Target") + facet_grid(~region) 
  
  
  train_result_df_means <- data.frame(train_result_df %>% group_by(target,model,region,Season) %>% summarize(mt=mean(total_log_prob,na.rm = T)))
  
  
  ggplot(train_result_df[(train_result_df$model == "M2" | train_result_df$model == "NONE" | train_result_df$model ==  "TRUE")& (train_result_df$target == "Season onset" | train_result_df$target == "Season peak percentage"),],aes(x=target,y=total_log_prob,col=model)) +  stat_summary(fun.y="mean", geom="point") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab("Average Log Probability") +  xlab("Target")  +
    geom_point(data=lag_df,aes(x=as.factor("Variance of X0"),y=min(X0,na.rm = T)),col='black') + facet_grid(~region)
  ggplot(train_result_df_means[train_result_df_means$target == "1 wk ahead",],aes(x=model,y=mt,group=1)) + geom_line()+ facet_grid(Season~region) + theme_bw()
  
  #print ("hello")
   # for (region in unique(train_result_df$region)){
   #   for (season in unique(train_result_df$Season)){
   #     for (week in unique(train_result_df$week)){
   #       for (model in unique(train_result_df$model)){
   #         for (target in unique(train_result_df$target)){
   # 
   #           formatted_season <- paste0(season,"/",as.numeric(season)+1)
   #           l_0 <- data[data$region == region_str_array_eval[match(region,region_str_data_set)] & data$season == formatted_season & data$week == week & data$lag == 0,]$weighted_ili
   #           l_infty <- data[data$region == region_str_array_eval[match(region,region_str_data_set)] & data$season == formatted_season & data$week == week & data$lag == max(data[data$region == region_str_array_eval[match(region,region_str_data_set)] & data$season == formatted_season & data$week == week,]$lag),]$weighted_ili
   #           train_result_df[train_result_df$Season ==season &train_result_df$target == target & train_result_df$week == week & train_result_df$region == region & train_result_df$model==model,]$ll <-l_0/l_infty
   #         }
   #       }
   #     }
   #   }
   # }
  
  
  ggplot(train_result_df[train_result_df$target == "1 wk ahead",],aes(x=ll,y=total_log_prob)) + geom_point() + geom_smooth() + theme_bw() + ylab("Log Prob") + xlab("Reporting Ratio") + ggtitle(label = "1 wk ahead for 2012/2013 & 2013/2014")
  
  
  
  
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
  
  for (week in unique(train_result_df$week)){
    print (train_result_df[train_result_df$week == week ,]$point)
  }
  
  train_result_df$point <- as.numeric(as.character(train_result_df$point))
  train_result_df$truth <- as.numeric(as.character(train_result_df$truth))
  ggplot(train_result_df[train_result_df$target == "1 wk ahead" & train_result_df$Season == 2016 ,],aes(x=season_week,y=point,col=model)) + geom_line() +
    geom_line(data=train_result_df[train_result_df$target == "1 wk ahead" & train_result_df$Season == 2016,],aes(x=season_week,y=truth),col='blue') + facet_grid(~region) + theme_bw()
}
