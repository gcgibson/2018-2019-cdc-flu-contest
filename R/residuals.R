library(EnvStats)
library(ggplot2)
library(nnet)
fully_observed_data <- readRDS("./data/fully_observed_data_formatted.rds")
#model_params <- read.csv("model_params.csv")

region_str_array_eval <- c("National",paste0("Region ",1:10))
region_str_data_set <- c("US National","HHS Region 1",  "HHS Region 2",  "HHS Region 3",  "HHS Region 4" ,"HHS Region 5",  "HHS Region 6",  "HHS Region 7",  "HHS Region 8",  "HHS Region 9","HHS Region 10")
region_str_true <- c("nat",paste0("hhs",1:10))


targets <- c("1 wk ahead","2 wk ahead","3 wk ahead","4 wk ahead","Season onset","Season peak week", "Season peak percentage")
#generate residuals for training seasons
for (target in targets){
  data_for_plot <- matrix(NA,ncol=3)
  
  for (top_level_season in c( "2011", "2012", "2013")){
    step_ahead <- as.numeric(substr(targets[1],1,1))
    for (test_region in region_str_data_set ){
        for (test_season in top_level_season){
        if (test_season == "2017"){
          end_week <- 12
        }else{
          end_week <- 20
        }
        for (test_week in c(seq(41,52),seq(end_week))){
          if (test_week < 40){
            test_season_formatted <- as.numeric(test_season) + 1
            if (test_week <= 9){
              test_week_formatted <- paste0("0",test_week)
            }else{
              test_week_formatted <- test_week
            }
          }else{
            test_season_formatted <- test_season
            test_week_formatted <- test_week
          }
          
          

          non_delay_adjusted_forecasts <- read.csv(paste0("./inst/submissions/region-sarima_seasonal_difference_TRUE/EW",test_week_formatted,"-",test_season_formatted,"-ReichLab_sarima_seasonal_difference_TRUE-delay-NONE.csv"))

          truth_time <- NULL
          if (as.numeric(test_week_formatted) + step_ahead > 52){
            tmp_week <- as.numeric(test_week_formatted) +step_ahead-52
            if (tmp_week <= 9){
              tmp_week_formatted <- paste0("0",tmp_week)
            }else{
              tmp_week_formatted <- tmp_week
            }
            truth_time <- paste0(as.numeric(test_season_formatted)+1,tmp_week_formatted)
          }else{
            tmp_week <- as.numeric(test_week_formatted) +step_ahead
            if (tmp_week <= 9){
              tmp_week_formatted <- paste0("0",tmp_week)
            }else{
              tmp_week_formatted <- tmp_week
            }
            truth_time <- paste0(test_season_formatted,tmp_week_formatted)
          }
          point_forecast_non_delayed <- non_delay_adjusted_forecasts[non_delay_adjusted_forecasts$Target == "1 wk ahead" & non_delay_adjusted_forecasts$Type=="Point" & non_delay_adjusted_forecasts$Location== test_region,]$Value
          point_forecast_true <-fully_observed_data[fully_observed_data$region == region_str_array_eval[match(test_region,region_str_data_set)] &fully_observed_data$epiweek == truth_time,]$wili

          data_for_plot <- rbind(data_for_plot,c(test_region,truth_time,point_forecast_true - point_forecast_non_delayed))
          }
        }
      }
  }
  data_for_plot <- data.frame(data_for_plot[2:nrow(data_for_plot),])
  colnames(data_for_plot) <- c("Region","epiweek","residual")
  data_for_plot$week <- unlist(lapply(data_for_plot$epiweek,function(x){as.numeric(as.character(substr(x,5,7)))}))
  data_for_plot$residual <- as.numeric(as.character(data_for_plot$residual))
  
  loess_fit <- nnet(residual~ week + Region, data_for_plot,size=20)
  predict_df <- data.frame(Region="US National",week = c(seq(40,52),seq(1,20))) 
  predict(loess_fit,predict_df)
  saveRDS(loess_fit, file = paste0(gsub(" ", "", target) ,"_week_ahead_residual_fit.rda"))
  
  }




hist(as.numeric(as.character(data_for_plot$residual)))


plot(data_for_plot$week,data_for_plot$residual)
library(ggplot2)

ggplot(data_for_plot,aes(week,residual)) +geom_point() +facet_wrap(~Region)



