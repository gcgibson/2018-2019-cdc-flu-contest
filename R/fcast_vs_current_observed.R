do_explore <- function(){

mse <- matrix(NA,ncol=6)
for (season in c("2012/2013","2013/2014","2015/2016","2016/2017")){
  for (region in c("National",paste0("Region ",1:10))){
    for (week in c(seq(42,52),paste0("0",1:9),seq(10,20))){
      
      if (week == "01"){
        previous_week <- 52
      } else if(week <= 10){
        previous_week <- paste0("0",as.numeric(week) - 1)
      } else{
        previous_week <- as.numeric(week) -1
      }
      if (week >= 40){
        formatted_season <- substr(season,1,4)
      } else{
        formatted_season <- substr(season,6,9)
      }
      
      if (week == "01"){
        previous_formatted_season <- as.numeric(formatted_season) -1
      }else{
        previous_formatted_season <- formatted_season
      }
      print (c(formatted_season,previous_formatted_season,previous_week,week,region,season))
      fcast <- read.csv(paste0("./inst/submissions/region-sarima_seasonal_difference_TRUE/EW",previous_week,"-",previous_formatted_season,"-ReichLab_sarima_seasonal_difference_TRUE-delay-NONE.csv"))
      fcast_point <- fcast[fcast$Location == region_str_data_set[match(region,region_str_array_eval)] & fcast$Target == "1 wk ahead" & fcast$Type == "Point", ]$Value
      
      current_observed_data_with_lags <- data[ data$issue <= as.numeric(paste0(formatted_season,week)),]
      current_observed_data <- current_observed_data_with_lags %>% group_by(region,epiweek) %>%
        filter(lag == max(lag))
      current_observed_data <- data.frame(current_observed_data[order(current_observed_data$epiweek),])
      
      current_observed_data_point <- current_observed_data[current_observed_data$season == season & current_observed_data$region == region & current_observed_data$week == as.numeric(week),]$weighted_ili
      
      
      truth <- fully_observed_data[fully_observed_data$season == season & fully_observed_data$region == "National" & fully_observed_data$week == as.numeric(week),]$weighted_ili
      
      mse <- rbind(mse,c(season,region,week,(fcast_point-truth)^2,(current_observed_data_point-truth)^2,truth))

      
    }
  }
}

mse <- mse[2:nrow(mse),]
colnames(mse) <- c("season","region" ,"week","fcast","current_observed","truth")
mse <- data.frame(mse)
mse$fcast <- as.numeric(as.character(mse$fcast))
mse$current_observed <-as.numeric(as.character(mse$current_observed))
mse$week <- as.numeric(as.character(mse$week))
library(reshape2)
mse_long <- melt(mse,id.vars= c("season","region","week","truth"))

ggplot(mse_long, mapping = aes(x = week,y=value, color = variable)) + geom_point(position = position_jitter(),alpha=.01) + geom_smooth()
}