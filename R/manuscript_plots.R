
make_plots <- function(){
  library(reshape2)
  library(dplyr)
  library(ggplot2)
  lag_df <- read.csv("./data/lag_df")
  lag_df0 <- read.csv("./data/lag_df0")
  fully_observed <- readRDS("./data/fully_observed_data_formatted.rds")
  delayed_data <- readRDS("./data/flu_data_with_backfill_edit.rds")
  
  
  data <- readRDS("./data/flu_data_with_backfill_edit.rds")
  current_observed_data <- data %>% group_by(region,epiweek) %>%
    filter(lag == 0)
  current_observed_data <- data.frame(current_observed_data[order(current_observed_data$epiweek),])
  
  
  data_for_plot <- rbind(current_observed_data,fully_observed_data)
  
  data_for_plot$Data <- c(rep("Unrevised",nrow(current_observed_data)),rep("Revised",nrow(fully_observed_data)))
  season_onsets <- c()
  for (row in 1:nrow(data_for_plot)){
    season_onsets <- c(season_onsets,get_onset_baseline(data_for_plot[row,]$region,data_for_plot[row,]$season))
  }
  data_for_plot$season_onsets <- season_onsets
  ggplot(data_for_plot[ (data_for_plot$season == "2012/2013" |data_for_plot$season == "2014/2015" |data_for_plot$season == "2013/2014") ,],aes(x=season_week,y=weighted_ili,col=Data)) + geom_line() + geom_hline(aes(yintercept=season_onsets),linetype="dashed",alpha=.4) + facet_grid(season~region) + theme_bw() +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab("wILI") + xlab("Season Week")
  ggplot(data_for_plot[data_for_plot$region == "Region 8" & data_for_plot$season == "2015/2016",],aes(x=season_week,y=weighted_ili,col=Data)) + geom_line() + geom_hline(aes(yintercept=season_onsets),linetype="dashed",alpha=.4) + facet_grid(season~region) + theme_bw() +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab("wILI") + xlab("Season Week")
  
  
  
  ggplot(fully_observed[fully_observed_data$season >= "2010/2011" & fully_observed_data$season < "2017/2018",],aes(x=season_week,y=weighted_ili,col=region)) + geom_line() + facet_grid(~season) + theme_bw() +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab("wILI") + xlab("Season Week") +ggtitle("Revised Data")+ theme(plot.title = element_text(hjust = 0.5))

  long <- melt(lag_df[lag_df$week <=201540, ], id.vars = c("Region","season_week"))
    
  levels(long$variable) <- c("X.2","X.1","X"  ,"week",0:51)
  
  ggplot(long[long$variable != "X.1" & long$variable != "X.2" & long$variable != "X" & long$variable != "week", ],aes(x=variable,value,group=1))  +geom_point(aes(x=variable,value),alpha=.1,col='cornflowerblue') + stat_summary(fun.y="mean", geom="line")+ 
  ylab("a_{r,s,w,l}") + xlab("Lag")+ theme_bw() + facet_wrap(~Region)  + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  
  

  ggplot(long[(long$variable == "X0" |long$variable == "X5" |long$variable == "X15" |long$variable == "X52"|long$variable == "X10"|long$variable == "X30") & (long$Region == "National" | long$Region == "Region 1" | long$Region == "Region 5" | long$Region == "Region 10")    ,],aes(x=season_week,y=value))  + geom_point() + facet_grid(Region~variable) + xlab("a_{r,s,w,0}") + theme_bw()
  
  ggplot(long[long$variable == "X0"   ,],aes(x=value))  + geom_histogram() + facet_grid(Region~variable) + xlab("a_{r,s,w,0}") + theme_bw()
  
  
  
  data_sub <- data[data$season == "2015/2016"  & data$epiweek >="201540",]
 
   data_frame_for_lag_plot <- matrix(NA,ncol=21)
   for (i in seq(41,52)){
     data_frame_for_lag_plot <- rbind(data_frame_for_lag_plot,as.matrix(data_sub[ data_sub$issue <= as.numeric(paste0("2015",i)),] %>% group_by(region,epiweek) %>%
                filter(lag == max(lag))))
   }
   for (i in c(paste0("0",1:9),seq(10,20))){
     data_frame_for_lag_plot <- rbind(data_frame_for_lag_plot,as.matrix(data_sub[ data_sub$issue <= as.numeric(paste0("2016",i)),] %>% group_by(region,epiweek) %>%
                                                                          filter(lag == max(lag))))
   }
   data_frame_for_lag_plot <- data.frame(data_frame_for_lag_plot[2:nrow(data_frame_for_lag_plot),])
   data_frame_for_lag_plot$weighted_ili <- as.numeric(as.character(data_frame_for_lag_plot$weighted_ili))
    #annotate("text", x = c(25), y=1.8, label = c("Y_{r,s,15,7}"))
   data_frame_for_lag_plot$lag <- as.factor(as.numeric(as.character(data_frame_for_lag_plot$lag)))
   data_frame_for_lag_plot$season_week <- as.numeric(as.character(data_frame_for_lag_plot$season_week))
   ggplot(data_frame_for_lag_plot,aes(x=season_week,y=weighted_ili,col=lag)) + geom_point(size=.5) + theme_bw() +facet_wrap(~region) +
     geom_line(data=fully_observed_data[fully_observed_data$season  == "2015/2016",],aes(x=season_week,y=weighted_ili,col=as.factor("Revised")),color="black") + ylab("wILI") + xlab("Season week")+ggtitle("Revisions to 2015/2016")+ theme(plot.title = element_text(hjust = 0.5))
  
  
    
  
}
  
  
  ###
  
