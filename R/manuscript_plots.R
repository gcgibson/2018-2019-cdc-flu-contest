
make_plots <- function(){
  lag_df <- read.csv("./data/lag_df")
  lag_df0 <- read.csv("./data/lag_df0")
  
  library(ggplot2)
  
  epw <- lag_df[lag_df$Region == "National" & lag_df$week == 201502,paste0("X",0:51)]
  
  plot_df <-data.frame(epw = as.numeric(t(epw)),x=0:51)
  ggplot(plot_df,aes(x=x,y=epw)) + geom_point() + theme_bw() + ylab("201502 wILI") + xlab("Lag")
  
  
  library(reshape2)
  
  epw <- lag_df[ lag_df$week == 201502,c(paste0("X",0:51),"Region")]
  epw <- epw[epw$Region == "National" |epw$Region =="Region 1"|epw$Region =="Region 5"|epw$Region =="Region 10",]
  long <- melt(epw, id.vars = c("Region"))
  
  plot_df <-data.frame(epw = as.numeric(t(epw)),x=rep(0:52,4))
  ggplot(long,aes(x=variable,y=value)) + geom_point(col='cornflowerblue') + theme_bw() + ylab("201502 wILI") +theme(axis.text.x = element_text(angle = 90, hjust = 1))+ xlab("Lag") +theme(axis.title.x=element_blank(),
                                                                                                                                                                                           axis.text.x=element_blank(),
                                                                                                                                                                                           axis.ticks.x=element_blank())+ facet_grid(. ~ Region)
  
  
  epw <- lag_df[ lag_df$week == 201440,c(paste0("X",0:51),"Region")]
  long <- melt(epw, id.vars = c("Region"))
  
  plot_df <-data.frame(epw = as.numeric(t(epw)),x=0:51)
  ggplot(long,aes(x=variable,y=value,col=Region)) + geom_point() + theme_bw() + ylab("201440 wILI") + xlab("Lag")+theme(axis.text.x = element_text(angle = 90, hjust = 1))

  lm_fit <- lm(X0~Incidence +season_week + Region, data=lag_df0)
  summary(lm_fit)
  lm_fit_hierarchical <- lmer(X0~Incidence +season_week + (1|Region), data=lag_df0)
  summary(lm_fit_hierarchical)
  nn_fit <- nnet(X0~Incidence +season_week + Region, subset_lag_df, size=12, maxit=500, linout=T, decay=0.01)
  
  
  library(ggplot2)
  
  ggplot(lag_df,aes(x=1:length(X0),y=X0,col=Region)) +geom_point() + theme(axis.title.x=element_blank(),
                                                                           axis.text.x=element_blank(),
                                                                           axis.ticks.x=element_blank())   +ylab("a_{r,*,*,0}")
  ggplot(lag_df,aes(x=1:length(X15),y=X15,col=Region)) +geom_point() + theme(axis.title.x=element_blank(),
                                                                             axis.text.x=element_blank(),
                                                                             axis.ticks.x=element_blank()) 
  
  ggplot(lag_df,aes(x=1:length(X0),y=X51,col=Region)) +geom_point() + theme(axis.title.x=element_blank(),
                                                                            axis.text.x=element_blank(),
                                                                            axis.ticks.x=element_blank()) 
  
  
  ggplot(lag_df,aes(x=season_week,y=X0,col=Region)) +geom_point() +theme_bw()
  
  
  ggplot(lag_df[lag_df$Region=="National",],aes(x=season_week,y=X0)) +geom_point(color='maroon',alpha=.5) +theme_bw()
  ggplot(lag_df0[lag_df0$Region=="National",],aes(x=Incidence,y=X0)) +geom_point(color='maroon',alpha=.5) +theme_bw()
  ggplot(lag_df0,aes(x=Region,y=X0)) +geom_point(color='maroon',alpha=.5) +theme_bw()
  
  
  predict_df <- data.frame(season_week = 1:53,Region=rep("National",53),week=1:53,X=1:53,Incidence=rep(3,53))
  lm.predict <- predict(lm_fit,newdata = predict_df)
  
  p <- ggplot(data.frame(season_week=1:53,X0=lm.predict), aes(x=season_week,y=X0)) 
  p <- p+ geom_point(data=lag_df0,aes(x=season_week,y=X0),alpha=.1,col="cornflowerblue") + theme_bw()+geom_line()
  p
  
  lm.predict <- predict(lm_fit_hierarchical,newdata = predict_df)
  
  p <- ggplot(data.frame(season_week=1:53,X0=lm.predict), aes(x=season_week,y=X0)) 
  p <- p+ geom_point(data=lag_df0,aes(x=season_week,y=X0),alpha=.1,col="cornflowerblue") + theme_bw()+geom_line()
  p
  
  lm.predict <- predict(nn_fit,newdata = predict_df)
  
  p <- ggplot(data.frame(season_week=1:53,X0=lm.predict), aes(x=season_week,y=X0)) 
  p <- p+ geom_point(data=lag_df0,aes(x=season_week,y=X0),alpha=.1,col="cornflowerblue") + theme_bw()+geom_line()
  p
  
  
  data <-readRDS("./data/flu_data_with_backfill_edit.rds")
  data_2015_subset <- data[data$epiweek <=201520 &data$epiweek >=201440 & data$region == "National", ]
  
  
  
  for (epiweek in c(seq(40,52),seq(20))){
    if (epiweek <= 20){
      year <- 2015
    }else{
      year <-2014
    }
    
    if (epiweek <= 9){
      formatted_week <- paste0("0",epiweek)
    } else{
      formatted_week <- epiweek
    }
    
    data_subset_2 <- data_2015_subset[data_2015_subset$issue <= paste0(year,formatted_week),]
    current_observed_data <- data_subset_2 %>% group_by(region,epiweek) %>%
      filter(lag == max(lag))
    current_observed_data <- data.frame(current_observed_data)
    current_observed_data <- current_observed_data[order(current_observed_data$epiweek),]
    #print (current_observed_data)
    p <-ggplot(current_observed_data,aes(x=1:length(weighted_ili),y=weighted_ili)) +  geom_point() + xlim(0,30) + ylim(0,8) + xlab("Season Week") +theme_bw()
    ggsave(paste0("2015_",epiweek,".png"),plot=p)
  }
  
  
}
  
  
  ###
  
