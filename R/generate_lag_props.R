download_backfill_data <- function(){
  library(plyr) # for rbind.fill
  library(dplyr)
  source("https://raw.githubusercontent.com/cmu-delphi/delphi-epidata/master/src/client/delphi_epidata.R")
  
  # Fetch data
  all_obs <- lapply(c("nat", paste0("hhs", 1:10)),
                    function(region_val) {
                      lapply(0:51,
                             function(lag_val) {
                               obs_one_lag <- Epidata$fluview(
                                 regions = list(region_val),
                                 epiweeks = list(Epidata$range(199740, 201815)),
                                 lag = list(lag_val))
                               
                               lapply(obs_one_lag$epidata,
                                      function(x) {
                                        x[sapply(x, function(comp) is.null(comp))] <- NA
                                        return(as.data.frame(x))
                                      }) %>%
                                 rbind.fill()
                             }) %>%
                        rbind.fill()
                    }) %>%
    rbind.fill()
  
  saveRDS(all_obs,
          file = "flu_data_with_backfill.rds")
  
}


create_lag_df <- function(){
  data <- readRDS("./data/flu_data_with_backfill_edit.rds")
  
  lag_df <- matrix(NA,ncol=4)
  
  for (region in unique(data$region)){
    for (week in unique(data[data$region == region,]$epiweek)){
      tmp_data <- data[data$region == region & data$epiweek == week,]
      tmp_row <- c()
      for (lag in c(0)){
        current_observed_data <- tmp_data[tmp_data$lag == lag,]$wili
        finally_observed_data <- tmp_data[tmp_data$lag == max(tmp_data$lag),]$wili
        prop <- current_observed_data/finally_observed_data
        tmp_row <- c(tmp_row,prop)
      }
      # while (length(tmp_row) < 52){
      #   tmp_row <- c(tmp_row, NA)
      # }
      if (length(prop) ){
        lag_df <- rbind(lag_df,c(region,week,tmp_row,current_observed_data))
      }
    }
  }
  
  lag_df <- as.data.frame(lag_df)
  lag_df <- lag_df[2:nrow(lag_df),]
  colnames(lag_df) <- c("Region","week",0,"Incidence")
  lag_df$season_week <- unlist(lapply(lag_df$week,function(x) {return (substr(x,5,7))}))
  
  write.csv(lag_df,"./data/lag_df0")
}


lag_df <- as.data.frame(read.csv("./data/lag_df0"))
data <- as.data.frame(readRDS("./data/flu_data_with_backfill_edit.rds"))

subset_lag_df <- lag_df[lag_df$week < 201540,]
validation_lag_df <- lag_df[lag_df$week >= 201540,]
library(nlme)
library(MASS)
library(nnet)


lme1 <- lme(fixed=X0 ~ Region +season_week, random = ~ +1|Region,data=subset_lag_df,na.action=na.exclude)
lmfit <- lme(X0 ~ Incidence +season_week,random = ~ +1|Region,data=subset_lag_df,na.action=na.exclude)
lmfit2 <- lm(X0 ~ 1,data=subset_lag_df)
lmfit3 <- lm(X0 ~ Region + season_week+ Incidence,data=subset_lag_df)
#loess_fit <- loess(X0~ Incidence +season_week , data=subset_lag_df)
loess_fit <- nnet(X0~Incidence +season_week + Region, subset_lag_df, size=12, maxit=500, linout=T, decay=0.01)
subset_lag_df$binary <- ifelse(subset_lag_df$X0 < 1,0,1)

glm_fit <- glm(binary ~ Incidence +season_week + Region,data=subset_lag_df)
mse_m1 <- c()
mse_m2 <- c()
mse_m3 <- c()
mse_m4 <- c()
mse_m5 <- c()
mse_m6 <- c()
mse_m7 <- c()


for (region in c("National",paste0("Region ",1:10))){
  for (week in c(seq(41,52),seq(20))){
    if (week <=9){
      data_week <- paste0("0",week)
    }else{
      data_week <- week
    }
    if (week < 40){
      test_season <- 2016
    }else{
      test_season <- 2015
    }
    test_dat <- data.frame(Region=region,season_week = week,Incidence=data[data$region==region & data$epiweek==paste0(test_season,data_week) & data$lag==0,]$weighted_ili)
    
    if (data_week == 1){
      prop1 <- data[data$region==region & data$epiweek==paste0(test_season-1,52) & data$lag==0,]$weighted_ili
      
      prop2 <- data[data$region==region & data$epiweek==paste0(test_season-1,52) & data$lag==1,]$weighted_ili
    } else if(data_week <=10){
      prop1 <- data[data$region==region & data$epiweek==paste0(test_season,paste0("0",as.numeric(data_week)-1)) & data$lag==0,]$weighted_ili
      
      prop2 <- data[data$region==region & data$epiweek==paste0(test_season,paste0("0",as.numeric(data_week)-1)) & data$lag==1,]$weighted_ili
    }else{
        prop2 <-data[data$region==region & data$epiweek==paste0(test_season,data_week-1) & data$lag==1,]$weighted_ili
        prop1 <-data[data$region==region & data$epiweek==paste0(test_season,data_week-1) & data$lag==0,]$weighted_ili
        
    }
    
    m1_pred <- predict(lmfit,newdata = test_dat)
    m2_pred <- predict(lme1,newdata = test_dat)
    m3_pred <- predict(lmfit2,newdata=test_dat)
    m4_pred <- predict(lmfit3,newdata=test_dat)
    m5_pred <-predict(loess_fit,newdata=test_dat)
    m6_pred <- ifelse(predict(glm_fit,newdata=test_dat)<.5,.99,1.01)
    
    
    truth <- lag_df[lag_df$Region==region & lag_df$week==paste0(test_season,week),"X0" ]
    mse_m1 <- c(mse_m1,(m1_pred-truth)^2)
    mse_m2 <- c(mse_m2,(m2_pred-truth)^2)
    mse_m3 <- c(mse_m3,(m3_pred-truth)^2)
    mse_m4 <- c(mse_m4,(m4_pred-truth)^2)
    mse_m5 <- c(mse_m5,(m5_pred-truth)^2)
    mse_m6 <- c(mse_m6,(m6_pred-truth)^2)
    mse_m7 <- c(mse_m7,((prop1/prop2+m5_pred)/2-truth)^2)
  }
}

print (mean(mse_m1))
print (mean(mse_m3))
print (mean(mse_m2))
print (mean(mse_m4))
print (mean(mse_m5))
print (mean(mse_m6))
print (mean(mse_m7))
library(ggplot2)
ggplot(subset_lag_df[subset_lag_df$Region=="Region 1",],aes(x=season_week,y=X0))  +geom_point()

plot(subset_lag_df[subset_lag_df$Region=="Region 1",]$X0)
