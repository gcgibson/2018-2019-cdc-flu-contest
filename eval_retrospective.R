library(EnvStats)
library(ggplot2)
fully_observed_data <- readRDS("./data/fully_observed_data_formatted.rds")
#model_params <- read.csv("model_params.csv")

lag_df0 <- read.csv("lag_df0")

region_str_array_eval <- c("National",paste0("Region ",1:10))
region_str_data_set <- c("US National","HHS Region 1",  "HHS Region 2",  "HHS Region 3",  "HHS Region 4" ,"HHS Region 5",  "HHS Region 6",  "HHS Region 7",  "HHS Region 8",  "HHS Region 9","HHS Region 10")
region_str_true <- c("nat",paste0("hhs",1:10))
#test_region <-toString(model_params$region)
#model_var <- model_params$model_variance

targets <- c("1 wk ahead","2 wk ahead","3 wk ahead","4 wk ahead","Season onset","Season peak week", "Season peak percentage")
data_for_plot <- matrix(NA,ncol=10)

for (target in targets){
  residual_fit <- readRDS(paste0(gsub(" ", "", target),"_week_ahead_residual_fit.rda"))
  for (top_level_season in c("2017")){
    
    step_ahead <- as.numeric(substr(targets[1],1,1))
    score <- "MULTIBIN"
    gm_mean = function(x, na.rm=TRUE){
      exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
    }
    
    
    
    delay_adjusted_total_prob_m1 <- c()
    delay_adjusted_total_prob_m2 <- c()
    delay_adjusted_total_prob_m3 <- c()
    delay_adjusted_total_prob_m4 <- c()
    delay_adjusted_total_prob_m5 <- c()
    delay_adjusted_total_prob_m6 <- c()
    delay_adjusted_total_prob_m7 <- c()
    
    
    non_delay_adjusted_total_prob <- c()
    
    for (test_region in region_str_data_set ){
    
      true_total_prob <- c()
      for (test_season in top_level_season){
        if (test_season == "2017"){
          end_week <- 20
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
          delay_adjusted_forecasts_m1 <- read.csv(paste0("./inst/submissions/region-sarima_seasonal_difference_TRUE/EW",test_week_formatted,"-",test_season_formatted,"-ReichLab_sarima_seasonal_difference_TRUE-delay-M1.csv"))
          
          if (test_season_formatted == "2017"){
            delay_adjusted_forecasts_m2 <-delay_adjusted_forecasts_m1
          } else{
            delay_adjusted_forecasts_m2 <- read.csv(paste0("./inst/submissions/region-sarima_seasonal_difference_TRUE/EW",test_week_formatted,"-",test_season_formatted,"-ReichLab_sarima_seasonal_difference_TRUE-delay-M2.csv"))
            
          }
          delay_adjusted_forecasts_m3 <- read.csv(paste0("./inst/submissions/region-sarima_seasonal_difference_TRUE/EW",test_week_formatted,"-",test_season_formatted,"-ReichLab_sarima_seasonal_difference_TRUE-delay-M3.csv"))
          delay_adjusted_forecasts_m4 <- read.csv(paste0("./inst/submissions/region-sarima_seasonal_difference_TRUE/EW",test_week_formatted,"-",test_season_formatted,"-ReichLab_sarima_seasonal_difference_TRUE-delay-M4.csv"))
          delay_adjusted_forecasts_m5 <- read.csv(paste0("./inst/submissions/region-sarima_seasonal_difference_TRUE/EW",test_week_formatted,"-",test_season_formatted,"-ReichLab_sarima_seasonal_difference_TRUE-delay-M5.csv"))
          delay_adjusted_forecasts_m6 <- read.csv(paste0("./inst/submissions/region-sarima_seasonal_difference_TRUE/EW",test_week_formatted,"-",test_season_formatted,"-ReichLab_sarima_seasonal_difference_TRUE-delay-M6.csv"))
          
          non_delay_adjusted_forecasts <- read.csv(paste0("./inst/submissions/region-sarima_seasonal_difference_TRUE/EW",test_week_formatted,"-",test_season_formatted,"-ReichLab_sarima_seasonal_difference_TRUE-delay-NONE.csv"))
          true_forecasts <- read.csv(paste0("./inst/submissions/region-sarima_seasonal_difference_TRUE/EW",test_week_formatted,"-",test_season_formatted,"-ReichLab_sarima_seasonal_difference_TRUE-delay-TRUE.csv"))
          truth_time <- NULL
          
          
          if (target == "1 wk ahead" || target == "2 wk ahead"|| target == "3 wk ahead" || target == "4 wk ahead"  ){
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
            truth <- get_inc_bin(fully_observed_data[fully_observed_data$region == region_str_array_eval[match(test_region,region_str_data_set)] &fully_observed_data$epiweek == truth_time,]$wili)
            truth_l <- max(0,as.numeric(truth)-.1)
            truth_ll <- max(0,as.numeric(truth)-.2)
            truth_r <- min(13,as.numeric(truth)+.1)
            truth_rr <- min(13,as.numeric(truth)+.2)
          } else if(target == "Season onset") {
            baseline <- get_onset_baseline(region = region_str_array_eval[match(test_region,region_str_data_set)], season = "2016/2017")
            vec_tru <- fully_observed_data[fully_observed_data$region == region_str_array_eval[match(test_region,region_str_data_set)] &fully_observed_data$season == "2016/2017",]$wili >baseline
            unformatted_week <- min(which(vec_tru==TRUE))
            if (unformatted_week <=21){
              unformatted_week <- unformatted_week+30
            }else{
              unformatted_week <- unformatted_week-21
            }
            truth <- unformatted_week
            truth_l <- unformatted_week -1
            truth_ll <- unformatted_week -2
            truth_r <- unformatted_week + 1
            truth_rr <- unformatted_week + 2
            
          } else if (target == "Season peak week"){
            baseline <- get_onset_baseline(region = region_str_array_eval[match(test_region,region_str_data_set)], season = "2016/2017")
            vec_tru <- fully_observed_data[fully_observed_data$region == region_str_array_eval[match(test_region,region_str_data_set)] &fully_observed_data$season == "2016/2017",]$wili
            unformatted_week <- which.max(vec_tru)
            if (unformatted_week <=21){
              unformatted_week <- unformatted_week+30
            }else{
              unformatted_week <- unformatted_week-21
            }
            truth <- unformatted_week
            truth_l <- unformatted_week -1
            truth_ll <- unformatted_week -2
            truth_r <- unformatted_week + 1
            truth_rr <- unformatted_week + 2
            
          } else if (target == "Season peak percentage"){
            vec_tru <- fully_observed_data[fully_observed_data$region == region_str_array_eval[match(test_region,region_str_data_set)] &fully_observed_data$season == "2016/2017",]$wili
            unformatted_max <- vec_tru[which.max(vec_tru)]
            
            truth <- round(unformatted_max,1)
            truth_l <- unformatted_max -.1
            truth_ll <- unformatted_max -.2
            truth_r <- unformatted_max + .1
            truth_rr <- unformatted_max + .2
          }
          
          
          
          
          pred_df <- data.frame(Region =test_region,week=as.numeric(truth))
          adjustment_factor <- mean(lag_df[lag_df$Region == region_str_array_eval[match(test_region,region_str_data_set)] &lag_df$season_week ==tmp_week_formatted,]$X0,na.rm = T) #predict(residual_fit,pred_df)
          if (is.nan(adjustment_factor )){
            adjustment_factor <- 1
          } 
         
          
          prob_delay_adjusted_center_m1 <- delay_adjusted_forecasts_m1[delay_adjusted_forecasts_m1$Target==target & delay_adjusted_forecasts_m1$Location == test_region& delay_adjusted_forecasts_m1$Bin_start_incl == truth,]$Value[-1]
          prob_delay_adjusted_center_m2 <- delay_adjusted_forecasts_m2[delay_adjusted_forecasts_m2$Target==target & delay_adjusted_forecasts_m2$Location == test_region & delay_adjusted_forecasts_m2$Bin_start_incl == truth,]$Value[-1]
          prob_delay_adjusted_center_m3 <- delay_adjusted_forecasts_m3[delay_adjusted_forecasts_m3$Target==target & delay_adjusted_forecasts_m3$Location == test_region & delay_adjusted_forecasts_m3$Bin_start_incl == truth,]$Value[-1]
          prob_delay_adjusted_center_m4 <- delay_adjusted_forecasts_m4[delay_adjusted_forecasts_m4$Target==target & delay_adjusted_forecasts_m4$Location == test_region & delay_adjusted_forecasts_m4$Bin_start_incl == truth,]$Value[-1]
          prob_delay_adjusted_center_m5 <- delay_adjusted_forecasts_m5[delay_adjusted_forecasts_m5$Target==target & delay_adjusted_forecasts_m5$Location == test_region & delay_adjusted_forecasts_m5$Bin_start_incl == truth,]$Value[-1]
          prob_delay_adjusted_center_m6 <- delay_adjusted_forecasts_m6[delay_adjusted_forecasts_m6$Target==target & delay_adjusted_forecasts_m6$Location == test_region & delay_adjusted_forecasts_m6$Bin_start_incl == truth,]$Value[-1]
          
          
          prob_delay_adjusted_center_m7 <- non_delay_adjusted_forecasts[non_delay_adjusted_forecasts$Target==target & non_delay_adjusted_forecasts$Location == test_region & non_delay_adjusted_forecasts$Bin_start_incl == min(max(round(as.numeric(truth)*adjustment_factor,1),0),13),]$Value[-1]
          
          
          prob_non_delay_adjusted_center <- non_delay_adjusted_forecasts[non_delay_adjusted_forecasts$Target==target & non_delay_adjusted_forecasts$Location ==test_region& non_delay_adjusted_forecasts$Bin_start_incl == truth,]$Value[-1]
          prob_true_center <- true_forecasts[true_forecasts$Target==target & true_forecasts$Location ==test_region& true_forecasts$Bin_start_incl == truth,]$Value[-1]
          
          prob_delay_adjusted_l_m1 <- delay_adjusted_forecasts_m1[delay_adjusted_forecasts_m1$Target==target & delay_adjusted_forecasts_m1$Location == test_region & delay_adjusted_forecasts_m1$Bin_start_incl == truth_l,]$Value[-1]
          prob_delay_adjusted_l_m2 <- delay_adjusted_forecasts_m2[delay_adjusted_forecasts_m2$Target==target & delay_adjusted_forecasts_m2$Location == test_region& delay_adjusted_forecasts_m2$Bin_start_incl == truth_l,]$Value[-1]
          prob_delay_adjusted_l_m3 <- delay_adjusted_forecasts_m3[delay_adjusted_forecasts_m3$Target==target & delay_adjusted_forecasts_m3$Location == test_region& delay_adjusted_forecasts_m3$Bin_start_incl == truth_l,]$Value[-1]
          prob_delay_adjusted_l_m4 <- delay_adjusted_forecasts_m4[delay_adjusted_forecasts_m4$Target==target & delay_adjusted_forecasts_m4$Location == test_region& delay_adjusted_forecasts_m4$Bin_start_incl == truth_l,]$Value[-1]
          prob_delay_adjusted_l_m5 <- delay_adjusted_forecasts_m5[delay_adjusted_forecasts_m5$Target==target & delay_adjusted_forecasts_m5$Location == test_region& delay_adjusted_forecasts_m5$Bin_start_incl == truth_l,]$Value[-1]
          prob_delay_adjusted_l_m6 <- delay_adjusted_forecasts_m6[delay_adjusted_forecasts_m6$Target==target & delay_adjusted_forecasts_m6$Location == test_region& delay_adjusted_forecasts_m6$Bin_start_incl == truth_l,]$Value[-1]
          prob_delay_adjusted_l_m7 <- non_delay_adjusted_forecasts[non_delay_adjusted_forecasts$Target==target & non_delay_adjusted_forecasts$Location == test_region & non_delay_adjusted_forecasts$Bin_start_incl == min(max(round(as.numeric(truth_l)*adjustment_factor,1),0),13),]$Value[-1]
          
          
          prob_non_delay_adjusted_l <- non_delay_adjusted_forecasts[non_delay_adjusted_forecasts$Target==target & non_delay_adjusted_forecasts$Location ==test_region& non_delay_adjusted_forecasts$Bin_start_incl == truth_l,]$Value[-1]
          prob_true_l <- true_forecasts[true_forecasts$Target==target & true_forecasts$Location ==test_region& true_forecasts$Bin_start_incl == truth_l,]$Value[-1]
          
          prob_delay_adjusted_ll_m1 <- delay_adjusted_forecasts_m1[delay_adjusted_forecasts_m1$Target==target & delay_adjusted_forecasts_m1$Location == test_region& delay_adjusted_forecasts_m1$Bin_start_incl == truth_ll,]$Value[-1]
          prob_delay_adjusted_ll_m2 <- delay_adjusted_forecasts_m2[delay_adjusted_forecasts_m2$Target==target & delay_adjusted_forecasts_m2$Location ==test_region& delay_adjusted_forecasts_m2$Bin_start_incl == truth_ll,]$Value[-1]
          prob_delay_adjusted_ll_m3 <- delay_adjusted_forecasts_m3[delay_adjusted_forecasts_m3$Target==target & delay_adjusted_forecasts_m3$Location ==test_region& delay_adjusted_forecasts_m3$Bin_start_incl == truth_ll,]$Value[-1]
          prob_delay_adjusted_ll_m4 <- delay_adjusted_forecasts_m4[delay_adjusted_forecasts_m4$Target==target & delay_adjusted_forecasts_m4$Location ==test_region& delay_adjusted_forecasts_m4$Bin_start_incl == truth_ll,]$Value[-1]
          prob_delay_adjusted_ll_m5 <- delay_adjusted_forecasts_m5[delay_adjusted_forecasts_m5$Target==target & delay_adjusted_forecasts_m5$Location ==test_region& delay_adjusted_forecasts_m5$Bin_start_incl == truth_ll,]$Value[-1]
          prob_delay_adjusted_ll_m6 <- delay_adjusted_forecasts_m6[delay_adjusted_forecasts_m6$Target==target & delay_adjusted_forecasts_m6$Location ==test_region& delay_adjusted_forecasts_m6$Bin_start_incl == truth_ll,]$Value[-1]
          prob_delay_adjusted_ll_m7 <- non_delay_adjusted_forecasts[non_delay_adjusted_forecasts$Target==target & non_delay_adjusted_forecasts$Location == test_region & non_delay_adjusted_forecasts$Bin_start_incl == min(max(round(as.numeric(truth_ll)*adjustment_factor,1),0),13),]$Value[-1]
          
          prob_non_delay_adjusted_ll <- non_delay_adjusted_forecasts[non_delay_adjusted_forecasts$Target==target & non_delay_adjusted_forecasts$Location ==test_region & non_delay_adjusted_forecasts$Bin_start_incl == truth_ll,]$Value[-1]
          prob_true_ll <- true_forecasts[true_forecasts$Target==target & true_forecasts$Location ==test_region& true_forecasts$Bin_start_incl == truth_ll,]$Value[-1]
          
          prob_delay_adjusted_r_m1 <- delay_adjusted_forecasts_m1[delay_adjusted_forecasts_m1$Target==target & delay_adjusted_forecasts_m1$Location == test_region & delay_adjusted_forecasts_m1$Bin_start_incl == truth_r,]$Value[-1]
          prob_delay_adjusted_r_m2 <- delay_adjusted_forecasts_m2[delay_adjusted_forecasts_m2$Target==target & delay_adjusted_forecasts_m2$Location == test_region & delay_adjusted_forecasts_m2$Bin_start_incl == truth_r,]$Value[-1]
          prob_delay_adjusted_r_m3 <- delay_adjusted_forecasts_m3[delay_adjusted_forecasts_m3$Target==target & delay_adjusted_forecasts_m3$Location == test_region & delay_adjusted_forecasts_m3$Bin_start_incl == truth_r,]$Value[-1]
          prob_delay_adjusted_r_m4 <- delay_adjusted_forecasts_m4[delay_adjusted_forecasts_m4$Target==target & delay_adjusted_forecasts_m4$Location == test_region & delay_adjusted_forecasts_m4$Bin_start_incl == truth_r,]$Value[-1]
          prob_delay_adjusted_r_m5 <- delay_adjusted_forecasts_m5[delay_adjusted_forecasts_m5$Target==target & delay_adjusted_forecasts_m5$Location == test_region & delay_adjusted_forecasts_m5$Bin_start_incl == truth_r,]$Value[-1]
          prob_delay_adjusted_r_m6 <- delay_adjusted_forecasts_m6[delay_adjusted_forecasts_m6$Target==target & delay_adjusted_forecasts_m6$Location == test_region & delay_adjusted_forecasts_m6$Bin_start_incl == truth_r,]$Value[-1]
          prob_delay_adjusted_r_m7 <- non_delay_adjusted_forecasts[non_delay_adjusted_forecasts$Target==target & non_delay_adjusted_forecasts$Location == test_region & non_delay_adjusted_forecasts$Bin_start_incl == min(max(round(as.numeric(truth_r)*adjustment_factor,1),0),13),]$Value[-1]
          
          
          prob_non_delay_adjusted_r <- non_delay_adjusted_forecasts[non_delay_adjusted_forecasts$Target==target & non_delay_adjusted_forecasts$Location ==test_region & non_delay_adjusted_forecasts$Bin_start_incl == truth_r,]$Value[-1]
          prob_true_r <- true_forecasts[true_forecasts$Target==target & true_forecasts$Location ==test_region& true_forecasts$Bin_start_incl == truth_r,]$Value[-1]
          
          prob_delay_adjusted_rr_m1 <- delay_adjusted_forecasts_m1[delay_adjusted_forecasts_m1$Target==target & delay_adjusted_forecasts_m1$Location == test_region & delay_adjusted_forecasts_m1$Bin_start_incl == truth_rr,]$Value[-1]
          prob_delay_adjusted_rr_m2 <- delay_adjusted_forecasts_m2[delay_adjusted_forecasts_m2$Target==target & delay_adjusted_forecasts_m2$Location == test_region & delay_adjusted_forecasts_m2$Bin_start_incl == truth_rr,]$Value[-1]
          prob_delay_adjusted_rr_m3 <- delay_adjusted_forecasts_m3[delay_adjusted_forecasts_m3$Target==target & delay_adjusted_forecasts_m3$Location == test_region & delay_adjusted_forecasts_m3$Bin_start_incl == truth_rr,]$Value[-1]
          prob_delay_adjusted_rr_m4 <- delay_adjusted_forecasts_m4[delay_adjusted_forecasts_m4$Target==target & delay_adjusted_forecasts_m4$Location == test_region & delay_adjusted_forecasts_m4$Bin_start_incl == truth_rr,]$Value[-1]
          prob_delay_adjusted_rr_m5 <- delay_adjusted_forecasts_m5[delay_adjusted_forecasts_m5$Target==target & delay_adjusted_forecasts_m5$Location == test_region & delay_adjusted_forecasts_m5$Bin_start_incl == truth_rr,]$Value[-1]
          prob_delay_adjusted_rr_m6 <- delay_adjusted_forecasts_m6[delay_adjusted_forecasts_m6$Target==target & delay_adjusted_forecasts_m6$Location == test_region & delay_adjusted_forecasts_m6$Bin_start_incl == truth_rr,]$Value[-1]
          prob_delay_adjusted_rr_m7 <- non_delay_adjusted_forecasts[non_delay_adjusted_forecasts$Target==target & non_delay_adjusted_forecasts$Location == test_region & non_delay_adjusted_forecasts$Bin_start_incl == min(max(round(as.numeric(truth_rr)*adjustment_factor,1),0),13),]$Value[-1]
          
          prob_non_delay_adjusted_rr <- non_delay_adjusted_forecasts[non_delay_adjusted_forecasts$Target==target & non_delay_adjusted_forecasts$Location == test_region & non_delay_adjusted_forecasts$Bin_start_incl == truth_rr,]$Value[-1]
          prob_true_rr <- true_forecasts[true_forecasts$Target==target& true_forecasts$Location ==test_region& true_forecasts$Bin_start_incl == truth_rr,]$Value[-1]
          
          
          if (score == "MULTIBIN"){
            prob_delay_adjusted_m1 <- sum(prob_delay_adjusted_center_m1,prob_delay_adjusted_l_m1,prob_delay_adjusted_ll_m1,prob_delay_adjusted_r_m1,prob_delay_adjusted_rr_m1)
            prob_delay_adjusted_m2 <- sum(prob_delay_adjusted_center_m2,prob_delay_adjusted_l_m2,prob_delay_adjusted_ll_m2,prob_delay_adjusted_r_m2,prob_delay_adjusted_rr_m2)
            prob_delay_adjusted_m3 <- sum(prob_delay_adjusted_center_m3,prob_delay_adjusted_l_m3,prob_delay_adjusted_ll_m3,prob_delay_adjusted_r_m3,prob_delay_adjusted_rr_m3)
            prob_delay_adjusted_m4 <- sum(prob_delay_adjusted_center_m4,prob_delay_adjusted_l_m4,prob_delay_adjusted_ll_m4,prob_delay_adjusted_r_m4,prob_delay_adjusted_rr_m4)
            prob_delay_adjusted_m5 <- sum(prob_delay_adjusted_center_m5,prob_delay_adjusted_l_m5,prob_delay_adjusted_ll_m5,prob_delay_adjusted_r_m5,prob_delay_adjusted_rr_m5)
            prob_delay_adjusted_m6 <- sum(prob_delay_adjusted_center_m6,prob_delay_adjusted_l_m6,prob_delay_adjusted_ll_m6,prob_delay_adjusted_r_m6,prob_delay_adjusted_rr_m6)
            prob_delay_adjusted_m7 <- sum(prob_delay_adjusted_center_m7,prob_delay_adjusted_l_m7,prob_delay_adjusted_ll_m7,prob_delay_adjusted_r_m7,prob_delay_adjusted_rr_m7)
            
            
            m7_plot <- c(prob_delay_adjusted_center_m7,prob_delay_adjusted_l_m7,prob_delay_adjusted_ll_m7,prob_delay_adjusted_r_m7,prob_delay_adjusted_rr_m7)
            non_delay_plot <- c(prob_non_delay_adjusted_center,prob_non_delay_adjusted_l,prob_non_delay_adjusted_ll,prob_non_delay_adjusted_r,prob_non_delay_adjusted_rr)
            truth_vec <- c(truth,truth_l,truth_ll,truth_r,truth_rr)
           # plot_df <- data.frame(true = rep(truth_vec,2),estimates = c(m7_plot,non_delay_plot),model=as.factor(c(rep("m7",5),rep("non-delay",5))))
            
          #  p <- ggplot(plot_df,aes(x=true,y=estimates,col=model))  + geom_jitter() + theme_bw() + ggtitle(paste0(test_week))
           # print (p)
            prob_non_delay_adjusted <- sum(prob_non_delay_adjusted_center,prob_non_delay_adjusted_l,prob_non_delay_adjusted_ll,prob_non_delay_adjusted_r,prob_non_delay_adjusted_rr)
            prob_true <- sum(prob_true_center,prob_true_r,prob_true_rr,prob_true_ll,prob_true_l)
          }else{
            prob_delay_adjusted_m1 <- prob_delay_adjusted_center_m1
            prob_delay_adjusted_m2 <- prob_delay_adjusted_center_m2
            prob_delay_adjusted_m3 <- prob_delay_adjusted_center_m3
            prob_delay_adjusted_m4 <- prob_delay_adjusted_center_m4
            prob_non_delay_adjusted <- prob_non_delay_adjusted_center
            prob_true <- prob_true_center
          }
          delay_adjusted_total_prob_m1 <- c(delay_adjusted_total_prob_m1,prob_delay_adjusted_m1)
          delay_adjusted_total_prob_m2 <- c(delay_adjusted_total_prob_m2,prob_delay_adjusted_m2)
          delay_adjusted_total_prob_m3 <- c(delay_adjusted_total_prob_m3,prob_delay_adjusted_m3)
          delay_adjusted_total_prob_m4 <- c(delay_adjusted_total_prob_m4,prob_delay_adjusted_m4)
          delay_adjusted_total_prob_m5 <- c(delay_adjusted_total_prob_m5,prob_delay_adjusted_m5)
          delay_adjusted_total_prob_m6 <- c(delay_adjusted_total_prob_m6,prob_delay_adjusted_m6)
          delay_adjusted_total_prob_m7 <- c(delay_adjusted_total_prob_m7,prob_delay_adjusted_m7)
          
          
          
          non_delay_adjusted_total_prob <- c(non_delay_adjusted_total_prob,prob_non_delay_adjusted)
          true_total_prob <- c(true_total_prob,prob_true)
        }
      }
    
      
    }
    
    m1_score <- gm_mean(delay_adjusted_total_prob_m1 )
    m2_score <- gm_mean(delay_adjusted_total_prob_m2 )
    m3_score <- gm_mean(delay_adjusted_total_prob_m3 )
    m4_score <- gm_mean(delay_adjusted_total_prob_m4 )
    m5_score <- gm_mean(delay_adjusted_total_prob_m5 )
    m6_score <- gm_mean(delay_adjusted_total_prob_m6 )
    m7_score <- gm_mean(delay_adjusted_total_prob_m7 )
    
    non_adjusted_score <- gm_mean(non_delay_adjusted_total_prob)
    true_score <- gm_mean(true_total_prob)
    



    data_for_plot <- rbind(data_for_plot, c(m1_score/non_adjusted_score, m2_score/non_adjusted_score,m3_score/non_adjusted_score,m4_score/non_adjusted_score,m5_score/non_adjusted_score,m6_score/non_adjusted_score,m7_score/non_adjusted_score,true_r=true_score/non_adjusted_score,target,test_season))
  }
}


library(ggplot2)
data_for_plot <- as.data.frame(data_for_plot[2:nrow(data_for_plot),])
colnames(data_for_plot) <- c("Empirical Mean","Sampling","Previous Lag","NN","LM","HLM","post", "true","target","season")

library(reshape2)
library(dplyr)
long <- melt(data_for_plot, id.vars = c("target","season"))
long$value <- as.numeric(as.character(long$value))
long <- long %>%
  group_by(target,variable) %>%
  dplyr::summarize(value = mean(value, na.rm=TRUE))
p <- ggplot(long,aes(x=target,y=value,col=variable)) + geom_point() + theme_bw()+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_abline(intercept = 1,slope = 0,alpha=.2)  
p
p <- ggplot(long,aes(x=target,y=value,col=variable)) + geom_point() + theme_bw()+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_abline(intercept = 1,slope = 0,alpha=.2)  + ylim(.9,1.1)
p
