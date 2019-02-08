library(EnvStats)
library(ggplot2)
fully_observed_data <- readRDS("./data/fully_observed_data_formatted.rds")
#model_params <- read.csv("model_params.csv")

region_str_array_eval <- c("National",paste0("Region ",1:10))
region_str_data_set <- c("US National","HHS Region 1",  "HHS Region 2",  "HHS Region 3",  "HHS Region 4" ,"HHS Region 5",  "HHS Region 6",  "HHS Region 7",  "HHS Region 8",  "HHS Region 9","HHS Region 10")
region_str_true <- c("nat",paste0("hhs",1:10))
#test_region <-toString(model_params$region)
#model_var <- model_params$model_variance
target <- "Season onset"
step_ahead <- 1
score <- "MULTIBIN"
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}




delay_adjusted_total_prob_m1 <- c()
delay_adjusted_total_prob_m2 <- c()
non_delay_adjusted_total_prob <- c()

for (test_region in region_str_data_set ){

  true_total_prob <- c()
  for (test_season in c("2015")){
    if (test_season == "2017"){
      end_week <- 12
    }else{
      end_week <- 20
    }
    for (test_week in c(seq(40,52),seq(end_week))){
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
      
      
      delay_adjusted_forecasts_m2 <- read.csv(paste0("./inst/submissions/region-sarima_seasonal_difference_TRUE/EW",test_week_formatted,"-",test_season_formatted,"-ReichLab_sarima_seasonal_difference_TRUE-delay-M3.csv"))
      delay_adjusted_forecasts_m1 <- read.csv(paste0("./inst/submissions/region-sarima_seasonal_difference_TRUE/EW",test_week_formatted,"-",test_season_formatted,"-ReichLab_sarima_seasonal_difference_TRUE-delay-M4.csv"))
      non_delay_adjusted_forecasts <- read.csv(paste0("./inst/submissions/region-sarima_seasonal_difference_TRUE/EW",test_week_formatted,"-",test_season_formatted,"-ReichLab_sarima_seasonal_difference_TRUE-delay-NONE.csv"))
      true_forecasts <- read.csv(paste0("./inst/submissions/region-sarima_seasonal_difference_TRUE/EW",test_week_formatted,"-",test_season_formatted,"-ReichLab_sarima_seasonal_difference_TRUE-delay-TRUE.csv"))
      #delay_adjusted_forecasts_m1 <- delay_adjusted_forecasts_m2
      truth_time <- NULL
      composite_forecast <- non_delay_adjusted_forecasts
      
      composite_forecast$Value <- .2*delay_adjusted_forecasts_m1$Value + .8*non_delay_adjusted_forecasts$Value
      
      delay_adjusted_forecasts_m2 <- composite_forecast
      
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
      if (target != "Season onset"){
        truth <- get_inc_bin(fully_observed_data[fully_observed_data$region == region_str_array_eval[match(test_region,region_str_data_set)] &fully_observed_data$epiweek == truth_time,]$wili)
        truth_l <- max(0,as.numeric(truth)-.1)
        truth_ll <- max(0,as.numeric(truth)-.2)
        truth_r <- min(13,as.numeric(truth)+.1)
        truth_rr <- min(13,as.numeric(truth)+.2)
      } else{
        baseline <- get_onset_baseline(region = "National", season = "2016/2017")
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
        
      }
      
      
      
      
      plot_df <- data.frame(true =true_forecasts[true_forecasts$Target==target &true_forecasts$Location==test_region ,]$Value[2:40], delay=delay_adjusted_forecasts_m2[delay_adjusted_forecasts_m2$Target==target &delay_adjusted_forecasts_m2$Location==test_region ,]$Value[2:40], original=non_delay_adjusted_forecasts[non_delay_adjusted_forecasts$Target==target & non_delay_adjusted_forecasts$Location==test_region,]$Value[2:40])
      p <- ggplot(plot_df,aes(x=1:length(delay),y=delay,col="delay")) + geom_line() + geom_line(aes(x=1:length(original),y=original,col="original"),alpha=.5) + ylab("P(wili in bin i)") + xlab("bin index") + ggtitle(paste0(test_season_formatted,"-",test_week_formatted ,"- ", step_ahead,"  step ahead"))
      p <- p + geom_vline(xintercept=as.double(truth)*10,linetype="dotted") + theme_bw()
      p <- p + geom_line(aes(x=1:length(true),y=true,col="true"),alpha=.5) 
      #print (p)
      
      prob_delay_adjusted_center_m1 <- delay_adjusted_forecasts_m1[delay_adjusted_forecasts_m1$Target==target & delay_adjusted_forecasts_m1$Location == test_region& delay_adjusted_forecasts_m1$Bin_start_incl == truth,]$Value[-1]
      prob_delay_adjusted_center_m2 <- delay_adjusted_forecasts_m2[delay_adjusted_forecasts_m2$Target==target & delay_adjusted_forecasts_m2$Location == test_region & delay_adjusted_forecasts_m2$Bin_start_incl == truth,]$Value[-1]
      prob_non_delay_adjusted_center <- non_delay_adjusted_forecasts[non_delay_adjusted_forecasts$Target==target & non_delay_adjusted_forecasts$Location ==test_region& non_delay_adjusted_forecasts$Bin_start_incl == truth,]$Value[-1]
      prob_true_center <- true_forecasts[true_forecasts$Target==target & true_forecasts$Location ==test_region& true_forecasts$Bin_start_incl == truth,]$Value[-1]
      
      prob_delay_adjusted_l_m1 <- delay_adjusted_forecasts_m1[delay_adjusted_forecasts_m1$Target==target & delay_adjusted_forecasts_m1$Location == test_region & delay_adjusted_forecasts_m1$Bin_start_incl == truth_l,]$Value[-1]
      prob_delay_adjusted_l_m2 <- delay_adjusted_forecasts_m2[delay_adjusted_forecasts_m2$Target==target & delay_adjusted_forecasts_m2$Location == test_region& delay_adjusted_forecasts_m2$Bin_start_incl == truth_l,]$Value[-1]
      prob_non_delay_adjusted_l <- non_delay_adjusted_forecasts[non_delay_adjusted_forecasts$Target==target & non_delay_adjusted_forecasts$Location ==test_region& non_delay_adjusted_forecasts$Bin_start_incl == truth_l,]$Value[-1]
      prob_true_l <- true_forecasts[true_forecasts$Target==target & true_forecasts$Location ==test_region& true_forecasts$Bin_start_incl == truth_l,]$Value[-1]
      
      prob_delay_adjusted_ll_m1 <- delay_adjusted_forecasts_m1[delay_adjusted_forecasts_m1$Target==target & delay_adjusted_forecasts_m1$Location == test_region& delay_adjusted_forecasts_m1$Bin_start_incl == truth_ll,]$Value[-1]
      prob_delay_adjusted_ll_m2 <- delay_adjusted_forecasts_m2[delay_adjusted_forecasts_m2$Target==target & delay_adjusted_forecasts_m2$Location ==test_region& delay_adjusted_forecasts_m2$Bin_start_incl == truth_ll,]$Value[-1]
      prob_non_delay_adjusted_ll <- non_delay_adjusted_forecasts[non_delay_adjusted_forecasts$Target==target & non_delay_adjusted_forecasts$Location ==test_region & non_delay_adjusted_forecasts$Bin_start_incl == truth_ll,]$Value[-1]
      prob_true_ll <- true_forecasts[true_forecasts$Target==target & true_forecasts$Location ==test_region& true_forecasts$Bin_start_incl == truth_ll,]$Value[-1]
      
      prob_delay_adjusted_r_m1 <- delay_adjusted_forecasts_m1[delay_adjusted_forecasts_m1$Target==target & delay_adjusted_forecasts_m1$Location == test_region & delay_adjusted_forecasts_m1$Bin_start_incl == truth_r,]$Value[-1]
      prob_delay_adjusted_r_m2 <- delay_adjusted_forecasts_m2[delay_adjusted_forecasts_m2$Target==target & delay_adjusted_forecasts_m2$Location == test_region & delay_adjusted_forecasts_m2$Bin_start_incl == truth_r,]$Value[-1]
      prob_non_delay_adjusted_r <- non_delay_adjusted_forecasts[non_delay_adjusted_forecasts$Target==target & non_delay_adjusted_forecasts$Location ==test_region & non_delay_adjusted_forecasts$Bin_start_incl == truth_r,]$Value[-1]
      prob_true_r <- true_forecasts[true_forecasts$Target==target & true_forecasts$Location ==test_region& true_forecasts$Bin_start_incl == truth_r,]$Value[-1]
      
      prob_delay_adjusted_rr_m1 <- delay_adjusted_forecasts_m1[delay_adjusted_forecasts_m1$Target==target & delay_adjusted_forecasts_m1$Location == test_region & delay_adjusted_forecasts_m1$Bin_start_incl == truth_rr,]$Value[-1]
      prob_delay_adjusted_rr_m2 <- delay_adjusted_forecasts_m2[delay_adjusted_forecasts_m2$Target==target & delay_adjusted_forecasts_m2$Location == test_region & delay_adjusted_forecasts_m2$Bin_start_incl == truth_rr,]$Value[-1]
      prob_non_delay_adjusted_rr <- non_delay_adjusted_forecasts[non_delay_adjusted_forecasts$Target==target & non_delay_adjusted_forecasts$Location == test_region & non_delay_adjusted_forecasts$Bin_start_incl == truth_rr,]$Value[-1]
      prob_true_rr <- true_forecasts[true_forecasts$Target==target& true_forecasts$Location ==test_region& true_forecasts$Bin_start_incl == truth_rr,]$Value[-1]
      
      
      if (score == "MULTIBIN"){
        prob_delay_adjusted_m1 <- sum(prob_delay_adjusted_center_m1,prob_delay_adjusted_l_m1,prob_delay_adjusted_ll_m1,prob_delay_adjusted_r_m1,prob_delay_adjusted_rr_m1)
        prob_delay_adjusted_m2 <- sum(prob_delay_adjusted_center_m2,prob_delay_adjusted_l_m2,prob_delay_adjusted_ll_m2,prob_delay_adjusted_r_m2,prob_delay_adjusted_rr_m2)
        prob_non_delay_adjusted <- sum(prob_non_delay_adjusted_center,prob_non_delay_adjusted_l,prob_non_delay_adjusted_ll,prob_non_delay_adjusted_r,prob_non_delay_adjusted_rr)
        prob_true <- sum(prob_true_center,prob_true_r,prob_true_rr,prob_true_ll,prob_true_l)
      }else{
        prob_delay_adjusted <- prob_delay_adjusted_center
        prob_non_delay_adjusted <- prob_non_delay_adjusted_center
        prob_true <- prob_true_center
      }
      delay_adjusted_total_prob_m1 <- c(delay_adjusted_total_prob_m1,prob_delay_adjusted_m1)
      delay_adjusted_total_prob_m2 <- c(delay_adjusted_total_prob_m2,prob_delay_adjusted_m2)
      non_delay_adjusted_total_prob <- c(non_delay_adjusted_total_prob,prob_non_delay_adjusted)
      true_total_prob <- c(true_total_prob,prob_true)
    }
  }

  
}

m1_score <- gm_mean(delay_adjusted_total_prob_m1 )
m2_score <- gm_mean(delay_adjusted_total_prob_m2 )
non_adjusted_score <- gm_mean(non_delay_adjusted_total_prob)
true_score <- gm_mean(true_total_prob)



library(ggplot2)

data_for_plot <- data.frame(scores =c( m1_score/non_adjusted_score, m2_score/non_adjusted_score,true_r=true_score/non_adjusted_score),models=c("M1","M2","TRUE"))
ggplot(data_for_plot,aes(x=models,y=scores)) + geom_point() +theme_bw()
ggplot(data_for_plot[1:2,],aes(x=models,y=scores)) + geom_point() +theme_bw()
