library(EnvStats)
library(ggplot2)
fully_observed_data <- read.csv("fully_observed_data")
#model_params <- read.csv("model_params.csv")


#test_region <-toString(model_params$region)
#model_var <- model_params$model_variance
step_ahead <- 1
score <- "MULTIBIN"


delay_adjusted_total_prob_m1 <- c()
delay_adjusted_total_prob_m2 <- c()
non_delay_adjusted_total_prob <- c()
true_total_prob <- c()


delay_adjusted_total_bias <- c()
non_delay_adjusted_total_bias <- c()
true_total_bias <- c()
for (test_region in c("nat") ){
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
      
      
      delay_adjusted_forecasts_m2 <- read.csv(paste0("./inst/submissions/region-sarima_seasonal_difference_TRUE/EW",test_week_formatted,"-2018","-ReichLab_sarima_seasonal_difference_TRUE-delay-M1.csv"))
      #delay_adjusted_forecasts_m2 <- read.csv(paste0("/Users/gcgibson/2018-2019-cdc-flu-contest/inst/submissions/region-sarima_seasonal_difference_TRUE/EW",test_region,"-",test_season_formatted,"-",test_week_formatted,"-","M2","-",model_var))
      non_delay_adjusted_forecasts <- read.csv(paste0("./inst/submissions/region-sarima_seasonal_difference_TRUE/EW",test_week_formatted,"-2018","-ReichLab_sarima_seasonal_difference_TRUE-delay-NONE.csv"))
      #true_forecasts <- read.csv(paste0("output/",test_region,"-",test_season_formatted,"-",test_week_formatted,"-","TRUTH","-",model_var))
      delay_adjusted_forecasts_m1 <- delay_adjusted_forecasts_m2
      true_forecasts <- delay_adjusted_forecasts_m2
      
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
      truth <- get_inc_bin(fully_observed_data[fully_observed_data$region == test_region &fully_observed_data$epiweek == truth_time,]$wili)
      truth_l <- max(0,as.numeric(truth)-.1)
      truth_ll <- max(0,as.numeric(truth)-.2)
      truth_r <- min(13,as.numeric(truth)+.1)
      truth_rr <- min(13,as.numeric(truth)+.2)
      
      
      plot_df <- data.frame(delay=delay_adjusted_forecasts_m2[delay_adjusted_forecasts_m2$Target=="1 wk ahead" &delay_adjusted_forecasts_m2$Location=="US National" ,]$Value[2:25], original=non_delay_adjusted_forecasts[non_delay_adjusted_forecasts$Target=="1 wk ahead" & non_delay_adjusted_forecasts$Location=="US National",]$Value[2:25])
      p <- ggplot(plot_df,aes(x=1:length(delay),y=delay,col="delay")) + geom_line() + geom_line(aes(x=1:length(original),y=original,col="original",alpha=.1)) + ylab("P(wili in bin i)") + xlab("bin index") + ggtitle(paste0(test_season_formatted,"-",test_week_formatted ,"- ", step_ahead,"  step ahead"))
      p <- p + geom_vline(xintercept=as.double(truth)*10,linetype="dotted") + theme_bw()
      print (p)
      
      prob_delay_adjusted_center_m1 <- delay_adjusted_forecasts_m1[delay_adjusted_forecasts_m1$Target=="1 wk ahead" & delay_adjusted_forecasts_m1$Location == "US National" & delay_adjusted_forecasts_m1$Bin_start_incl == truth,]$Value[-1]
      prob_delay_adjusted_center_m2 <- delay_adjusted_forecasts_m2[delay_adjusted_forecasts_m2$Target=="1 wk ahead" & delay_adjusted_forecasts_m2$Location == "US National" & delay_adjusted_forecasts_m2$Bin_start_incl == truth,]$Value[-1]
      prob_non_delay_adjusted_center <- non_delay_adjusted_forecasts[non_delay_adjusted_forecasts$Target=="1 wk ahead" & non_delay_adjusted_forecasts$Location == "US National" & non_delay_adjusted_forecasts$Bin_start_incl == truth,]$Value[-1]
      prob_true_center <- true_forecasts[true_forecasts$X == truth,paste0("V",step_ahead)]
      
      prob_delay_adjusted_l_m1 <- delay_adjusted_forecasts_m1[delay_adjusted_forecasts_m1$Target=="1 wk ahead" & delay_adjusted_forecasts_m1$Location == "US National" & delay_adjusted_forecasts_m1$Bin_start_incl == truth_l,]$Value[-1]
      prob_delay_adjusted_l_m2 <- delay_adjusted_forecasts_m2[delay_adjusted_forecasts_m2$Target=="1 wk ahead" & delay_adjusted_forecasts_m2$Location == "US National" & delay_adjusted_forecasts_m2$Bin_start_incl == truth_l,]$Value[-1]
      prob_non_delay_adjusted_l <- non_delay_adjusted_forecasts[non_delay_adjusted_forecasts$Target=="1 wk ahead" & non_delay_adjusted_forecasts$Location == "US National" & non_delay_adjusted_forecasts$Bin_start_incl == truth_l,]$Value[-1]
      prob_true_l <- true_forecasts[true_forecasts$X == truth,paste0("V",step_ahead)]
      
      prob_delay_adjusted_ll_m1 <- delay_adjusted_forecasts_m1[delay_adjusted_forecasts_m1$Target=="1 wk ahead" & delay_adjusted_forecasts_m1$Location == "US National" & delay_adjusted_forecasts_m1$Bin_start_incl == truth_ll,]$Value[-1]
      prob_delay_adjusted_ll_m2 <- delay_adjusted_forecasts_m2[delay_adjusted_forecasts_m2$Target=="1 wk ahead" & delay_adjusted_forecasts_m2$Location == "US National" & delay_adjusted_forecasts_m2$Bin_start_incl == truth_ll,]$Value[-1]
      prob_non_delay_adjusted_ll <- non_delay_adjusted_forecasts[non_delay_adjusted_forecasts$Target=="1 wk ahead" & non_delay_adjusted_forecasts$Location == "US National" & non_delay_adjusted_forecasts$Bin_start_incl == truth_ll,]$Value[-1]
      prob_true_ll <- true_forecasts[true_forecasts$X == truth,paste0("V",step_ahead)]
      
      prob_delay_adjusted_r_m1 <- delay_adjusted_forecasts_m1[delay_adjusted_forecasts_m1$Target=="1 wk ahead" & delay_adjusted_forecasts_m1$Location == "US National" & delay_adjusted_forecasts_m1$Bin_start_incl == truth_r,]$Value[-1]
      prob_delay_adjusted_r_m2 <- delay_adjusted_forecasts_m2[delay_adjusted_forecasts_m2$Target=="1 wk ahead" & delay_adjusted_forecasts_m2$Location == "US National" & delay_adjusted_forecasts_m2$Bin_start_incl == truth_r,]$Value[-1]
      prob_non_delay_adjusted_r <- non_delay_adjusted_forecasts[non_delay_adjusted_forecasts$Target=="1 wk ahead" & non_delay_adjusted_forecasts$Location == "US National" & non_delay_adjusted_forecasts$Bin_start_incl == truth_r,]$Value[-1]
      prob_true_r <- true_forecasts[true_forecasts$X == truth,paste0("V",step_ahead)]
      
      prob_delay_adjusted_rr_m1 <- delay_adjusted_forecasts_m1[delay_adjusted_forecasts_m1$Target=="1 wk ahead" & delay_adjusted_forecasts_m1$Location == "US National" & delay_adjusted_forecasts_m1$Bin_start_incl == truth_rr,]$Value[-1]
      prob_delay_adjusted_rr_m2 <- delay_adjusted_forecasts_m2[delay_adjusted_forecasts_m2$Target=="1 wk ahead" & delay_adjusted_forecasts_m2$Location == "US National" & delay_adjusted_forecasts_m2$Bin_start_incl == truth_rr,]$Value[-1]
      prob_non_delay_adjusted_rr <- non_delay_adjusted_forecasts[non_delay_adjusted_forecasts$Target=="1 wk ahead" & non_delay_adjusted_forecasts$Location == "US National" & non_delay_adjusted_forecasts$Bin_start_incl == truth_rr,]$Value[-1]
      prob_true_rr <- true_forecasts[true_forecasts$X == truth,paste0("V",step_ahead)]
      
      
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
      delay_adjusted_total_prob_m1 <- c(delay_adjusted_total_prob_m2,prob_delay_adjusted_m1)
      delay_adjusted_total_prob_m2 <- c(delay_adjusted_total_prob_m2,prob_delay_adjusted_m2)
      non_delay_adjusted_total_prob <- c(non_delay_adjusted_total_prob,prob_non_delay_adjusted)
      true_total_prob <- c(true_total_prob,prob_true)
    }
  }
}
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
print (c(paste0(step_ahead, " step ahead "),score))
print ("Delay Adjusted - Model 1 (Mean)")
print (log(gm_mean(delay_adjusted_total_prob_m1 +.000000000000000000001)))
print ("Delay Adjusted - Model 2 (Sampling)")
print (log(gm_mean(delay_adjusted_total_prob_m2 +.000000000000000000001)))
print ("Non-delay Adjusted")
print (log(gm_mean(non_delay_adjusted_total_prob+.000000000000000000001)))
print ("True")
print (log(gm_mean(true_total_prob+.000000000000000000001)))

