data <- readRDS("./data/flu_data_with_backfill_edit.rds")

extreme_lag_df <- matrix(NA,ncol=7)
for (region in unique(data$region)){
  for (epiweek  in unique(data[data$region == region,]$epiweek)){
     l0 <- data[data$region == region & data$epiweek==epiweek & data$lag == 0, ]$weighted_ili
     linfinity <- data[data$region == region & data$epiweek==epiweek & data$lag == max(data[data$region == region & data$epiweek==epiweek,]$lag), ]$weighted_ili
      
     out <- tryCatch({
      if( abs(l0-linfinity)  > .5  ){
        previous_data_point <- data[data$region == region & data$epiweek==epiweek-1 & data$lag == 1, ]$weighted_ili
        extreme_lag_df<- rbind(extreme_lag_df,c(region,epiweek,l0-linfinity,l0,previous_data_point,linfinity,l0-previous_data_point))
      }
       },error=function(cond){
        
      }
     )
      
     
     
  }
}

extreme_lag_df <- data.frame(extreme_lag_df[2:nrow(extreme_lag_df),])
colnames(extreme_lag_df) <- c("Region","Epiweek","reporting_ratio","l0","previous","linfty","l0/previos")
extreme_lag_df_fromatted <- cbind.data.frame(extreme_lag_df[,1],sapply( extreme_lag_df[,2:ncol(extreme_lag_df)],function(x){ as.numeric(as.character(x))} ))
extreme_lag_df_fromatted <- data.frame(extreme_lag_df_fromatted)
colnames(extreme_lag_df_fromatted) <- c("Region","Epiweek","reporting_ratio","l0","previous","linfty","l0/previos")
