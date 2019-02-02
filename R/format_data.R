data <- readRDS("data/flu_data_with_backfill.rds")
data$year <- unlist(lapply(data$epiweek,function(x){substr(x,1,4)}))
data$week<- unlist(lapply(data$epiweek,function(x){substr(x,5,7)}))
for (row in 1:nrow(data)){
  if (data[row,]$week >20){
      current_year <- data[row,]$year
      current_year_1 <- as.numeric(current_year) + 1
      data[row,"season"] <- paste0(current_year,"/",current_year_1)
  }else{
    current_year <- data[row,]$year
    current_year_1 <- as.numeric(current_year) - 1
    data[row,"season"] <- paste0(current_year_1,"/",current_year)
  }
}
levels(data$region) <- c("National",paste0("Region ",1:10))


week_to_season_week <- function(w){
  w <- as.numeric(w)
  if (sum(w==40:53) ==1){
    return (w-30)
  }else if(sum(w==31:39)==1){
      return (w-30)  
  }else {
    return (w+22)
  }
}

data$season_week <- unlist(lapply(data$week,week_to_season_week))
data$weighted_ili <- data$wili
data$week <- as.numeric(data$week)
saveRDS(data,"data/flu_data_with_backfill_edit.rds")


