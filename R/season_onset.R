
data <- readRDS("./data/flu_data_with_backfill_edit.rds")
fully_observed_data <- data.frame(readRDS("./data/fully_observed_data_formatted.rds"))

test_regions <- c("National")

for (test_region in test_regions){
  for (s in c("2010/2011")){
    peak_week <- which.max(fully_observed_data[fully_observed_data$season == s & fully_observed_data$region ==test_region ,]$weighted_ili)
    if (peak_week > 12){
      peak_week_formatted <- peak_week - 22
    } else{
      peak_week_formatted <- peak_week+ 31
    }
    print (peak_week_formatted)
  }
}
