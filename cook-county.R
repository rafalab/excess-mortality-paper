library(tidyverse)
library(lubridate)
library(excessdeaths)
library(directlabels)
dslabs::ds_theme_set()
data("cook_records")

flu_season <- seq(make_date(2017, 12, 16), make_date(2018, 1, 16), by = "day")
exclude_dates <- c(flu_season, seq(make_date(2020, 1, 1), max(cdc_state_counts$date, na.rm = TRUE), by = "day"))

the_breaks <- c(0, 20, 40, 60, 75, Inf)
all_counts <- filter(cook_records, race %in% c("white", "black")) %>%
  compute_counts(demo = cook_demographics, by = c("agegroup", "race"), breaks = the_breaks)

tmp <- expand.grid(race = unique(all_counts$race), agegroup = unique(all_counts$agegroup))
fits <- lapply(1:nrow(tmp), function(i){
  print(tmp[i,])
  res <- all_counts %>% filter(race == tmp$race[i] & agegroup == tmp$agegroup[i]) %>%
    excess_model(exclude = exclude_dates, start = min(all_counts$date), end = max(all_counts$date)) 
  res$race <- tmp$race[i]
  res$agegroup <- tmp$agegroup[i]
  return(res)
})

tmp <- map_df(fits, function(f){
  with(f, tibble(race = race,
                 agegroup = agegroup,
                 date = date, 
                 expected =  expected, 
                 fitted = 100* fitted, 
                 se = 100* se,
                 sd = 100* sd))
})


# Figur 3 - difference in groups ------------------------------------------

tmp %>%
  filter(year(date) == 2020 & agegroup %in% c("40-59", "60-74", "75-Inf")) %>%
  ggplot(aes(date, fitted, color = race, fill = race)) +
  geom_ribbon(aes(ymin = fitted - 2 * se, ymax = fitted + 2*se), alpha = 0.5) +
  geom_line() +
  facet_wrap( ~ agegroup)
  