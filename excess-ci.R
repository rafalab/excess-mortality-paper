### -- ------ ------------------------------------------------------------------
### -- Set up ------------------------------------------------------------------
### -- ------ ------------------------------------------------------------------
# -- Libraries
library(tidyverse)
library(lubridate)
library(excessmort)
dslabs::ds_theme_set()

source("pr-init.R")
### -- ------------------------------ ------------------------------------------------------------------
### -- Figure 2A: Excess deaths in PR ------------------------------------------------------------------
### -- ------------------------------ ------------------------------------------------------------------
# -- Set up to be used below
knots  <- rep(4, 3)
disc   <- c(TRUE, TRUE, FALSE)
before <- c(365, 365, 365) 
after  <- c(365, 365, 365) * 2

interval_start <- c(hurricane_dates[2], # Georges
                    hurricane_dates[3], # Maria
                    Chikungunya = make_date(2014,08,01))
                    

# Get the fits ------------------------------------------------------------
fs <- map_df(seq_along(interval_start), function(i){
  cat("\nEvent:", names(interval_start)[i], "\n")
  # -- Fitting the model to each age group
  message("Fitting model to estimate period of effect")
  groups <- unique(all_counts$agegroup)
  tmp <- map_df(groups, function(x){
    cat(".")
    ev <- NULL
    if(disc[i]) ev <- interval_start[i]
    f <- all_counts %>% 
      filter(agegroup == x) %>%
      excess_model(event          = ev,
                   start          = interval_start[i] - before[i],
                   end            = interval_start[i] + after[i], 
                   exclude        = exclude_dates,
                   control.dates  = control_dates,
                   knots.per.year = knots[i],
                   weekday.effect = TRUE,
                   model          = "correlated",
                   discontinuity  = disc[i], 
                   verbose = FALSE)
    tibble(date = f$date, expected = f$expected, fitted = f$fitted, se = f$se, agegroup = x)
  })
  
  # -- Computing marginal effect
  fit <- tmp %>% 
    group_by(date) %>% 
    summarize(fitted = sum(expected * fitted) / sum(expected), 
              se     = sqrt(sum(expected^2 * se^2)) / sum(expected)) %>%
    ungroup() %>% 
    mutate(event = names(interval_start)[i])
  fit
})

## eye test
fs %>% 
  ggplot(aes(date,fitted)) + 
  geom_ribbon(aes(ymin = fitted - 2*se, ymax = fitted + 2*se), alpha = 0.25)+
  geom_line() + 
  facet_wrap(~event, scales = "free")
  
## Look for 0s and max
tmp <- fs %>% 
  group_by(event) %>%
  mutate(max = fitted == max(fitted),
         max_date = date[which.max(fitted)], 
         days_to_max = as.numeric(date - max_date),
         sign_change = abs(sign(fitted) != lag(sign(fitted)))>0) %>%
  ungroup() %>%
  filter(sign_change | max) 

  
## after goerges, first day crossing 0
tmp %>% 
  filter(event == "Georges" & days_to_max > 0) %>%
  top_n(1, -days_to_max)  

## after maria, first day crossing 0
tmp %>% 
  filter(event == "Maria" & days_to_max > 0) %>%
  top_n(1, -days_to_max)  

## Chikungunya start, max, start of decrease, end
tmp %>% 
  filter(event == "Chikungunya") %>%
  top_n(5, -abs(days_to_max)) %>%
  slice(-1)

f