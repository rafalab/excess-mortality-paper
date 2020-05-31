### -- ------ ------------------------------------------------------------------
### -- Set up ------------------------------------------------------------------
### -- ------ ------------------------------------------------------------------
# -- Libraries
library(scales)
library(ggpubr)
library(tidyverse)
library(lubridate)
library(excessmort)
library(directlabels)
dslabs::ds_theme_set()

# -- Hurricanes information
hurricane_dates        <- as.Date(c("1989-09-18","1998-09-21","2017-09-20"))
hurricane_effect_ends  <- as.Date(c("1990-03-18","1999-03-21","2018-03-20"))
names(hurricane_dates) <- c("Hugo", "Georges", "Maria")

# -- Exclude dates
exclude_dates <- c(seq(hurricane_dates[1], hurricane_effect_ends[1], by = "day"),
                   seq(hurricane_dates[2], hurricane_effect_ends[2], by = "day"),
                   seq(hurricane_dates[3], hurricane_effect_ends[3], by = "day"),
                   seq(as.Date("2004-12-01"), as.Date("2015-03-21"), by = "day"),
                   seq(as.Date("2004-12-01"), as.Date("2005-12-01"), by = "day"),
                   seq(as.Date("2001-01-01"), as.Date("2001-01-15"), by = "day"),
                   seq(as.Date("2020-01-01"), lubridate::today(), by = "day"))

# -- Control dates
control_dates <- seq(as.Date("2006-01-01"), as.Date("2013-12-31"), by = "day")

# -- Loading PR data and creating age groups
data("puerto_rico_counts")
the_breaks <- c(0, 5, 20, 40, 60, 75, Inf)
all_counts <- collapse_counts_by_age(puerto_rico_counts, the_breaks)
### -- ---------- ------------------------------------------------------------------
### -- END Set up ------------------------------------------------------------------
### -- ---------- ------------------------------------------------------------------

### -- ------------------------------ ------------------------------------------------------------------
### -- Figure 2A: Excess deaths in PR ------------------------------------------------------------------
### -- ------------------------------ ------------------------------------------------------------------
# -- Set up to be used below
knots  <- rep(12, 5)
disc   <- c(TRUE, TRUE, FALSE, FALSE, FALSE)
before <- c(365, 365, 365, 365, 365) 
after  <- c(365, 365, 365, 365, 70)

interval_start <- c(hurricane_dates[2], # Georges
                    hurricane_dates[3], # Maria
                    "Flu-2005"  = make_date(2004,11,01),
                    Chikungunya = make_date(2014,08,01),
                    "Covid-19"  = make_date(2020,01,01))


# Get the fits ------------------------------------------------------------
fs <- map_df(seq_along(interval_start), function(i){
  cat("\nEvent:", names(interval_start)[i], "\n")
  # -- Fitting the model to each age group
  message("Fitting model to estimate period of effect")
  if(names(interval_start)[i] %in% c("Georges", "Maria")) groups <- unique(all_counts$agegroup) else groups <- unique(all_counts$agegroup)[5:6]
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

fs %>% ggplot(aes(date,fitted))+geom_line() + facet_wrap(~event, scales = "free")
  
## Look for 0s and max
tmp <- fs %>% group_by(event) %>%
  mutate(max = fitted == max(fitted),
         max_date = date[which.max(fitted)], days_to_max = as.numeric(date - max_date),
  sign_change = abs(sign(fitted) != lag(sign(fitted)))>0) %>%
  ungroup %>%
  filter(sign_change | max) 
  
## after goerges, first day crossing 0
tmp %>% filter(event == "Georges" & days_to_max > 0) %>%
  top_n(1, -days_to_max)  

## after maria, first day crossing 0
tmp %>% filter(event == "Maria" & days_to_max > 0) %>%
  top_n(1, -days_to_max)  

## flu 2005 start max, and end of bump
tmp %>% filter(event == "Flu-2005") %>%
  top_n(3, -abs(days_to_max))

## Chikungunya start, max, start of decrease, end
tmp %>% filter(event == "Chikungunya") %>%
  top_n(6, -abs(days_to_max)) %>% 
  slice(-c(1:2))

