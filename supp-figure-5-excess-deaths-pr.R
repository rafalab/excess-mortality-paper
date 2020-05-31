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
                   seq(as.Date("2004-09-01"), as.Date("2005-12-31"), by = "day"),
                   seq(as.Date("2014-09-01"), as.Date("2015-03-21"), by = "day"),
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

### -- ---------------------------------- ------------------------------------------------------------------
### -- Supp figure 5: Excess deaths in PR ------------------------------------------------------------------
### -- ---------------------------------- ------------------------------------------------------------------
# -- Set up to be used below
ndays  <- 365
knots  <- c(4, 4, 4, 4, 4)
disc   <- c(TRUE, TRUE, FALSE, FALSE, FALSE)
before <- c(365, 365, 365, 365, 365) 
after  <- c(365, 365, 365, 365, 105)
interval_start <- c(hurricane_dates[2], # Georges
                    hurricane_dates[3], # Maria
                    Chikungunya = make_date(2014,08,01),
                    "Flu-2005"  = make_date(2004,11,01),
                    "Covid-19"  = make_date(2020,01,01))

# -- Outer loop that goes through events in PR
excess_deaths_pr <- map_df(seq_along(interval_start), function(i)
{
  cat("\nEvent:", names(interval_start)[i], "\n")
  
  # -- Fitting the model to each age group
  message("Fitting model to estimate period of effect")
  tmp <- map_df(unique(all_counts$agegroup), function(x){
    cat(".")
    f <- suppressMessages(all_counts %>% 
                            filter(agegroup == x) %>%
                            excess_model(event          = interval_start[i],
                                         start          = interval_start[i] - before[i],
                                         end            = interval_start[i] + after[i], 
                                         exclude        = exclude_dates,
                                         control.dates  = control_dates,
                                         knots.per.year = knots[i],
                                         weekday.effect = TRUE,
                                         model          = "correlated",
                                         discontinuity  = disc[i], 
                                         verbose = FALSE))
    
    tibble(date = f$date, expected = f$expected, fitted = f$fitted, se = f$se, agegroup = x)
  })
  
  # -- Computing marginal effect
  fit <- tmp %>% 
    group_by(date) %>% 
    summarize(fitted = sum(expected * fitted) / sum(expected), 
              se     = sqrt(sum(expected^2 * se^2)) / sum(expected)) %>%
    ungroup()
  
  # -- Determine period of effect
  # dates <- filter(fit, date >= interval_start[i], fitted-1.96*se >= 0)$date
  dates <- filter(fit, date >= interval_start[i], fitted >= 0)$date
  ind   <- suppressWarnings(min(which(diff(dates) > 1)))
  
  # -- Find last day
  if(is.infinite(ind)) {
    last_day <- ymd(last(dates))
    if(is.na(last_day)) {last_day <- interval_start[i]+1}
  } else {
    last_day <- ymd(dates[ind])
  }
  
  p <- fit %>%
    ggplot(aes(date, fitted)) +
    geom_hline(yintercept = 0, color="gray") +
    geom_ribbon(aes(ymin=fitted-1.96*se, ymax=fitted+1.96*se), alpha=0.50) +
    geom_line() +
    ggtitle(names(interval_start)[i]) +
    geom_vline(xintercept = last_day, color="red")
  print(p)
  
  # -- Period of indirect effect
  message(paste0("\nPeriod of indirect effect: ",last_day - interval_start[i], " days"))
  
  # -- Now fit model to compute cumulative excess deaths
  message("\nFitting model to estimate cumulative excess deaths")
  tmp <- map_df(unique(all_counts$agegroup), function(x){
    cat(".")
    f <- suppressMessages(all_counts %>% 
                            filter(agegroup == x) %>%
                            excess_model(event          = interval_start[i],
                                         start          = interval_start[i] - before[i],
                                         end            = interval_start[i] + after[i], 
                                         exclude        = exclude_dates,
                                         control.dates  = control_dates,
                                         knots.per.year = knots[i],
                                         weekday.effect = TRUE,
                                         model          = "correlated",
                                         discontinuity  = disc[i], 
                                         verbose = FALSE))
    
    ndays <- 365*2
    excess_cumulative(f, 
                      start = interval_start[i], 
                      end   = ymd(interval_start[i]) + ndays) %>%
      mutate(agegroup = x, event_day = interval_start[i], event = names(interval_start)[i]) %>%
      as_tibble()
    
  })
  tmp %>% 
    group_by(date) %>% 
    summarize(fitted = sum(fitted),
              observed = sum(observed),
              sd = sqrt(sum(sd^2)),
              se = sqrt(sum(se^2)),
              event_day = event_day[1], 
              event = event[1]) %>%
    ungroup()
})

# -- Supp figure 5
supp_fig5 <- excess_deaths_pr %>%
  mutate(lwr = fitted-1.96*se,
         upr = fitted+1.96*se) %>%
  filter(event != "Hugo") %>%
  mutate(day = as.numeric(date - event_day)) %>%
  mutate(event = factor(event, levels = c("Maria", "Georges", "Hugo", "Chikungunya", "Covid-19", "Flu-2005"))) %>% 
  ggplot(aes(day, fitted)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha=0.50) +
  geom_line(show.legend = F) +
  ylab("Cumulative excess deaths") +
  xlab("Days after the event") +
  facet_wrap(~event) +
  scale_x_continuous(limits = c(0, ndays),
                     breaks = seq(0, ndays, by=100)) +
  scale_y_continuous(limits = c(-400, 4000),
                     breaks = seq(0, 4000, by=1000),
                     labels = scales::comma) +
  theme(axis.title = element_text(size=13),
        axis.text  = element_text(size=13),
        legend.title      = element_blank())

# -- Save supp figure 5
ggsave("figs/supp-figure-5.pdf",
       plot   = supp_fig5,
       dpi    = 300, 
       height = 4,
       width  = 6)
### -- -------------------------------------- ------------------------------------------------------------------
### -- END Supp figure 5: Excess deaths in PR ------------------------------------------------------------------
### -- -------------------------------------- ------------------------------------------------------------------
