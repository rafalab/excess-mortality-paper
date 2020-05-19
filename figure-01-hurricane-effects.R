### -- Set up --------------------------------------------------------
# -- Libraries
library(tidyverse)
library(lubridate)
library(excessdeaths)

# source("pr-init.R")
dslabs::ds_theme_set()

data("puerto_rico_counts")
data("florida_counts")
data("new_jersey_counts")
data("louisiana_counts")

## remove the outlier from louisana
louisiana_counts$outcome[which.max(louisiana_counts$outcome)] <- 126

the_breaks <- c(0, 5, 20, 40, 60, 75, Inf)
counts <- list(florida_counts, louisiana_counts, new_jersey_counts, 
               collapse_counts_by_age(puerto_rico_counts, the_breaks))
               


# -- Hurricane dates
hurricane_dates <- list(irma = as.Date(c("2017-09-10")),
                        katrina = as.Date(c("2005-08-29")),
                        sandy = as.Date(c("2012-10-29")),
                        hugo = as.Date(c("1989-09-18")),
                        georges = as.Date(c("1998-09-21")),
                        maria = as.Date("2017-09-20"))
# -- Hurricane names
hurricane_names <-  list(irma = "FL: Irma",
                         katrina = "LA: Katrina",
                         sandy  = "NJ: Sandy",
                         hugo = "PR: Hugo", 
                         georges = "PR: Georges",
                         maria = "PR: Maria")

count_index <- c(1, 2, 3, 4, 4, 4)

control_dates <- list(irma = seq(make_date(2015,1,1), make_date(2016, 12, 31), by = "day"),
                      katrina = seq(make_date(2003,1,1), make_date(2005, 8, 1), by = "day"),
                      sandy = seq(make_date(2007,1,1), make_date(2012, 10, 1), by = "day"),
                      hugo = seq(make_date(2002, 1, 1), make_date(2013, 12, 31), by = "day"),
                      georges = seq(make_date(2002, 1, 1), make_date(2013, 12, 31), by = "day"),
                      maria = seq(make_date(2002, 1, 1), make_date(2013, 12, 31), by = "day"))

puerto_rico_hurricane_dates <- as.Date(c("1989-09-18","1998-09-21","2017-09-20"))
puerto_rico_hurricane_effect_ends <- as.Date(c("1990-03-18","1999-03-21","2018-03-20"))

puerto_rico_out_dates <- c(
  seq(puerto_rico_hurricane_dates[1], puerto_rico_hurricane_effect_ends[1], by = "day"),
  seq(puerto_rico_hurricane_dates[2], puerto_rico_hurricane_effect_ends[2], by = "day"),
  seq(puerto_rico_hurricane_dates[3], puerto_rico_hurricane_effect_ends[3], by = "day"),
  seq(as.Date("2014-09-01"), as.Date("2015-03-21"), by = "day"),
  seq(as.Date("2001-01-01"), as.Date("2001-01-15"), by = "day"),
  seq(as.Date("2020-01-01"), lubridate::today(), by = "day"))

exclude_dates  <- list(irma = hurricane_dates[["irma"]] + 0:180,
                       katrina = hurricane_dates[["katrina"]]  + 0:180,
                       sandy = hurricane_dates[["sandy"]]  + 0:180,
                       hugo = puerto_rico_out_dates,
                       georges = puerto_rico_out_dates,
                       maria = puerto_rico_out_dates)

before <-days(122)
after <- days(244)
nknots <- 6


tmp <- map_df(seq_along(count_index), function(i){
  h <- count_index[i]
  if(i < 4){ #not PR
    fit <-  excess_model(counts[[h]], 
                         event = hurricane_dates[[i]],
                         start =  hurricane_dates[[i]] - before,
                         end =  hurricane_dates[[i]] + after, 
                         exclude = exclude_dates[[i]],
                         control.dates = control_dates[[i]],
                         nknots = nknots,
                         model = "correlated")
    ret <- tibble(date = fit$date, fitted = fit$fitted, se = fit$se, 
                  hurricane = hurricane_names[[i]],
                  hurricane_date = hurricane_dates[[i]])
  } else{
    tmp <- map_df(unique(counts[[h]]$agegroup), function(x){
      fit <- counts[[h]] %>% filter(agegroup == x) %>%
        excess_model(event = hurricane_dates[[i]],
                     start =  hurricane_dates[[i]] - before,
                     end =  hurricane_dates[[i]] + after, 
                     exclude = exclude_dates[[i]],
                     control.dates = control_dates[[i]],
                     nknots = nknots, 
                     model = "correlated")
      tibble(date = fit$date, expected = fit$expected, fitted = fit$fit, se = fit$se)
    })
    ret <- tmp %>% group_by(date) %>% 
      summarize(fitted = sum(expected * fitted) / sum(expected), 
                se = sqrt(sum(expected^2 * se^2)) / sum(expected)) %>%
      mutate(hurricane = hurricane_names[[i]],
             hurricane_date = hurricane_dates[[i]])
  }
  return(ret)
}) %>%
  mutate(hurricane = reorder(hurricane, date, min))


##Figure 1
tmp %>% mutate(day = as.numeric(date - hurricane_date)) %>%
  ggplot(aes(day, fitted*100, color = hurricane)) +
  xlab("Days since the hurricane") +
  ylab("Percent death rate increase") +
  geom_line() 

## Supp Figure
tmp %>% 
  ggplot(aes(date, fitted*100)) +
  geom_ribbon(aes(ymin = fitted*100 - 2*se*100, ymax = fitted*100 + 2*se*100), alpha = 0.5, color = NA) +
  geom_line() +
  facet_wrap( ~ hurricane, scale= "free_x")
