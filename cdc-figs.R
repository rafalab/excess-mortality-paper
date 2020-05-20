library(tidyverse)
library(lubridate)
library(excessdeaths)
library(directlabels)
dslabs::ds_theme_set()
data(cdc_state_counts)

flu_season <- seq(make_date(2017, 12, 16), make_date(2018, 1, 16), by = "day")
exclude_dates <- c(flu_season, seq(make_date(2020, 1, 1), max(cdc_state_counts$date, na.rm = TRUE), by = "day"))
max_date <- make_date(2020, 4, 25)

## remover the latest
counts <- cdc_state_counts %>% filter(date <= max_date)
  
states <- setdiff(cdc_state_counts$state, c("Pennsylvania", "North Carolina", "Connecticut"))
fits <- map_df(states, function(x){
  print(x)
  fit <- counts %>% filter(state == x) %>%
    excess_model(exclude = exclude_dates,
                 start =  min(counts$date),
                 end = max_date,
                 weekday.effect = FALSE,
                 nknots = 36)
  tibble(state = x, 
         date = fit$date, 
         expected =  fit$expected, 
         fitted = fit$fitted, 
         se = fit$se,
         sd = sqrt(diag(fit$cov)))
})


## US plot, compare to flue
fits %>% group_by(date) %>% 
  summarize(fitted = sum(expected * fitted) / sum(expected), 
            se = sqrt(sum(expected^2 * se^2)) / sum(expected),
            sd = sqrt(sum(expected^2 * sd^2)) / sum(expected)) %>%
  ggplot(aes(date, fitted, ymin = fitted - 2.54*se, ymax = fitted + 2.54*se)) +
  geom_ribbon(alpha = 0.5) +
  geom_ribbon(aes(date, ymin = -2*sd, ymax = 2*sd), color = 1, fill = NA, lty = 2) +
  geom_line()

fits %>% group_by(date) %>%
  ggplot(aes(date, fitted, color = state))+
  geom_line()
### table

intervals <- list(seq(make_date(2017, 12, 16), make_date(2018, 2, 10), by = "day"),
                  seq(make_date(2020, 3, 14), max_date, by = "day"))
e <- map_df(states, function(x){
  print(x)
  fit <- counts %>% filter(state == x) %>%
    excess_model(exclude = exclude_dates,
                 intervals = intervals,
                 weekday.effect = FALSE,
                 nknots = 36) %>%
    mutate(state = x)
})

tate.name.2
e %>% mutate(virus = ifelse(year(start)==2020, "COVID_19", "Flu_18")) %>%
  select(state, virus, excess) %>%
  mutate(state = str_remove(state, "\\(not including NYC\\)")) %>%
  mutate(abb = state.abb[match(state, state.name)]) %>%
  mutate(abb = state.abb[match(state, state.name)]) %>%
  
  spread(virus, excess) %>%
  ggplot(aes(Flu_18, COVID_19, label = abb)) +
  geom_point() +
  ggrepel::geom_text_repel() +
  geom_abline(intercept = 0, slope = 1) #+
  scale_x_continuous(limits = c(1,25000), trans = "log10")+
  scale_y_continuous(limits = c(1,25000), trans = "log10")
  
  



