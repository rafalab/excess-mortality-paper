library(tidyverse)
library(lubridate)
library(excessdeaths)
library(directlabels)
dslabs::ds_theme_set()
data(cdc_state_counts)


# Expand state abbrevaition objects ---------------------------------------
state.name.2 <- c(state.name, "New York City", "Puerto Rico", "District of Columbia")
state.abb.2 <- c(state.abb, "NYC", "PR", "DC")

# Import covd-19 reported deaths data -------------------------------------
covid_nyc <- read_csv("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv")
covid_nyc <- filter(covid_nyc, county =="New York City")
covid_states <- read_csv("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv") %>%
  mutate(abb = state.abb.2[match(state, state.name.2)]) %>%
  mutate(date = ymd(date)) %>%
  filter(!is.na(state)) %>%
  arrange(state)
# split new york into nyc and the rest of ny
ny <- filter(covid_states, state == "New York")
ny <- left_join(ny, covid_nyc, by = "date") %>%
  mutate(death = deaths.x - deaths.y, state = "Rest of New York") %>%
  select(date, state, death)
covid_states <- filter(covid_states, state!="New York") %>%
  rename(death = deaths) %>%
  select(date, state, death)
covid_nyc <- covid_nyc %>%
  mutate(state = "New York City", death = deaths)%>%
  select(date, state, death)
covid_states <- bind_rows(covid_states, ny, covid_nyc) %>%
  filter(!is.na(death))

rm(covid_nyc, ny)

# Fit excess deaths models ------------------------------------------------
flu_season <- seq(make_date(2017, 12, 16), make_date(2018, 1, 16), by = "day")
exclude_dates <- c(flu_season, seq(make_date(2020, 1, 1), max(cdc_state_counts$date, na.rm = TRUE), by = "day"))
max_date <- make_date(2020, 4, 25)

## remover the latest
counts <- cdc_state_counts %>% filter(date <= max_date)
  
fits <- map_df(states, function(x){
  if(x == "Puerto Rico"){
    exclude_dates <- unique(sort(c(exclude_dates, seq(make_date(2017, 9, 20), make_date(2018, 3, 31), by = "day"))))
  }
  fit <- counts %>% filter(state == x) %>%
    excess_model(exclude = exclude_dates,
                 start =  min(counts$date),
                 end = max_date,
                 weekday.effect = FALSE,
                 nknots = 20,
                 verbose = FALSE)
  tibble(state = x, 
         date = fit$date, 
         expected =  fit$expected, 
         fitted = 100*fit$fitted, 
         se = 100*fit$se,
         sd = 100*sqrt(diag(fit$cov)))
})


# Supplemental Figure - Show f for worse 12 -------------------------------
show <- fits %>% group_by(state) %>%
  fitler(state %in% "Puerto Rico") %>% # remove due to María effect.
  summarize(max = max(fitted)) %>%
  ungroup() %>%
  arrange(desc(max)) %>%
  top_n(12, max) %>%
  pull(state)

fits %>% 
  filter(state %in% show) %>%
  mutate(state = factor(state, levels = show)) %>%
  ggplot(aes(date, fitted, ymin = fitted - 2*se, ymax = fitted + 2*se)) +
  geom_line() +
  geom_hline(yintercept = 0) +
  geom_ribbon(alpha = 0.5) +
  geom_ribbon(aes(date, ymin = -2*sd, ymax = 2*sd), color = 1, fill = NA, lty = 2) +
  facet_wrap(~state, ncol = 3, scale = "free_y") +
  ylab("Percent change") +
  xlab("Date")


# Figure 1C - US effect ---------------------------------------------------
fits %>% 
  #filter(!states %in% c("Pennsylvania", "North Carolina", "Connecticut")) %>%
  group_by(date) %>% 
  summarize(fitted = sum(expected * fitted) / sum(expected), 
            se = sqrt(sum(expected^2 * se^2)) / sum(expected),
            sd = sqrt(sum(expected^2 * sd^2)) / sum(expected)) %>%
  ggplot(aes(date, fitted, ymin = fitted - 2.54*se, ymax = fitted + 2.54*se)) +
  geom_ribbon(aes(date, ymin = -2*sd, ymax = 2*sd), color = 1, fill = NA, lty = 2) +
  geom_hline(yintercept = 0) +
  geom_ribbon(alpha = 0.5) +
  geom_line()


intervals <- list(seq(make_date(2017, 12, 16), make_date(2018, 2, 10), by = "day"),
                  seq(make_date(2020, 3, 14), max_date, by = "day"))
e <- map_df(states, function(x){
  print(x)
  if(x == "Puerto Rico"){
    exclude_dates <- unique(sort(c(exclude_dates, seq(make_date(2017, 9, 20), make_date(2018, 3, 31), by = "day"))))
  }
  fit <- counts %>% filter(state == x) %>%
    excess_model(exclude = exclude_dates,
                 intervals = intervals,
                 weekday.effect = FALSE,
                 nknots = 36) %>%
    mutate(state = x)
})


pop1 <- counts %>%
  filter(date %in% intervals[[1]]) %>%
  group_by(state) %>%
  summarize(pop_flu = mean(population, na.rm=TRUE))
pop2 <- counts %>%
  filter(date %in% intervals[[2]]) %>%
  group_by(state) %>%
  summarize(pop_covid = mean(population, na.rm=TRUE))
pop <- left_join(pop1, pop2, by = 'state')


e %>% mutate(virus = ifelse(year(start)==2020, "COVID_19", "Flu_18")) %>%
  select(state, virus, excess) %>%
  mutate(excess = pmax(0.5, excess)) %>%
  left_join(pop, by = "state") %>%
  mutate(state = str_remove(state, " \\(not including NYC\\)")) %>%
  mutate(abb = state.abb.2[match(state, state.name.2)]) %>%
  spread(virus, excess) %>%
  ggplot(aes(Flu_18/pop_flu*10^6, COVID_19/pop_covid*10^6, label = abb)) +
  geom_point() +
  ggrepel::geom_text_repel() +
  geom_abline(intercept = 0, slope = 1) 
