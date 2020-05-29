library(tidyverse)
library(lubridate)
library(excessmort)
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
max_date <- max(cdc_state_counts$date) - 7 
## remover the latest
counts <- cdc_state_counts %>% filter(date <= max_date)


## take out states with incomplete data
states <- counts %>% 
  filter(!is.na(outcome)) %>% 
  group_by(state) %>% summarize(max = max(date)) %>% 
  filter(max == max(counts$date)) %>%
  pull(state) 

## Only states, and also weird data
states <- setdiff(states, "Puerto Rico")

setdiff(counts$state, states)

### check
check <- FALSE
if(check){
  sapply(unique(counts$state), function(x){
    print(x)
    if(x == "Puerto Rico"){
      exclude_dates <- unique(sort(c(exclude_dates, seq(make_date(2017, 9, 20), make_date(2018, 3, 31), by = "day"))))
    }
    ret <- counts %>% filter(state == x) %>%
      compute_expected(exclude = exclude_dates,
                       verbose = FALSE)
    print(attr(ret, "frequency"))
    print(expected_plot(ret, title = x))
    
    print(qplot(date, 100*(population-min(population))/min(population), data = filter(counts, state == x), main = x))
    NULL
  })
}

kpy <- 24
fits <- lapply(states, function(x){
  if(x == "Puerto Rico"){
    exclude_dates <- unique(sort(c(exclude_dates, seq(make_date(2017, 9, 20), make_date(2018, 3, 31), by = "day"))))
  }
  ret <- counts %>% filter(state == x) %>%
    excess_model(exclude = exclude_dates,
                 start = min(counts$date),
                 end = max_date,
                 knots.per.year = kpy,
                 verbose = FALSE)
  ret$state <- x
  return(ret)
})
names(fits) <- states
  
# Supplemental Figure - Show f for worse 12 -------------------------------
df <- map_df(fits, function(f)
  with(f, tibble(state = state, 
                 date = date, 
                 expected =  expected, 
                 fitted = 100* fitted, 
                 se = 100* se,
                 sd = 100* sd)))

show <- df %>% group_by(state) %>%
#  filter(!state %in% "Puerto Rico") %>% # remove due to María effect.
  summarize(max = max(fitted)) %>%
  ungroup() %>%
  arrange(desc(max)) %>%
  top_n(12, max) %>%
  pull(state)

df %>% 
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
df %>% 
  group_by(date) %>% 
  summarize(fitted = sum(expected * fitted) / sum(expected), 
            se = sqrt(sum(expected^2 * se^2)) / sum(expected),
            sd = sqrt(sum(expected^2 * sd^2)) / sum(expected)) %>%
  ggplot(aes(date, fitted, ymin = fitted - 2.54*se, ymax = fitted + 2.54*se)) +
  geom_ribbon(aes(date, ymin = -2*sd, ymax = 2*sd), color = 1, fill = NA, lty = 2) +
  geom_hline(yintercept = 0) +
  geom_ribbon(alpha = 0.5) +
  geom_line()


# Excess mortality --------------------------------------------------------

e <- map_df(fits, function(f){
  ret <- excess_cumulative(f, 
                    start = make_date(2020, 3, 14), 
                    end = max_date) 
  ret$state <- f$state
  return(ret)
})

tmp <- covid_states %>% 
  filter(date >= min(e$date)) %>%
  group_by(date) %>%
  summarize(covid = sum(death))

us <- e %>% group_by(date) %>%
  summarize(observed = sum(observed),
            sd = sqrt(sum(sd^2)),
            fitted = sum(fitted),
            se = sqrt(sum(se^2))) %>%
  ungroup()
  
us %>% 
  ggplot(aes(date)) +
  geom_ribbon(aes(ymin = observed- 2*sd, ymax = observed + 2*sd), alpha = 0.5, fill = "blue") +
  geom_point(aes(y = observed), col = "blue") +
  geom_line(aes(x=date, y=covid), data =tmp) +
  geom_hline(yintercept = 100000, lty = 2) +
  scale_y_continuous(breaks = seq(0,120000,25000), labels = scales::comma) +
  annotate("text",  x = make_date(2020, 05, 10), y = 110000, label = "Excess mortality\nestimate", hjust = 0, color = "blue") +
  annotate("text",  x = make_date(2020, 05, 10), y = 70000, label = "Total reported\nCOVID-19 Deaths", hjust = 0) +
  ylab("Total Deaths") + xlab("Date") +
  ggtitle("Excess Mortality Estimate for USA Based on CDC data") +
  labs(caption = "Excess mortality estimate excludes Connecticut and North Carolina due to missing data.\nIncludes all other 48 states.")+
  theme(plot.caption = element_text(hjust = 0))

ggsave("~/Desktop/excess-mort.jpg", width = 6, height = 6/1.6)


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
                  knots.per.year =npy) %>%
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
