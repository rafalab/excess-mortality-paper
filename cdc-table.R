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
  mutate(death = deaths.x - deaths.y, state = "New York (not including NYC)") %>%
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


# Excess mortality --------------------------------------------------------
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
                 knots.per.year = npy) %>%
    mutate(state = x)
})

o <- e %>% select(state, end, excess) %>%
  spread(end, excess) %>%
  setNames(c("state", "flu18_excess", "covid19_excess"))
s <- e  %>% select(state, end, sd) %>%
  spread(end, sd) %>%
  setNames(c("state", "covid19_sd", "flu18_sd")) 
r <- covid_states %>% filter(date == max_date) %>% select(-date) %>% setNames(c("state", "reported_covid_19"))
p <- counts %>% filter(date == max_date) %>% select(state, population)

tab <- o %>% left_join(s, by = "state") %>% left_join(r, by = "state") %>% left_join(p, by = "state") %>%
  mutate( covid19_excess_rate =  covid19_excess/ population*10^6, flu18_excess_rate = flu18_excess/population*10^6) %>%
  select(state, covid19_excess_rate, covid19_excess,  covid19_sd,
         reported_covid_19,  
         flu18_excess_rate,
         flu18_excess,
         flu18_sd) %>%
  arrange(desc( covid19_excess_rate)) %>%
  mutate_if(is.numeric, round)
