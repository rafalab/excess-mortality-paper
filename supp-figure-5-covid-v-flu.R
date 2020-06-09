# -- Libraries
library(scales)
library(ggpubr)
library(tidyverse)
library(lubridate)
library(excessmort)
library(directlabels)
dslabs::ds_theme_set()

# -- Loading data
data("cdc_state_counts")
state.name.2 <- c(state.name, "New York City", "Puerto Rico", "District of Columbia")
state.abb.2  <- c(state.abb, "NYC", "PR", "DC")

# -- Import covid-19 reported deaths data and wrangle
covid_nyc    <- read_csv("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv")
covid_nyc    <- filter(covid_nyc, county =="New York City")
covid_states <- read_csv("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv") %>%
  mutate(abb = state.abb.2[match(state, state.name.2)]) %>%
  mutate(date = ymd(date)) %>%
  filter(!is.na(state)) %>%
  arrange(state)

# -- Split New York into NYC and the rest of NY
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

# -- Define regions of interest
flu_season    <- seq(make_date(2017, 12, 16), make_date(2018, 1, 16), by = "day")
exclude_dates <- c(flu_season, seq(make_date(2020, 1, 1), max(cdc_state_counts$date, na.rm = TRUE), by = "day"))
max_date      <- max(cdc_state_counts$date)-7

# -- Remove data after the max date
counts <- cdc_state_counts %>% filter(date <= max_date)

# -- Taking some states out due to low sample size
states <- sort(unique(counts$state))
states <- setdiff(states, c("Connecticut", "North Carolina", "Puerto Rico"))

### -- ------------------------------------- ------------------------------------------------------------------
### -- Supp Figure 5: Covid19 vs Flu18 rates ------------------------------------------------------------------
### -- ------------------------------------- ------------------------------------------------------------------
# -- Intervals for Flu 18 and Covid 19, respectively
intervals <- list(seq(make_date(2017, 12, 10), make_date(2018, 2, 24), by = "day"),
                  seq(make_date(2020, 03, 01), max_date, by = "day"))

# -- Fitting model to each state 
fits <- map_df(states, function(x){
  print(x)
  if(x == "Puerto Rico"){
    exclude_dates <- unique(sort(c(exclude_dates, seq(make_date(2017, 9, 20), make_date(2018, 3, 31), by = "day"))))
  }
  
  fit <- suppressMessages(counts %>% 
                            filter(state == x) %>%
                            excess_model(exclude   = exclude_dates,
                                         intervals = intervals) %>%
                            mutate(state = x))
})

# -- Average population size during Flu 2018
pop1 <- counts %>%
  filter(date %in% intervals[[1]]) %>%
  group_by(state) %>%
  summarize(pop_flu = mean(population, na.rm=TRUE))

# -- Average population size during Covid-19
pop2 <- counts %>%
  filter(date %in% intervals[[2]]) %>%
  group_by(state) %>%
  summarize(pop_covid = mean(population, na.rm=TRUE))

# -- Putting population data together
pop <- left_join(pop1, pop2, by = 'state')

# -- Wrangling before the viz
fits <- fits %>% 
  as_tibble() %>%
  mutate(virus = ifelse(year(start)==2020, "COVID_19", "Flu_18")) %>%
  select(state, virus, excess) %>%
  mutate(excess = pmax(0.5, excess)) %>%
  left_join(pop, by = "state") %>%
  mutate(state = str_remove(state, " \\(not including NYC\\)")) %>%
  mutate(abb = state.abb.2[match(state, state.name.2)]) %>%
  filter(!abb %in% c("CT", "PR","NC")) %>%
  spread(virus, excess) 

# -- Supp fig rates
supp_fig <- fits %>%
  ggplot(aes(Flu_18/pop_flu*100000, COVID_19/pop_covid*100000, label = abb)) +
  geom_abline(intercept = 0, slope = 1, color="#cb181d", lty=2) +
  geom_point(alpha=0.50) +
  geom_point(pch=1) +
  ggrepel::geom_text_repel(force = 0.06) +
                           # data = filter(fits, Flu_18 >= 1000 | COVID_19 >= 5000)) +
  scale_y_continuous(breaks = seq(0, 300, by=50)) +
  scale_x_continuous(breaks = seq(0, 25, by=5)) +
  ylab("Covid-19 excess deaths per 100,000") +
  xlab("Seasonal flu (2018) excess deaths per 100,000") +
  theme(axis.text  = element_text(size=18),
        axis.title = element_text(size=18))

# -- Save figure 2C
ggsave("figs/supp-figure-5-covid-v-flu.pdf",
       plot   = supp_fig,
       dpi    = 300, 
       height = 5,
       width  = 7)
### -- --------------------------------------------- ------------------------------------------------------------------
### -- END Figure 2C: Covid19 vs Flu18 excess deaths ------------------------------------------------------------------
### -- --------------------------------------------- ------------------------------------------------------------------
