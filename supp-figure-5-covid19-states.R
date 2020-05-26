### -- ------ ------------------------------------------------------------------
### -- Set up ------------------------------------------------------------------
### -- ------ ------------------------------------------------------------------
# -- Libraries
library(ggpubr)
library(tidyverse)
library(lubridate)
library(excessmort)
dslabs::ds_theme_set()

# -- Loading state mortality data 
data(cdc_state_counts)

# -- Expand state abbrevaition objects 
state.name.2 <- c(state.name, "New York City", "Puerto Rico", "District of Columbia")
state.abb.2  <- c(state.abb, "NYC", "PR", "DC")

# -- Importing covd-19 reported deaths data 
covid_nyc    <- read_csv("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv")
covid_nyc    <- filter(covid_nyc, county =="New York City")
covid_states <- read_csv("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv") %>%
  mutate(abb = state.abb.2[match(state, state.name.2)]) %>%
  mutate(date = ymd(date)) %>%
  filter(!is.na(state)) %>%
  arrange(state)

# -- Subsetting new york data
ny <- filter(covid_states, state == "New York")

# -- Covid 19 data for the rest of new york
ny <- left_join(ny, covid_nyc, by = "date") %>%
  mutate(death = deaths.x - deaths.y, 
         state = "Rest of New York") %>%
  select(date, state, death)

# -- Covid 19 data for states
covid_states <- filter(covid_states, state!="New York") %>%
  rename(death = deaths) %>%
  select(date, state, death)

# -- Covid 19 for NYC
covid_nyc <- covid_nyc %>%
  mutate(state = "New York City", death = deaths)%>%
  select(date, state, death)

# -- Putting data together
covid_states <- bind_rows(covid_states, ny, covid_nyc) %>%
  filter(!is.na(death))
rm(covid_nyc, ny)

# -- Denoting periods of interest
flu_season    <- seq(make_date(2017, 12, 16), make_date(2018, 1, 16), by = "day")
exclude_dates <- c(flu_season, seq(make_date(2020, 1, 1), max(cdc_state_counts$date, na.rm = TRUE), by = "day"))
max_date      <- make_date(2020, 5, 2)

# -- Remove last dates
counts <- cdc_state_counts %>% filter(date <= max_date)
states <- unique(counts$state)

# -- Number of knots per year
nknots <- 6
### -- ------ ------------------------------------------------------------------
### -- END Set up ------------------------------------------------------------------
### -- ------ ------------------------------------------------------------------

### -- ------------------------------------------ ------------------------------------------------------------------
### -- Supp Figure 5: Fhat for worse 12 US states ------------------------------------------------------------------
### -- ------------------------------------------ ------------------------------------------------------------------
# -- Fitting the model to each state
fits <- lapply(states, function(x){
  if(x == "Puerto Rico"){
    exclude_dates <- unique(sort(c(exclude_dates, seq(make_date(2017, 9, 20), make_date(2018, 3, 31), by = "day"))))
  }
  print(x)
  ret <- counts %>% 
    filter(state == x) %>%
    excess_model(exclude        = exclude_dates,
                 start          = min(counts$date),
                 end            = max_date,
                 weekday.effect = FALSE,
                 knots.per.year = nknots,
                 verbose        = FALSE)
  ret$state <- x
  return(ret)
})
names(fits) <- states

# -- Extracting components from model fit
df <- map_df(fits, function(f)
  with(f, tibble(state = state, 
                 date = date, 
                 expected =  expected, 
                 fitted = 100* fitted, 
                 se = 100* se,
                 sd = 100* sd)))

# -- Extract worse 12 states (including NYC)
show <- df %>% 
        group_by(state) %>%
        filter(!state %in% "Puerto Rico") %>% # remove due to María effect.
        summarize(max = max(fitted)) %>%
        ungroup() %>%
        arrange(desc(max)) %>%
        top_n(12, max) %>%
        pull(state)

# -- Supp figure 5
supp_fig5 <- df %>% 
  filter(state %in% show) %>%
  mutate(state = factor(state, levels = show)) %>%
  ggplot(aes(date, fitted, ymin = fitted - 2*se, ymax = fitted + 2*se)) +
  geom_line() +
  geom_hline(yintercept = 0) +
  geom_ribbon(alpha = 0.5) +
  geom_ribbon(aes(date, ymin = -2*sd, ymax = 2*sd), color = 1, fill = NA, lty = 2) +
  ylab("Percent change in mortality") +
  xlab("") +
  facet_wrap(~state, scale = "free_y")

# -- Saving figure c
ggsave("figs/supp-figure-5.pdf",
       plot   = supp_fig5,
       dpi    = 300, 
       height = 4,
       width  = 8)
### -- ---------------------------------------------- ------------------------------------------------------------------
### -- END Supp Figure 5: Fhat for worse 12 US states ------------------------------------------------------------------
### -- ---------------------------------------------- ------------------------------------------------------------------