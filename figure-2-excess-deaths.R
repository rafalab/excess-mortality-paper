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
library(directlabels)
dslabs::ds_theme_set()

# -- Hurricanes information
hurricane_dates        <- as.Date(c("1989-09-18","1998-09-21","2017-09-20"))
hurricane_effect_ends  <- as.Date(c("1990-03-18","1999-03-21","2018-03-20"))
names(hurricane_dates) <- c("Hugo", "Georges", "Maria")

## define periods not used to compute the expected counts
## remove 6 months after hurricane
## we also remove:
##  - chikingunya fall+winder in 2014,
##  - a couple of weeks of clear outliers in 2001
##  - 2020 and beyond as the data might be incomplete
exclude_dates <- c(seq(hurricane_dates[1], hurricane_effect_ends[1], by = "day"),
                   seq(hurricane_dates[2], hurricane_effect_ends[2], by = "day"),
                   seq(hurricane_dates[3], hurricane_effect_ends[3], by = "day"),
                   seq(as.Date("2014-09-01"), as.Date("2015-03-21"), by = "day"),
                   seq(as.Date("2001-01-01"), as.Date("2001-01-15"), by = "day"),
                   seq(as.Date("2020-01-01"), lubridate::today(), by = "day"))

# -- Control dates
control_dates <- seq(as.Date("2002-01-01"), as.Date("2013-12-31"), by = "day")

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
nknots <- 6
ndays  <- 365
disc   <- c(TRUE, TRUE, TRUE, FALSE, FALSE)
before <- c(365, 365, 365, 365, 548) 
after  <- c(365, 365, 365, 365, 90)
interval_start <- c(hurricane_dates[1], # Hugo
                    hurricane_dates[2], # Georges
                    hurricane_dates[3], # Maria
                    Chikungunya = make_date(2014, 8, 1),
                    "Covid-19"  = make_date(2020, 1, 1))

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
                   # knots.per.year = round(npy*(before[i] + after[i])/365), 
                   knots.per.year = nknots,
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
  dates <- filter(fit, date >= interval_start[i], fitted-1.96*se >= 0)$date
  ind   <- suppressWarnings(min(which(diff(dates) > 1)))
  
  # -- Find last day
  if(is.infinite(ind)) {
    last_day <- ymd(last(dates))
    if(is.na(last_day)) {last_day <- interval_start[i]+1}
  } else {
    last_day <- ymd(dates[ind])
  }

  
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
                   # knots.per.year = round(npy*(before[i] + after[i])/365), 
                   knots.per.year = nknots,
                   model          = "correlated",
                   discontinuity  = disc[i], 
                   verbose = FALSE))
    
    excess_cumulative(f, 
                      start = interval_start[i], 
                      end   = last_day) %>%
      mutate(agegroup = x, event_day = interval_start[i], event = names(interval_start)[i])
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

# -- Figure 2A
fig2a <- excess_deaths_pr %>%
  mutate(day = as.numeric(date - event_day)) %>%
  mutate(event = factor(event, levels = c("Maria", "Georges", "Hugo", "Chikungunya", "Covid-19"))) %>% 
  ggplot(aes(color = event, fill = event)) +
  geom_ribbon(aes(day, ymin = fitted - 2*se, ymax = fitted + 2*se), alpha = 0.25, show.legend = F) + 
  geom_point(aes(day, observed), size=1, alpha = 0.25, show.legend = F) +
  geom_line(aes(day, fitted), size=1) +
  ylab("Cumulative excess deaths") +
  xlab("Days after the event") +
  scale_color_manual(name="",
                     values = c("#D55E00","#0571b0","#009E73","#56B4E9","#CC79A7","#E69F00","#ca0020","gray")) +
  scale_fill_manual(name="",
                     values = c("#D55E00","#0571b0","#009E73","#56B4E9","#CC79A7","#E69F00","#ca0020","gray")) +
  theme(axis.title = element_text(face="bold", color="black"),
        axis.text  = element_text(face="bold", color="black"),
        legend.title    = element_blank(),
        legend.text     = element_text(face="bold", color="black", size=8),
        legend.background = element_rect(color    = "black",
                                         fill     = "white",
                                         linetype = "solid"),
        legend.position = c(0.50, 0.05),
        legend.key.size = unit(0.1, "cm"),
        legend.direction = "horizontal")

# -- Save figure 2A
ggsave("figs/figure-2a.pdf",
       plot   = fig2a,
       dpi    = 300, 
       height = 8.50,
       width  = 11.0)
### -- ---------------------------------- ------------------------------------------------------------------
### -- END Figure 2A: Excess deaths in PR ------------------------------------------------------------------
### -- ---------------------------------- ------------------------------------------------------------------

### -- ------------------------------ ------------------------------------------------------------------
### -- Figure 2B: Excess deaths in US ------------------------------------------------------------------
### -- ------------------------------ ------------------------------------------------------------------
nknots <- 20
# -- Load US state data and and expand abbreviations
data(cdc_state_counts)
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
max_date      <- make_date(2020, 5, 2)

# -- Remove data after the max date
counts <- cdc_state_counts %>% filter(date <= max_date)

# -- Taking some states out due to low sample size
states <- sort(unique(counts$state))
states <- setdiff(states, c("Connecticut", "North Carolina"))

# -- Fitting model to each date. Using 6 knots as before
fits   <- lapply(states, function(x){
  print(x)
  if(x == "Puerto Rico"){
    exclude_dates <- unique(sort(c(exclude_dates, seq(make_date(2017, 9, 20), make_date(2018, 3, 31), by = "day"))))
  }
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

# -- Assigning a name to each element of the fits list
names(fits) <- states

# -- Calculating excess deaths for each state
excess_deaths_us <- map_df(fits, function(f){
  ret <- excess_cumulative(f, 
                           start = make_date(2020, 03, 01), 
                           end   = max_date) 
  ret$state <- f$state
  return(ret)
})

# -- Cumulative deaths by Covid-19 in USA
covid_us <- covid_states %>% 
  filter(date >= min(excess_deaths_us$date)) %>%
  group_by(date) %>%
  summarize(covid = sum(death))

# -- Cumulative deaths in USA
us <- excess_deaths_us %>%
  as_tibble() %>%
  group_by(date) %>%
  summarize(observed = sum(observed),
            sd       = sqrt(sum(sd^2)),
            fitted   = sum(fitted),
            se       = sqrt(sum(se^2))) %>%
  ungroup()

# -- Figure 2B
fig2b <- us %>%
  left_join(covid_us, by="date") %>%
  ggplot() +
  
  geom_line(aes(date, fitted, color="Fitted")) +
  geom_point(aes(date, fitted), size=3, color="white", show.legend = F) +
  geom_point(aes(date, fitted, color="Fitted"), size=3, alpha=0.50, show.legend = F) +
  
  geom_line(aes(date, observed, color="Observed")) +
  geom_point(aes(date, observed), size=3, color="white", show.legend = F) +
  geom_point(aes(date, observed, color="Observed"), size=3, alpha=0.50, show.legend = F) +
  
  geom_line(aes(date, covid, color="Covid-19")) +
  geom_point(aes(date, covid), size=3, color="white", show.legend = F) +
  geom_point(aes(date, covid, color="Covid-19"), size=3, alpha=0.50, show.legend = F) +
  
  ylab("Cumulative excess deaths") +
  xlab("") +
  
  scale_y_continuous(limits = c(-1500, 100000),
                     breaks = seq(0, 100000, by=10000),
                     labels = scales::comma) +
  scale_x_date(date_breaks = "1 week", date_labels = "%b %y") +
  
  scale_color_manual(name="",
                    values = c("#D55E00","#0571b0","#009E73","#56B4E9","#CC79A7","#E69F00","#ca0020","gray")) +
  theme(axis.title = element_text(face="bold", color="black"),
        axis.text  = element_text(face="bold", color="black"),
        legend.title    = element_blank(),
        legend.text     = element_text(face="bold", color="black", size=8),
        legend.background = element_rect(color    = "black",
                                         fill     = "white",
                                         linetype = "solid"),
        legend.position = c(0.50, 0.05),
        legend.key.size = unit(0.1, "cm"),
        legend.direction = "horizontal")

# -- Save figure 2B
ggsave("figs/figure-2b.pdf",
       plot   = fig2b,
       dpi    = 300, 
       height = 8.50,
       width  = 11.0)
### -- ---------------------------------- ------------------------------------------------------------------
### -- END Figure 2B: Excess deaths in US ------------------------------------------------------------------
### -- ---------------------------------- ------------------------------------------------------------------

### -- ----------------------------------------- ------------------------------------------------------------------
### -- Figure 2C: Covid19 vs Flu18 excess deaths ------------------------------------------------------------------
### -- ----------------------------------------- ------------------------------------------------------------------
# -- Intervals for Flu 18 and Covid 19, respectively
intervals <- list(seq(make_date(2017, 12, 16), make_date(2018, 2, 10), by = "day"),
                  seq(make_date(2020, 03, 14), max_date, by = "day"))

# -- Fitting model to each state 
fits <- map_df(states, function(x){
  print(x)
  if(x == "Puerto Rico"){
    exclude_dates <- unique(sort(c(exclude_dates, seq(make_date(2017, 9, 20), make_date(2018, 3, 31), by = "day"))))
  }
  
  fit <- suppressMessages(counts %>% 
    filter(state == x) %>%
    excess_model(exclude        = exclude_dates,
                 intervals      = intervals,
                 weekday.effect = FALSE,
                 knots.per.year = nknots) %>%
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

# -- Figure 2C
fig2c <- fits %>% 
  as_tibble() %>%
  mutate(virus = ifelse(year(start)==2020, "COVID_19", "Flu_18")) %>%
  select(state, virus, excess) %>%
  mutate(excess = pmax(0.5, excess)) %>%
  left_join(pop, by = "state") %>%
  mutate(state = str_remove(state, " \\(not including NYC\\)")) %>%
  mutate(abb = state.abb.2[match(state, state.name.2)]) %>%
  spread(virus, excess) %>%
  ggplot(aes(Flu_18/pop_flu*10000, COVID_19/pop_covid*10000, label = abb)) +
  geom_point(size=2, alpha=0.50) +
  geom_point(size=2, pch=1) +
  ggrepel::geom_text_repel(fontface="bold") +
  scale_y_continuous(limits = c(0, 30),
                     breaks = seq(0, 30, by=5)) +
  scale_x_continuous(limits = c(0, 3),
                     breaks = seq(0, 3, by=0.50)) +
  geom_abline(intercept = 0, slope = 1, color="#cb181d", lty=2) +
  ylab("Covid-19 mortality rate per 10,000") +
  xlab("Seasonal flu (2018) mortality rate per 10,000") +
  theme(axis.title = element_text(face="bold", color="black"),
        axis.text  = element_text(face="bold", color="black"),
        legend.title    = element_blank(),
        legend.text     = element_text(face="bold", color="black", size=8),
        legend.background = element_rect(color    = "black",
                                         fill     = "white",
                                         linetype = "solid"),
        legend.position = c(0.50, 0.05),
        legend.direction = "horizontal")

# -- Save figure 2C
ggsave("figs/figure-2c.pdf",
       plot   = fig2c,
       dpi    = 300, 
       height = 8.50,
       width  = 11.0)
### -- --------------------------------------------- ------------------------------------------------------------------
### -- END Figure 2C: Covid19 vs Flu18 excess deaths ------------------------------------------------------------------
### -- --------------------------------------------- ------------------------------------------------------------------

### -- --------------------------------------------- ------------------------------------------------------------------
### -- Figure 2D: Cook county Covid19 white vs black ------------------------------------------------------------------
### -- --------------------------------------------- ------------------------------------------------------------------
# -- Loading cook county data
data("cook_records")

# -- Wrangling mortality & population
the_breaks <- c(seq(0, 85, by=5), Inf)
counts     <- compute_counts(cook_records, by=c("sex", "race", "agegroup"), breaks = the_breaks)
counts     <- left_join(counts, cook_demographics, by=c("date", "race", "sex", "agegroup"))

# -- Creating new agegroups and wrangle
the_breaks <- c(0, 5, 20, 40, 60, 75, Inf)
counts     <- collapse_counts_by_age(counts, the_breaks) %>%
  filter(sex  %in% c("male", "female"),
         race %in% c("white", "black")) %>%
  mutate(agegroup = factor(agegroup, levels = c("0-4", "5-19", "20-39", "40-59", "60-74", "75-Inf"))) %>%
  mutate(group = paste0(sex," | ",race," | ",agegroup)) %>%
  filter(!agegroup %in% c("0-4", "5-19"))

# -- Fitting models
exclude_dates <- seq( ymd("2020-01-01"), max(counts$date), by="days")
control_dates <- seq(min(counts$date), ymd("2019-12-31"), by="days")
excess_deaths_cook <- map_df(unique(counts$group), function(x){
  
  print(x)
  tmp_counts <- filter(counts, group == x)
  f          <- excess_model(counts         = tmp_counts,
                             exclude        = exclude_dates,
                             start          = ymd("2020-01-01"),
                             end            = max(counts$date),
                             weekday.effect = FALSE,
                             knots.per.year = nknots,
                             verbose        = FALSE)
  
  ret <- excess_cumulative(f, 
                           start = make_date(2020, 03, 01), 
                           end   = max(counts$date)) %>%
    mutate(sex = tmp_counts$sex[1], race = tmp_counts$race[1], agegroup = tmp_counts$agegroup[1])
}) %>% as_tibble()

# -- Figure 2D
fig2d <- excess_deaths_cook %>%
  group_by(date, sex, race) %>%
  summarize(observed = sum(observed),
            sd       = sqrt(sum(sd^2)),
            fitted   = sum(fitted),
            se       = sqrt(sum(se^2))) %>%
  ungroup() %>%
  mutate(group = paste0(race," ",sex,"s"),
         group = str_to_sentence(group)) %>%
  ggplot(aes(x=date, y=fitted, color=group, fill=group)) +
  geom_ribbon(aes(ymin=fitted-1.96*se, ymax=fitted+1.96*se), alpha=0.20, show.legend = F) +
  geom_point(aes(date, observed), size=0.80, alpha=0.40, show.legend = F) +
  geom_line(size=1) +
  ylab("Cumulative excess deaths") +
  xlab("") +
  scale_x_date(date_breaks = "10 days", date_labels = "%b %d") +
  scale_color_manual(name="",
                     values = c("#D55E00","#0571b0","#009E73","#56B4E9","#CC79A7","#E69F00","#ca0020","gray")) +
  scale_fill_manual(name="",
                    values = c("#D55E00","#0571b0","#009E73","#56B4E9","#CC79A7","#E69F00","#ca0020","gray")) +
  theme(axis.title = element_text(face="bold", color="black"),
        axis.text  = element_text(face="bold", color="black"),
        legend.title    = element_blank(),
        legend.text     = element_text(face="bold", color="black", size=8),
        legend.background = element_rect(color    = "black",
                                         fill     = "white",
                                         linetype = "solid"),
        legend.position = c(0.50, 0.05),
        legend.key.size = unit(0.1, "cm"),
        legend.direction = "horizontal")

# -- Save figure 2D
ggsave("figs/figure-2d.pdf",
       plot   = fig2d,
       dpi    = 300, 
       height = 8.50,
       width  = 11.0)
### -- ------------------------------------------------- ------------------------------------------------------------------
### -- END Figure 2D: Cook county Covid19 white vs black ------------------------------------------------------------------
### -- ------------------------------------------------- ------------------------------------------------------------------

# -- Create Figure 2
fig2 <- ggarrange(fig2a,fig2b,fig2c,fig2d,
          nrow=2, ncol=2, 
          labels=LETTERS[1:4])

# -- Save figure 2
ggsave("figs/figure-2.pdf",
       plot   = fig2,
       dpi    = 300, 
       height = 8.50,
       width  = 11.0)
