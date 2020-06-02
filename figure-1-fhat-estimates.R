### -- ------ ------------------------------------------------------------------
### -- Set up ------------------------------------------------------------------
### -- ------ ------------------------------------------------------------------
source("pr-init.R")

# -- Loading data
data("new_jersey_counts")
data("louisiana_counts")
data("florida_counts")

# -- Remove the outlier from louisana
louisiana_counts$outcome[which.max(louisiana_counts$outcome)] <- 126

# -- Age groups
the_breaks <- c(0, 5, 20, 40, 60, 75, Inf)

# -- List of data
counts <- list(florida_counts, louisiana_counts, new_jersey_counts, 
               collapse_counts_by_age(puerto_rico_counts, the_breaks))

# -- Hurricane dates
hurricane_dates <- list(irma    = as.Date(c("2017-09-10")),
                        katrina = as.Date(c("2005-08-29")),
                        sandy   = as.Date(c("2012-10-29")),
                        hugo    = as.Date(c("1989-09-18")),
                        georges = as.Date(c("1998-09-21")),
                        maria   = as.Date("2017-09-20"))

# -- Hurricane names
hurricane_names <-  list(irma    = "FL: Irma",
                         katrina = "LA: Katrina",
                         sandy   = "NJ: Sandy",
                         hugo    = "PR: Hugo", 
                         georges = "PR: Georges",
                         maria   = "PR: Maria")

# -- To be used below
count_index <- c(1, 2, 3, 4, 4, 4)

# -- Control dates for each hurricane
control_dates <- list(irma    = seq(make_date(2015, 01, 01), make_date(2016, 12, 31), by = "day"),
                      katrina = seq(make_date(2003, 01, 01), make_date(2005, 08, 01), by = "day"),
                      sandy   = seq(make_date(2007, 01, 01), make_date(2012, 10, 01), by = "day"),
                      hugo    = seq(make_date(2002, 01, 01), make_date(2013, 12, 31), by = "day"),
                      georges = seq(make_date(2002, 01, 01), make_date(2013, 12, 31), by = "day"),
                      maria   = seq(make_date(2002, 01, 01), make_date(2013, 12, 31), by = "day"))

# -- Period of effect for hurricanes in Puerto Rico
puerto_rico_hurricane_dates       <- as.Date(c("1989-09-18","1998-09-21","2017-09-20"))
puerto_rico_hurricane_effect_ends <- as.Date(c("1990-03-18","1999-03-21","2018-03-20"))

# -- Dates to exclude for hurricanes in Puerto Rico
puerto_rico_out_dates <- exclude_dates

# -- Dates to exclude for each hurricane
exclude_dates  <- list(irma    = hurricane_dates[["irma"]] + 0:180,
                       katrina = hurricane_dates[["katrina"]]  + 0:180,
                       sandy   = hurricane_dates[["sandy"]]  + 0:180,
                       hugo    = puerto_rico_out_dates,
                       georges = puerto_rico_out_dates,
                       maria   = puerto_rico_out_dates)

# -- To be used as parameters in model fitting below
before <- days(365)
after  <- days(365)
### -- ---------- ------------------------------------------------------------------
### -- END Set up ------------------------------------------------------------------
### -- ---------- ------------------------------------------------------------------

### -- ---------------------------------------- ------------------------------------------------------------------
### -- Figure 1A: Fhat estimates for hurricanes ------------------------------------------------------------------
### -- ---------------------------------------- ------------------------------------------------------------------
# -- Number of knots per year
nknots <- 6

# -- Loop to fit models to hurricanes
fits <- map_df(seq_along(count_index), function(i){
  
  print(hurricane_names[[i]])
  # -- Current index
  h <- count_index[i]
  
  if(i < 4){ # -- Here we fit model to hurricanes outside of PR
    
    # -- Fitting model
    fit <- suppressMessages(excess_model(counts[[h]], 
                         event          = hurricane_dates[[i]],
                         start          = hurricane_dates[[i]] - before,
                         end            = hurricane_dates[[i]] + after, 
                         exclude        = exclude_dates[[i]],
                         control.dates  = control_dates[[i]],
                         knots.per.year = nknots,
                         weekday.effect = TRUE,
                         verbose        = FALSE,
                         model          = "correlated"))
    
    # -- Dataframe with fit results
    ret <- tibble(date           = fit$date, 
                  observed       = fit$observed,
                  expected       = fit$expected,
                  fitted         = fit$fitted, 
                  se             = fit$se, 
                  hurricane      = hurricane_names[[i]],
                  hurricane_date = hurricane_dates[[i]])
    
  } else { # -- Here we fit model to hurricanes in PR
    
    # -- This loops through the demographic groups for PR
    tmp <- map_df(unique(counts[[h]]$agegroup), function(x){
      
      # -- Fitting model to each demographic group
      fit <- suppressMessages(counts[[h]] %>% 
        filter(agegroup == x) %>%
        excess_model(event          = hurricane_dates[[i]],
                     start          = hurricane_dates[[i]] - before,
                     end            = hurricane_dates[[i]] + after, 
                     exclude        = exclude_dates[[i]],
                     weekday.effect = TRUE,
                     control.dates  = control_dates[[i]],
                     knots.per.year = nknots, 
                     verbose        = FALSE,
                     model          = "correlated"))
      
      # -- Dataframe with fit results for each demographics
      tibble(date     = fit$date, 
             observed = fit$observed,
             expected = fit$expected, 
             fitted   = fit$fit, 
             se       = fit$se)
    })
    
    # -- Putting PR results together
    ret <- tmp %>% 
      group_by(date) %>% 
      summarize(fitted   = sum(expected * fitted) / sum(expected), 
                se       = sqrt(sum(expected^2 * se^2)) / sum(expected),
                observed  = sum(observed),
                expected = sum(expected)) %>%
      mutate(hurricane      = hurricane_names[[i]],
             hurricane_date = hurricane_dates[[i]])
  }
  return(ret)
}) %>%
  mutate(hurricane = reorder(hurricane, date, min))

# -- Figure 1A
fig1a <- fits %>% 
  mutate(day       = as.numeric(date - hurricane_date),
         hurricane = factor(hurricane, levels=rev(unlist(hurricane_names)))) %>%
  ggplot(aes(day, fitted*100, color = hurricane)) +
  geom_hline(yintercept = 0, lty=2, color="gray") +
  geom_line() +
  xlab("Days since the event") +
  ylab("Percent increase from expected mortality") +
  scale_x_continuous(limits = c(-120, 360),
                     breaks = seq(-100, 300, by=100)) +
  scale_y_continuous(limits = c(-6, 75),
                     breaks = seq(0, 75, by=10)) +
  theme(axis.title = element_text(size=17),
        axis.text  = element_text(size=18),
        legend.title      = element_blank(),
        legend.text       = element_text(size=13),
        legend.background = element_rect(color="black"),
        legend.position   = c(0.80, 0.60))
fig1a

# -- Save figure 1A
ggsave("figs/figure-1a.pdf",
       plot   = fig1a,
       dpi    = 300, 
       height = 5,
       width  = 7)
### -- -------------------------------------------- ------------------------------------------------------------------
### -- END Figure 1A: Fhat estimates for hurricanes ------------------------------------------------------------------
### -- -------------------------------------------- ------------------------------------------------------------------

### -- ----------------------------------------- ------------------------------------------------------------------
### -- Figure 1B: Fhat estimates for Chikungunya ------------------------------------------------------------------
### -- ----------------------------------------- ------------------------------------------------------------------
# -- Number of knots per year
nknots <- 6

# -- Creating breaks for age groups and collapsing data
the_breaks <- c(0, 5, 20, 40, 60, Inf)
all_counts <- collapse_counts_by_age(puerto_rico_counts, the_breaks)

# -- Event dates
event_dates <- ymd("2014-07-14")
events      <- "Chikungunya"

# -- Fitting model only to the groups of interest
res <- map_df(levels(all_counts$agegroup), function(x){
  
  print(x)
  
  tmp_counts <- filter(all_counts, agegroup == x)
  tmp_fit    <- suppressMessages(excess_model(counts          = tmp_counts, 
                                               start          = event_dates - before,
                                               end            = event_dates + after,
                                               exclude        = exclude_dates$maria,
                                               control.dates  = control_dates$maria,
                                               weekday.effect = TRUE,
                                               verbose        = FALSE,
                                               knots.per.year = nknots,
                                               model          = "correlated",
                                               discontinuity  = FALSE))
  
  
  tibble(date       = tmp_fit$date, 
         fitted     = tmp_fit$fitted, 
         observed   = tmp_fit$observed, 
         expected   = tmp_fit$expected, 
         se         = tmp_fit$se, 
         event      = events, 
         event_date = event_dates,
         agegroup   = x)
})

# -- Figure 1B
fig1b <- res %>%
  group_by(date) %>%
  summarize(fitted     = sum(fitted * expected) / sum(expected),
            se         = sqrt(sum(se^2 * expected^2)) / sum(expected),
            event_date = event_date[1]) %>%
  mutate(day = as.numeric(date - event_date),
         lwr = fitted-1.96*se,
         upr = fitted+1.96*se) %>%
  filter(day >= -120, day <= 360) %>%
  ggplot(aes(date, 100*fitted)) +
  geom_hline(yintercept = 0, lty=2, color="gray") +
  geom_ribbon(aes(ymin=100*lwr, ymax=100*upr), alpha=0.50) +
  geom_line(size=1, show.legend = F) +
  xlab("Days since the event") +
  ylab("Percent increase from expected mortality") +
  scale_x_date(date_breaks = "4 months", date_labels = "%b %Y") +
  scale_y_continuous(limits = c(-13, 30),
                     breaks = seq(-10, 30, by=10)) +
  theme(axis.title = element_text(size=17),
        axis.text  = element_text(size=18))

fig1b

# -- Save figure 1B
ggsave("figs/figure-1b.pdf",
       plot   = fig1b,
       dpi    = 300, 
       height = 5,
       width  = 7)
### -- --------------------------------------------- ------------------------------------------------------------------
### -- END Figure 1B: Fhat estimates for Chikungunya ------------------------------------------------------------------
### -- --------------------------------------------- ------------------------------------------------------------------

### -- -------------------------------- ------------------------------------------------------------------
### -- Figure 1C: Fhat estimate for USA ------------------------------------------------------------------
### -- -------------------------------- ------------------------------------------------------------------
# -- Loading state mortality data 
data("cdc_state_counts")

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
max_date      <- make_date(2020, 5, 9)

# -- Remove last dates
counts <- cdc_state_counts %>% filter(date <= max_date)
states <- unique(counts$state)

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

# -- Figure 1C
fig1c <- df %>% 
  filter(!state %in% c("Puerto Rico", "North Carolina", "Connecticut")) %>%
  group_by(date) %>% 
  summarize(fitted = sum(expected * fitted) / sum(expected), 
            se = sqrt(sum(expected^2 * se^2)) / sum(expected),
            sd = sqrt(sum(expected^2 * sd^2)) / sum(expected)) %>%
  ggplot(aes(date, fitted, ymin = fitted - 2.54*se, ymax = fitted + 2.54*se)) +
  geom_hline(yintercept = 0, lty=2, color="gray") +
  geom_ribbon(aes(date, ymin = -2*sd, ymax = 2*sd), color = 1, fill = NA, lty = 2) +
  geom_ribbon(alpha = 0.5) +
  geom_line() +
  xlab("") +
  ylab("Percent increase from expected mortality") +
  theme(axis.text.y  = element_text(size=18),
        axis.text.x  = element_text(size=18),
        axis.title = element_text(size=17))
fig1c

# -- Saving figure 1c
ggsave("figs/figure-1c.pdf",
       plot   = fig1c,
       dpi    = 300,
       height = 5,
       width  = 7)
### -- -------------------------------- ------------------------------------------------------------------
### -- Figure 1C: Fhat estimate for USA ------------------------------------------------------------------
### -- -------------------------------- ------------------------------------------------------------------
