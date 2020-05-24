### -- ------ ------------------------------------------------------------------
### -- Set up ------------------------------------------------------------------
### -- ------ ------------------------------------------------------------------
# -- Libraries
library(ggpubr)
library(tidyverse)
library(lubridate)
library(excessmort)
dslabs::ds_theme_set()

# -- Loading data
data("puerto_rico_counts")
data("florida_counts")
data("new_jersey_counts")
data("louisiana_counts")

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
puerto_rico_out_dates <- c(
  seq(puerto_rico_hurricane_dates[1], puerto_rico_hurricane_effect_ends[1], by = "day"),
  seq(puerto_rico_hurricane_dates[2], puerto_rico_hurricane_effect_ends[2], by = "day"),
  seq(puerto_rico_hurricane_dates[3], puerto_rico_hurricane_effect_ends[3], by = "day"),
  seq(as.Date("2014-09-01"), as.Date("2015-03-21"), by = "day"),
  seq(as.Date("2001-01-01"), as.Date("2001-01-15"), by = "day"),
  seq(as.Date("2020-01-01"), lubridate::today(), by = "day"))

# -- Dates to exclude for each hurricane
exclude_dates  <- list(irma    = hurricane_dates[["irma"]] + 0:180,
                       katrina = hurricane_dates[["katrina"]]  + 0:180,
                       sandy   = hurricane_dates[["sandy"]]  + 0:180,
                       hugo    = puerto_rico_out_dates,
                       georges = puerto_rico_out_dates,
                       maria   = puerto_rico_out_dates)

# -- To be used as parameters in model fitting below
before <- days(122)
after  <- days(244)

# -- Number of knots per year
nknots <- 6
### -- ---------- ------------------------------------------------------------------
### -- END Set up ------------------------------------------------------------------
### -- ---------- ------------------------------------------------------------------

### -- ---------------------------------------- ------------------------------------------------------------------
### -- Figure 1A: Fhat estimates for hurricanes ------------------------------------------------------------------
### -- ---------------------------------------- ------------------------------------------------------------------
# -- Loop to fit models to hurricanes
tmp <- map_df(seq_along(count_index), function(i){
  
  # -- Current index
  h <- count_index[i]
  
  if(i < 4){ # -- Here we fit model to hurricanes outside of PR
    
    # -- Fitting model
    fit <-  excess_model(counts[[h]], 
                         event          = hurricane_dates[[i]],
                         start          = hurricane_dates[[i]] - before,
                         end            = hurricane_dates[[i]] + after, 
                         exclude        = exclude_dates[[i]],
                         control.dates  = control_dates[[i]],
                         knots.per.year = nknots,
                         model          = "correlated")
    
    # -- Dataframe with fit results
    ret <- tibble(date           = fit$date, 
                  fitted         = fit$fitted, 
                  se             = fit$se, 
                  hurricane      = hurricane_names[[i]],
                  hurricane_date = hurricane_dates[[i]])
    
  } else { # -- Here we fit model to hurricanes in PR
    
    # -- This loops through the demographic groups for PR
    tmp <- map_df(unique(counts[[h]]$agegroup), function(x){
      
      # -- Fitting model to each demographic group
      fit <- counts[[h]] %>% 
        filter(agegroup == x) %>%
        excess_model(event          = hurricane_dates[[i]],
                     start          = hurricane_dates[[i]] - before,
                     end            = hurricane_dates[[i]] + after, 
                     exclude        = exclude_dates[[i]],
                     control.dates  = control_dates[[i]],
                     knots.per.year = nknots, 
                     model          = "correlated")
      
      # -- Dataframe with fit results for each demographics
      tibble(date = fit$date, expected = fit$expected, fitted = fit$fit, se = fit$se)
    })
    
    # -- Putting PR results together
    ret <- tmp %>% 
      group_by(date) %>% 
      summarize(fitted = sum(expected * fitted) / sum(expected), 
                se     = sqrt(sum(expected^2 * se^2)) / sum(expected)) %>%
      mutate(hurricane      = hurricane_names[[i]],
             hurricane_date = hurricane_dates[[i]])
  }
  return(ret)
}) %>%
  mutate(hurricane = reorder(hurricane, date, min))

# -- Figure 1A
fig1a <- tmp %>% 
  group_by(hurricane) %>%
  mutate(date = seq(make_date(2017,07,01), make_date(2018,07,02), by="days")) %>%
  ungroup() %>%
  mutate(day       = as.numeric(date - hurricane_date),
         hurricane = factor(hurricane, levels=rev(unlist(hurricane_names)))) %>%
  ggplot(aes(date, fitted*100, color = hurricane)) +
  geom_hline(yintercept = 0, lty=2, color="gray") +
  geom_line(size=1, show.legend = F) +
  geom_dl(aes(color=hurricane, label=hurricane), 
          method=list(fontface="bold", "smart.grid")) +#"last.qp
  xlab("") +
  ylab("Percent change in mortality") +
  scale_x_date(date_breaks = "1 month", 
               date_labels = "%b",
               limits = c(ymd("2017-07-01"), ymd("2018-08-01"))) +
  scale_y_continuous(limits = c(-5, 75),
                     breaks = seq(0, 75, by=10)) +
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
fig1a

# -- Save figure 1A
ggsave("figs/figure-1a.pdf",
       plot   = fig1a,
       dpi    = 300, 
       height = 8.50,
       width  = 11.0)
### -- -------------------------------------------- ------------------------------------------------------------------
### -- END Figure 1A: Fhat estimates for hurricanes ------------------------------------------------------------------
### -- -------------------------------------------- ------------------------------------------------------------------

### -- ----------------------------------------- ------------------------------------------------------------------
### -- Figure 1B: Fhat estimates for Chikungunya ------------------------------------------------------------------
### -- ----------------------------------------- ------------------------------------------------------------------
# -- Creating breaks for age groups and collapsing data
the_breaks <- c(0, 5, 20, 40, 60, Inf)
all_counts <- collapse_counts_by_age(puerto_rico_counts, the_breaks)

# -- Determining period of interest
start <- lubridate::make_date(2014, 01, 01)
end   <- lubridate::make_date(2015, 12, 31)
chick_exclude <- unique(c(exclude_dates$maria, seq(start, end, by="days")))

# -- Fitting model only to the groups of interest
res <- map_df(c("0-4", "60-Inf"), function(x){
  
  tmp_counts <- filter(all_counts, agegroup == x)
  tmp_fit    <- excess_model(counts         = tmp_counts, 
                             start          = start,
                             end            = end,
                             exclude        = chick_exclude,
                             control.dates  = control_dates$maria,
                             knots.per.year = nknots,
                             model          = "correlated",
                             discontinuity  = FALSE)
  
  if(x == "0-4") { flag <- "0 to 4 years" } else {flag <- "Above 65 years" }
  tibble(date = tmp_fit$date, fhat = tmp_fit$fitted, se = tmp_fit$se, agegroup = flag)
})

# -- Figure 1B
fig1b <- res %>%
  filter(date >= "2014-05-01", date <= "2015-05-01") %>%
  ggplot(aes(date, 100*fhat, color=agegroup)) +
  geom_hline(yintercept = 0, lty=2, color="gray") +
  geom_line(size=1, show.legend = F) +
  xlab("") +
  ylab("Percent change in mortality") +
  geom_dl(aes(color=agegroup, label=agegroup), 
          method=list(fontface="bold", "smart.grid")) +#"last.qp"
  scale_x_date(date_breaks = "1 month", 
               date_labels = "%b") +
  scale_y_continuous(limits = c(-80, 70),
                     breaks = seq(-80, 70, by=20)) +
  scale_color_manual(name="",
                     values = c("#D55E00","#0571b0","#009E73","#56B4E9","#CC79A7","#E69F00","#ca0020","gray")) +
  theme(axis.title = element_text(face="bold", color="black"),
        axis.text  = element_text(face="bold", color="black"),
        legend.title    = element_blank(),
        legend.text     = element_text(face="bold", color="black"),
        legend.background = element_rect(color    ="black",
                                         fill     ="white",
                                         linetype = "solid"),
        legend.position = c(0.50, 0.05),
        legend.key.size = unit(0.1, "cm"),
        legend.direction = "horizontal")
fig1b

# -- Save figure 1B
ggsave("figs/figure-1b.pdf",
       plot   = fig1b,
       dpi    = 300, 
       height = 8.50,
       width  = 11.0)
### -- --------------------------------------------- ------------------------------------------------------------------
### -- END Figure 1B: Fhat estimates for Chikungunya ------------------------------------------------------------------
### -- --------------------------------------------- ------------------------------------------------------------------

### -- -------------------------------- ------------------------------------------------------------------
### -- Figure 1C: Fhat estimate for USA ------------------------------------------------------------------
### -- -------------------------------- ------------------------------------------------------------------
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

# -- Figure 1C
fig1c <- df %>% 
  filter(!state %in% c("North Carolina", "Connecticut")) %>%
  group_by(date) %>% 
  summarize(fitted = sum(expected * fitted) / sum(expected), 
            se = sqrt(sum(expected^2 * se^2)) / sum(expected),
            sd = sqrt(sum(expected^2 * sd^2)) / sum(expected)) %>%
  ggplot(aes(date, fitted, ymin = fitted - 2.54*se, ymax = fitted + 2.54*se)) +
  geom_hline(yintercept = 0, lty=2, color="gray") +
  geom_ribbon(aes(date, ymin = -2*sd, ymax = 2*sd), color = 1, fill = NA, lty = 2) +
  geom_ribbon(alpha = 0.5) +
  geom_line(size=1) +
  xlab("") +
  ylab("Percent change in mortality") +
  scale_x_date(date_breaks = "4 month", date_labels = "%b %y") +
  theme(axis.title = element_text(face="bold", color="black"),
        axis.text  = element_text(face="bold", color="black"),
        legend.title    = element_text(face="bold", color="black"),
        legend.text     = element_text(face="bold", color="black"),
        legend.background = element_rect(color    ="black",
                                         fill     ="white",
                                         linetype = "solid"),
        legend.position = c(0.90,0.80))
fig1c

# -- Saving figure c
ggsave("figs/figure-1c.pdf",
       plot   = fig1c,
       dpi    = 300, 
       height = 8.50,
       width  = 11.0)
### -- -------------------------------- ------------------------------------------------------------------
### -- Figure 1C: Fhat estimate for USA ------------------------------------------------------------------
### -- -------------------------------- ------------------------------------------------------------------

# -- Figure 1
fig1 <- ggarrange(fig1a, fig1b, fig1c, 
                  labels = c("A", "B", "C"), 
                  nrow = 3, ncol = 1)

# -- Saving figure 1
ggsave("figs/figure-1.pdf",
       plot   = fig1,
       dpi    = 300, 
       height = 8.50,
       width  = 11.0)
