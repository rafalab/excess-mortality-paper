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

# -- Number of knots per year
nknots <- 6
### -- ---------- ------------------------------------------------------------------
### -- END Set up ------------------------------------------------------------------
### -- ---------- ------------------------------------------------------------------

### -- ------------------------------------------- ------------------------------------------------------------------
### -- Supp figure 2: Individual hurricane effects ------------------------------------------------------------------
### -- ------------------------------------------- ------------------------------------------------------------------
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

# -- Supp Figure 2
supp_fig <- tmp %>% 
  mutate(day       = as.numeric(date - hurricane_date),
         hurricane = factor(hurricane, levels=rev(unlist(hurricane_names)))) %>%
  filter(day >= -125, day <= 245) %>%
  ggplot(aes(date, fitted*100)) +
  geom_hline(yintercept = 0, lty=2, color="gray") +
  geom_ribbon(aes(ymin = fitted*100 - 2*se*100, ymax = fitted*100 + 2*se*100), alpha = 0.5, color = NA) +
  geom_line() +
  xlab("Date") +
  ylab("Percent increase from expected mortality") +
  facet_wrap( ~ hurricane, scale= "free_x") +
  scale_x_date(date_breaks = "3 month", date_labels = "%b") +
  theme(axis.text  = element_text(size=15),
        axis.title = element_text(size=17))

# -- Save supp figure 1
ggsave("figs/supp-figure-3-hurricane-effects.pdf",
       plot   = supp_fig,
       dpi    = 300,
       height = 5,
       width  = 7)
### -- ----------------------------------------------- ------------------------------------------------------------------
### -- END Supp figure 2: Individual hurricane effects ------------------------------------------------------------------
### -- ----------------------------------------------- ------------------------------------------------------------------