### -- Set up --------------------------------------------------------

# -- Libraries
library(tidyverse)
library(lubridate)
library(excessdeaths)
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
nknots <- 6

# -- Loop to fit models
tmp <- map_df(seq_along(count_index), function(i){
  
  # -- Current index
  h <- count_index[i]
  
  if(i < 4){ # -- Here we fit model to hurricanes outside of PR
    
    # -- Fitting model
    fit <-  excess_model(counts[[h]], 
                         event         = hurricane_dates[[i]],
                         start         = hurricane_dates[[i]] - before,
                         end           = hurricane_dates[[i]] + after, 
                         exclude       = exclude_dates[[i]],
                         control.dates = control_dates[[i]],
                         nknots        = nknots,
                         model         = "correlated")
    
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
        excess_model(event         = hurricane_dates[[i]],
                     start         = hurricane_dates[[i]] - before,
                     end           = hurricane_dates[[i]] + after, 
                     exclude       = exclude_dates[[i]],
                     control.dates = control_dates[[i]],
                     nknots        = nknots, 
                     model         = "correlated")
      
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

# -- Figure 1
fig1 <- tmp %>% 
  group_by(hurricane) %>%
  mutate(date = seq(make_date(2017,07,01), make_date(2018,07,02), by="days")) %>%
  ungroup() %>%
  mutate(day       = as.numeric(date - hurricane_date),
         hurricane = factor(hurricane, levels=rev(unlist(hurricane_names)))) %>%
  ggplot(aes(date, fitted*100, color = hurricane)) +
  geom_hline(yintercept = 0, lty=2, color="gray") +
  geom_line(size=1) +
  xlab("") +
  ylab("Percent change in rate of mortality") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  scale_y_continuous(limits = c(-5, 75),
                     breaks = seq(0, 75, by=10)) +
  scale_color_manual(name="",
                     values = c("#D55E00","#0571b0","#009E73","#56B4E9","#CC79A7","#E69F00","#ca0020","gray")) +
  theme(axis.title = element_text(face="bold", color="black"),
        axis.text  = element_text(face="bold", color="black"),
        legend.title = element_text(face="bold", color="black"),
        legend.text  = element_text(face="bold", color="black"))
fig1

# -- Save figure 1
ggsave("figs/figure-1.pdf",
       plot   = fig1,
       dpi    = 300, 
       height = 6,
       width  = 8)

# -- Supp Figure 1
supp_fig1 <- tmp %>% 
  mutate(hurricane = factor(hurricane, levels=rev(unlist(hurricane_names)))) %>%
  ggplot(aes(date, fitted*100)) +
  geom_hline(yintercept = 0, lty=2, color="gray") +
  geom_ribbon(aes(ymin = fitted*100 - 2*se*100, ymax = fitted*100 + 2*se*100), alpha = 0.5, color = NA) +
  geom_line() +
  xlab("") +
  ylab("Percent change in rate of mortality") +
  facet_wrap( ~ hurricane, scale= "free_x") +
  scale_x_date(date_breaks = "3 month", date_labels = "%b") +
  theme(axis.title = element_text(face="bold", color="black"),
        axis.text  = element_text(face="bold", color="black"),
        strip.text   = element_text(face="bold", color="black"),
        legend.title = element_text(face="bold", color="black"),
        legend.text  = element_text(face="bold", color="black"))

# -- Save supp figure 1
ggsave("figs/supp-figure-1.pdf",
       plot   = supp_fig1,
       dpi    = 300,
       height = 6,
       width  = 8)
