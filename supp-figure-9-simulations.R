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

# -- Hurricanes information
hurricane_dates        <- as.Date(c("1989-09-18","1998-09-21","2017-09-20"))
hurricane_effect_ends  <- as.Date(c("1990-03-18","1999-03-21","2018-03-20"))
names(hurricane_dates) <- c("Hugo", "Georges", "Maria")

# -- Control & exclude periods
control_dates <- seq(as.Date("2002-01-01"), as.Date("2013-12-31"), by = "day")
exclude_dates <- c(seq(hurricane_dates[1], hurricane_effect_ends[1], by = "day"),
                   seq(hurricane_dates[2], hurricane_effect_ends[2], by = "day"),
                   seq(hurricane_dates[3], hurricane_effect_ends[3], by = "day"),
                   seq(as.Date("2014-09-01"), as.Date("2015-03-21"), by = "day"),
                   seq(as.Date("2001-01-01"), as.Date("2001-01-15"), by = "day"),
                   seq(as.Date("2020-01-01"), lubridate::today(), by = "day"))

# -- Daily mortality
counts <- puerto_rico_counts %>%
            group_by(date) %>%
            summarize(outcome    = sum(outcome),
                      population = sum(population)) %>%
            ungroup()

# -- Number of knots per year
# nknots <- 6
### -- ---------- ------------------------------------------------------------------
### -- END Set up ------------------------------------------------------------------
### -- ---------- ------------------------------------------------------------------

### -- ----------------------- ------------------------------------------------------------------------
### -- Function simulated data ------------------------------------------------------------------------
### -- ----------------------- ------------------------------------------------------------------------
# -- Function to simulated data
excess_simulation <- function(counts,
                           effect        = 0.60,
                           period_effect = 30,
                           add_years     = 1,
                           event_day     = ymd("1998-09-21"),
                           num_sims      = 1, 
                           nknots        = 10,
                           type          = c("correlated", "overdispersed", "poisson"),
                           model         = "correlated",
                           phi           = 2, 
                           normal        = FALSE,
                           discontinuity = TRUE)
{
  # counts        : A group-specific dataset from compute_counts
  # effect        : theoretical effect size. This should be in [0,1]
  # period_effect : SD of normal density used to generate simualted data.
  #                 This roughly translates to the period of indirect effect
  # add_years     : Deternmines the span around the event day
  # event_day     : Date of disaster
  # num_sims      : Number of simulated datasets
  # discontinuity : If true death counts will jump drastically
  
  # -- Computing expected mortality
  res <- suppressMessages(compute_expected(counts = counts, exclude = exclude_dates))
  
  # -- Make sure event_day is in date format
  event_day <- ymd(event_day)
  
  # -- Parameters
  dates      <- sort(unique(counts$date))
  test_dates <- seq(ymd(event_day)-(365*add_years), ymd(event_day)+(365*add_years), "days")
  test_dates <- test_dates[!(month(test_dates) == 2 & day(test_dates) == 29)]
  
  # -- Normal distribution that is centered at the event date
  f <- dnorm(as.numeric(dates - event_day), 0, period_effect/3)
  
  # -- Adding effect size. This is the true % change
  f <- f*effect/max(f)
  
  # -- Adding a jump if discontinuity == TRUE
  if(discontinuity){
    idx <- which(dates==event_day)
    f[(idx-4*period_effect):(idx-1)] <- 0
  }
  
  # -- If normal==TRUE then no change in percent change
  if(normal){ f <- rep(0, length(f)) }
  
  # -- Computing true excess deaths
  true_excess <- tibble(date = dates, true_excess = res$expected*f)
  # true_excess <- tibble(date = dates, true_excess = res$expected*f) %>%
  #   filter(true_excess != 0)
  
  # -- Simulating data
  message("Simulating data")
  sims    <- lapply(1:num_sims, function(i){
    
    # -- Expected deaths
    mu     <- res$expected
    change <- mu*(1+f)
    
    # -- Simulated dataset
    sim <- counts
    
    if(type == "correlated"){
      
      # -- Syntethic output from ar
      md <- list(order = c(5,0,0), ar = seq(0.1, 0.01, len = 5))
      
      # -- Simulated valus from an ARIMA
      sim_ar <- arima.sim(model = md, 
                          sd = 0.05, 
                          n = length(dates))
      
      # -- Correlated errors (Centered at more or less 1)
      e <- pmax(0, 1 - sim_ar)
      
      # -- Generating Poisson correlated data
      y <- rpois(length(dates), change * e)
    } else if(type == "overdispersed"){
      
      # -- Generating Poisson overdispersed data using NB trick
      y <- rnbinom(length(dates), mu=change, size=(change^2)/(change*(phi-1)))
    } else {
      # -- Generating Poisson data
      y <- rpois(length(dates), change)
    }
    
    # -- Adding y to the dataset
    sim$deaths <- y
    sim$sim    <- i
    return(sim)
  })
  sim_res <- do.call(rbind, sims)
  
  # -- Fitting model to fake data
  message("Fitting excess model")
  fits <- lapply(sims, function(sim){
    
    cat(".")
    
    # -- Temporarily changing name of columns
    sim <- sim %>%
      rename(tmp     = outcome,
             outcome = deaths)
    
    # -- Fitting model
    tmp <- suppressMessages(excess_model(
      counts         = sim,
      event          = event_day,
      start          = first(test_dates),
      end            = last(test_dates),
      model          = model,
      control.dates  = control_dates,
      exclude        = test_dates,#out_dates
      discontinuity  = discontinuity,
      knots.per.year = nknots))
    
    
    # -- End of effect period
    end_day <- with(tmp, 
                    tibble(date, fitted, se) %>%
                      filter(date >= event_day) %>%
                      filter(fitted-1.96*se <= 0) %>%
                      slice(1) %>%
                      .$date)
    # end_day <- tmp$detected_intervals$end[nrow(tmp$detected_intervals)]
    # if(is.na(end_day)){end_day <- event_day + 30}
    # if(end_day < event_day){end_day <- event_day + 1}
    
    # -- Index for period of indirect effect
    # idx_period <- seq(as.Date(event_day), as.Date(end_day), by=1)
    idx_period <- seq(as.Date(event_day, format="%y-%m-%d"), as.Date(end_day, format="%y-%m-%d"), by="1 day")
    
    # -- True % change in the period of indirect effect
    true_change <- f[which(dates %in% idx_period)]
    
    # -- Getting expected counts for the period of indirect effect
    tmp_expected <- tmp$expected[which(tmp$date %in% idx_period)]
    
    # -- Adding last day of period effect
    tmp$end_day <- end_day
    
    return(tmp)
  })
  
  message("\nCreating figures")
  
  # -- Finding correct dates
  index <- match(fits[[1]]$date, dates)
  
  # -- Getting expected, fitted, and se estimates for each iteration
  fitted <- do.call(rbind, lapply(1:length(fits), function(i){ 
    tibble(date=fits[[i]]$date, expected=fits[[i]]$expected, fitted=fits[[i]]$fitted, se=fits[[i]]$se, end_day = fits[[i]]$end_day, sim=i)
  })) %>%
    left_join(tibble(date = dates[index], true_f = f[index]), by="date")
  
  # -- Computing avg and sd of fitted values 
  avg_fitted <- fitted %>% group_by(date) %>% summarize(se = sd(fitted), fitted = mean(fitted))
  
  # -- Viz effect size
  fig_effect_size <- ggplot() +
    geom_line(aes(date, 100*fitted, group=sim), alpha=0.10, data=filter(fitted, sim!=1)) +
    geom_line(aes(date, 100*fitted, group=sim, color="Estimate"), alpha=0.50, data=filter(fitted, sim==1)) +
    geom_line(aes(date, 100*fitted, color="black"), lty=2, data=avg_fitted) +
    geom_line(aes(date, 100*true_f, color="Truth"), size=1, data=fitted) +
    scale_color_manual(name="",
                       values=c("black", "gray", "red2"),
                       labels=c("Avg.", "Estimate", "Truth")) +
    ylab("% change in mortality") +
    xlab("") +
    ggtitle("True vs. Estimated Effect Size") +
    theme(axis.title = element_text(face="bold"),
          axis.text  = element_text(face="bold"),
          legend.position = "bottom",
          legend.text     = element_text(face="bold"))
  
  # -- Viz bias
  fig_bias <- ggplot() +
    geom_hline(yintercept = 0, lty=2) +
    geom_line(aes(date, f[index] - fitted, group=sim), alpha=0.10, data=filter(fitted, sim!=1)) +
    geom_line(aes(date, f[index] - fitted), size=1, color="red2", data=avg_fitted) +
    ylab("Bias") +
    xlab("") +
    ggtitle("Bias: True - Estiamted Effect") +
    theme(axis.title = element_text(face="bold"),
          axis.text  = element_text(face="bold"),
          legend.position = "bottom",
          legend.text     = element_text(face="bold"))  
  
  # -- Variability estimates
  ses <- fitted %>%
    group_by(date) %>%
    mutate(sed = sd(fitted)) %>%
    ungroup()
  
  # -- Viz variability
  fig_se <- ggplot() +
    geom_line(aes(date, 100*se, group = sim), alpha=0.50, data = filter(ses, sim!=1)) +
    geom_line(aes(date, 100*se, color="Estimated standard erros"), alpha=0.50, data = filter(ses, sim==1)) +
    geom_line(aes(date, 100*sed, color="Standard deviation of fits"), size=1, data=filter(ses, sim==1)) +
    ylab("Variability estimates") +
    xlab("") +
    scale_color_manual(name="",
                       values=c("black", "red2")) +
    ggtitle("Asssesing Variability Estimates") +
    theme(axis.title = element_text(face="bold"),
          axis.text  = element_text(face="bold"),
          legend.position = "bottom",
          legend.text     = element_text(face="bold"))
  
  # -- Computing excess deaths
  excess_deaths <- map_df(seq_along(fits), function(i){
    
    fit     <- fits[[i]]
    end_day <- fit$end_day
    
    
    excess <- excess_cumulative(fit, event_day, end_day) %>%
      as_tibble() %>%
      mutate(days = as.numeric(date-min(date)+1)) %>%
      select(days, observed, fitted, sd, se) %>%
      bind_rows(tibble(days=0, observed=0, fitted=0, sd=0, se=0)) %>%
      arrange(days) %>%
      mutate(sim = i)
  })
  
  # -- True excess deaths
  true_excess <- true_excess %>%
    filter(date >= event_day) %>%
    mutate(days = as.numeric(date - min(date) + 1)) %>%
    select(days, true_excess) %>%
    mutate(true_excess = cumsum(true_excess)) %>%
    bind_rows(tibble(days=0, true_excess=0)) %>%
    arrange(days)
  
  # -- Average excess deaths
  avg_excess <- excess_deaths %>%
    group_by(days) %>%
    summarize(avg_fitted    = mean(fitted),
              avg_observed  = mean(observed),
              se.observed = sd(observed),
              se.fitted   = sd(fitted),
              mean.se.observed = mean(sd),
              mean.se.fitted   = mean(se)) %>%
    ungroup() %>%
    left_join(true_excess, by="days")
  
  # -- Just for output purposes
  true_excess <- avg_excess %>%
    select(days, true_excess)
  
  # -- Viz excess deaths
  fig_excess <- ggplot() +
    geom_vline(xintercept = period_effect, color="#525252", alpha=0.70) +
    geom_step(aes(days, observed, group=sim), color="gray", alpha=0.50, data=excess_deaths) +
    geom_step(aes(days, fitted, group=sim), color="gray", alpha=0.50, data=excess_deaths) +
    geom_step(aes(days, avg_observed, color="Avg. Obs"),size=0.70, data=avg_excess) +
    geom_step(aes(days, avg_fitted, color="Avg. Fit"), size=0.70, data=avg_excess) +
    geom_step(aes(days, true_excess, color="Truth"), size=0.80, data=avg_excess) +
    scale_color_manual(name="",
                       values=c("#6baed6", "#fb6a4a", "black")) +
    ylab("Cumulative excess deaths") +
    xlab("Days since the event day") +
    labs(caption = "Vertical line corresponds to true period") +
    ggtitle("True vs. Estimated Excess Deaths") +
    theme(axis.title = element_text(face="bold"),
          axis.text  = element_text(face="bold"),
          legend.position = "bottom",
          legend.text     = element_text(face="bold"))
  
  # -- Viz variability estimates for excess deaths
  fig_excess_se <- avg_excess %>%
    ggplot() +
    geom_vline(xintercept = period_effect, color="#525252", alpha=0.70) +
    geom_line(aes(days, se.observed, color="Observed", lty="1")) +
    geom_line(aes(days, mean.se.observed, color="Observed", lty="2")) +
    geom_line(aes(days, se.fitted, color="Fitted", lty="1")) +
    geom_line(aes(days, mean.se.fitted, color="Fitted", lty="2")) +
    ylab("Variability estimates") +
    xlab("") +
    labs(caption = "Vertical line corresponds to true period") +
    scale_color_manual(name="",
                       values=c("#6baed6", "#fb6a4a")) +
    scale_linetype_manual(name="",
                          values = 1:2,
                          labels   = c("True", "Estimate")) +
    ggtitle("Asssesing Variability Estimates") +
    theme(axis.title = element_text(face="bold"),
          axis.text  = element_text(face="bold"),
          legend.position = "bottom",
          legend.text     = element_text(face="bold"))
  
  return(list("true_excess"     = true_excess, 
              "sim"             = sim_res, 
              "fitted"          = fitted, 
              "avg_fitted"      = avg_fitted, 
              "excess_deaths"   = excess_deaths, 
              "fig_effect_size" = fig_effect_size,
              "fig_bias"        = fig_bias,
              "fig_se"          = fig_se,
              "fig_excess"      = fig_excess,
              "fig_excess_se"   = fig_excess_se))
}
### -- --------------------------- ------------------------------------------------------------------------
### -- END Function simulated data ------------------------------------------------------------------------
### -- --------------------------- ------------------------------------------------------------------------

### -- ------------------------------- ------------------------------------------------------------------------
### -- Supp figure 9: Simulation study ------------------------------------------------------------------------
### -- ------------------------------- ------------------------------------------------------------------------
set.seed(1)
# -- Scenario 1: Hurricane like event
hurricane <- excess_simulation(counts        = counts,
                               effect        = 0.75, 
                               period_effect = 150,
                               num_sims      = 100,
                               nknots        = 6,
                               model         = "correlated",
                               type          = "correlated",
                               discontinuity = TRUE)

# -- Scenario 2: Epidemic like event
epidemic <- excess_simulation(counts        = counts,
                              effect        = 0.35, 
                              period_effect = 75,
                              num_sims      = 100,
                              nknots        = 10,
                              model         = "correlated",
                              type          = "correlated",
                              discontinuity = FALSE)

# -- Scenario 3: No event
normal <- excess_simulation(counts        = counts,
                            effect        = 0.35, 
                            period_effect = 75,
                            num_sims      = 100,
                            nknots        = 6,
                            model         = "correlated",
                            type          = "correlated",
                            normal        = TRUE,
                            discontinuity = FALSE)

# -- Supp figure 9a
supp_fig9a <- hurricane$fitted %>%
  filter(date >= "1998-04-01", date <= "1999-04-01") %>%
  ggplot(aes(date, 100*fitted, group=sim)) +
  geom_line(color="gray", alpha=0.50) +
  geom_line(aes(date, 100*true_f)) +
  xlab("") +
  ylab("Percent increase from expected mortality") +
  scale_x_date(date_breaks = "1 month", 
               date_labels = "%b") +
  theme(axis.text  = element_text(size=12),
        axis.title = element_text(size=13))

# -- Save Supp figure 9a
ggsave("figs/supp-figure-9a.pdf",
       plot   = supp_fig9a,
       dpi    = 300, 
       height = 6,
       width  = 8)

# -- Supp figure 9b
supp_fig9b <- epidemic$fitted %>%
  filter(date >= "1998-04-01", date <= "1999-04-01") %>%
  ggplot(aes(date, 100*fitted, group=sim)) +
  geom_line(color="gray", alpha=0.50) +
  geom_line(aes(date, 100*true_f)) +
  xlab("") +
  ylab("Percent increase from expected mortality") +
  scale_x_date(date_breaks = "1 month", 
               date_labels = "%b") +
  theme(axis.text  = element_text(size=12),
        axis.title = element_text(size=13))

# -- Save Supp figure 9b
ggsave("figs/supp-figure-9b.pdf",
       plot   = supp_fig9b,
       dpi    = 300, 
       height = 6,
       width  = 8)

# -- Supp figure 9c
supp_fig9c <- normal$fitted %>%
  filter(date >= "1998-04-01", date <= "1999-04-01") %>%
  ggplot(aes(date, 100*fitted, group=sim)) +
  geom_line(color="gray", alpha=0.50) +
  geom_line(aes(date, 100*true_f)) +
  xlab("") +
  ylab("Percent increase from expected mortality") +
  scale_x_date(date_breaks = "1 month", 
               date_labels = "%b") +
  theme(axis.text  = element_text(size=12),
        axis.title = element_text(size=13))

# -- Save Supp figure 9c
ggsave("figs/supp-figure-9c.pdf",
       plot   = supp_fig9c,
       dpi    = 300, 
       height = 6,
       width  = 8)

# -- Supp figure 9d
supp_fig9d <-hurricane$fitted %>%
  filter(date >= "1998-04-01", date <= "1999-04-01") %>%
  group_by(date) %>%
  mutate(sd_fitted = sd(fitted)) %>%
  ungroup() %>%
  ggplot(aes(date, 100*se, group = sim)) +
  geom_line(color="gray", alpha=0.50) +
  geom_line(aes(date, 100*sd_fitted)) +
  scale_x_date(date_breaks = "1 month", 
               date_labels = "%b") +
  ylab("Variability estimate") +
  xlab("") +
  theme(axis.text  = element_text(size=12),
        axis.title = element_text(size=13))

# -- Save Supp figure 9d
ggsave("figs/supp-figure-9d.pdf",
       plot   = supp_fig9d,
       dpi    = 300, 
       height = 6,
       width  = 8)

# -- Supp figure 9e
supp_fig9e <- epidemic$fitted %>%
  filter(date >= "1998-04-01", date <= "1999-04-01") %>%
  group_by(date) %>%
  mutate(sd_fitted = sd(fitted)) %>%
  ungroup() %>%
  ggplot(aes(date, 100*se, group = sim)) +
  geom_line(color="gray", alpha=0.50) +
  geom_line(aes(date, 100*sd_fitted)) +
  scale_x_date(date_breaks = "1 month", 
               date_labels = "%b") +
  ylab("Variability estimate") +
  xlab("") +
  theme(axis.text  = element_text(size=12),
        axis.title = element_text(size=13))

# -- Save Supp figure 9e
ggsave("figs/supp-figure-9e.pdf",
       plot   = supp_fig9e,
       dpi    = 300, 
       height = 6,
       width  = 8)

# -- Supp figure 9f
supp_fig9f <- normal$fitted %>%
  filter(date >= "1998-04-01", date <= "1999-04-01") %>%
  group_by(date) %>%
  mutate(sd_fitted = sd(fitted)) %>%
  ungroup() %>%
  ggplot(aes(date, 100*se, group = sim)) +
  geom_line(color="gray", alpha=0.50) +
  geom_line(aes(date, 100*sd_fitted)) +
  scale_x_date(date_breaks = "1 month", 
               date_labels = "%b") +
  ylab("Variability estimate") +
  xlab("") +
  theme(axis.text  = element_text(size=12),
        axis.title = element_text(size=13))

# -- Save Supp figure 9f
ggsave("figs/supp-figure-9f.pdf",
       plot   = supp_fig9f,
       dpi    = 300, 
       height = 6,
       width  = 8)
### -- ----------------------------------- ------------------------------------------------------------------------
### -- END Supp figure 9: Simulation study ------------------------------------------------------------------------
### -- ----------------------------------- ------------------------------------------------------------------------