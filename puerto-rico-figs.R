library(tidyverse)
library(lubridate)
library(excessdeaths)
library(directlabels)

hurricane_dates <- as.Date(c("1989-09-18","1998-09-21","2017-09-20"))
hurricane_effect_ends <- as.Date(c("1990-03-18","1999-03-21","2018-03-20"))
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

control_dates <- seq(as.Date("2002-01-01"), as.Date("2013-12-31"), by = "day")


data("puerto_rico_counts")
the_breaks <- c(0, 5, 20, 40, 60, 75, Inf)
all_counts <- collapse_counts_by_age(puerto_rico_counts, the_breaks)

# Supp Figure - Trend esimtate --------------------------------------------
res <- map_df(unique(all_counts$agegroup), function(x){
  tmp <- all_counts %>% filter(agegroup == x) %>%
    compute_expected(exclude = exclude_dates, keep.components = TRUE)
  data.frame(agegroup = x, date = tmp$counts$date, trend= tmp$trend)
})
library(directlabels)
res %>%
  mutate(agegroup = str_replace(agegroup, "-Inf", "+")) %>%
  ggplot(aes(date, trend, color = agegroup)) +
  geom_line(show.legend = FALSE) +
  geom_dl(aes(color =  agegroup, label =agegroup), method = list("first.points")) +
  scale_x_date(limit = c(make_date(1982,1,1), make_date(2020,1,1))) +
  ylab("Percent of population") +
  facet_wrap(~agegroup, scale = "free_y")

# Supp Figure - Sesonal esimtate --------------------------------------------

res <- map_df(unique(all_counts$agegroup), function(x){
  all_counts %>% filter(agegroup == x) %>%
    compute_expected(exclude = exclude_dates, keep.components = TRUE) %>%
    .$seasonal %>%
  mutate(agegroup = x)
})

res %>%
  mutate(agegroup = str_replace(agegroup, "-Inf", "+")) %>%
  ggplot(aes(day, s, color = agegroup)) +
  geom_line(show.legend = FALSE) +
  geom_dl(aes(color =  agegroup, label = agegroup), method = list("first.points")) +
  scale_x_continuous(limit = c(-5,365))+
  ylab("Percent of population") +
  facet_wrap(~agegroup, scale = "free_y")


# Supp Figure: Correlated errors ------------------------------------------
counts <- filter(all_counts, agegroup == "75-Inf")
counts <- compute_expected(counts, exclude = exclude_dates)
example_dates <- control_dates[year(control_dates) <= 2005]
r <- tibble(date = counts$date, observed = counts$outcome, expected = counts$expected) %>%
  filter(date %in% example_dates) %>%
  mutate(r = (observed - expected)/sqrt(expected)) %>%
  pull(r)

## remake these nice in same ggplot2 style
## A - show that poisson model doesn't quite fit
qqnorm(r); abline(0,1)
# B - there is correlation
acf(r)
## And the SD is not 1
sd(r)

## Now show that we if "prewhiten" the errors things get better
tmp <-counts %>%
  filter(date %in% example_dates) %>%
  mutate(r = (outcome - expected)/expected)
r <- tmp$r
mu <- tmp$expected
n <- length(r)
arfit<- excessdeaths:::fit_ar(tmp,  control.dates = control_dates)
rhos <- ARMAacf(ar = arfit$ar, ma = 0, lag.max = n)
s <- arfit$sigma
mu <- counts$expected[counts$date %in% example_dates]
Sigma <- apply(abs(outer(1:n, 1:n, "-")) + 1, 1, function(i) rhos[i]) *
  outer(sqrt(s^2 + 1/mu), sqrt(s^2 + 1/mu))
V <- chol(Sigma)
V_inv = backsolve(r = V, x = diag(ncol(V))) ## V is upper triangular so backsolve faster

## whitened residuals
qqnorm(V_inv %*% r); abline(0,1)
acf(V_inv %*% r)
sd(V_inv %*% r)


# Figure 1B - Chikingunya -------------------------------------------------
the_breaks <- c(0, 5, 20, 40, 60, Inf)
all_counts <- collapse_counts_by_age(puerto_rico_counts, the_breaks)
fit <- lapply(c("0-4","60-Inf"), function(x){
  ret <- all_counts %>% 
    filter(agegroup == x) %>%
    excess_model(start = make_date(2014, 1, 1),
                 end = make_date(2015, 12, 31),
                 exclude = exclude_dates,
                 control.dates = control_dates,
                 nknots = 12, 
                 model = "correlated")
  return(ret)
})

excess_plot(fit[[1]])
excess_plot(fit[[2]])



# Figure 2 - Excess deaths ------------------------------------------------
the_breaks <- c(0, 5, 20, 40, 60, 75, Inf)
ndays <- 365
all_counts <- collapse_counts_by_age(puerto_rico_counts, the_breaks)
interval_start <- c(hurricane_dates[2],
                  hurricane_dates[3],
                  Chikungunya = make_date(2014, 8, 1),
                  Covid_19 = make_date(2020, 1, 1))

disc <- c(TRUE, TRUE, FALSE, FALSE)
before <-c(365, 365, 365, 548) 
after <-c(365, 365, 365, 90)
npy <- 12

e <- map_df(seq_along(interval_start), function(i){
  tmp <- map_df(unique(all_counts$agegroup), function(x){
    f <- all_counts %>% filter(agegroup == x) %>%
      excess_model(event = interval_start[i],
                   start = interval_start[i] - before[i],
                   end = interval_start[i] + after[i], 
                   exclude = exclude_dates,
                   control.dates = control_dates,
                   nknots = round(npy*(before[i] + after[i])/365), 
                   discontinuity = disc[i],
                   model = "correlated")
    excess_cumulative(f, 
                      start = interval_start[i], 
                      end = pmin(make_date(2020, 3, 31), interval_start[i] + ndays)) %>%
      mutate(agegroup = x, event_day = interval_start[i], event = names(interval_start)[i])
  })
  tmp %>% group_by(date) %>% 
    summarize(fitted = sum(fitted),
              observed = sum(observed),
              sd = sqrt(sum(sd^2)),
              se = sqrt(sum(se^2)),
              event_day = event_day[1], 
              event = event[1])
})

e %>%
  mutate(day = as.numeric(date - event_day)) %>%
  ggplot(aes(color = event, fill = event)) +
  geom_ribbon(aes(day, ymin = fitted - 2*se, ymax = fitted + 2*se), alpha = 0.25) + 
  geom_point(aes(day, observed), alpha = 0.25, cex = 1) +
  geom_line(aes(day, fitted))
  
  
  
  






