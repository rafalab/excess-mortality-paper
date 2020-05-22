# -- Libraries
library(ggpubr)
library(tidyverse)
library(lubridate)
library(excessdeaths)
library(directlabels)
library(directlabels)

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

### -- ------------------------------------------ -----------------------------------------------------
### -- Supp Figure 1 (Trend & Seasonal Estimates) -----------------------------------------------------
### -- ------------------------------------------ -----------------------------------------------------
# -- Getting trend estimates for each age group
res <- map_df(unique(all_counts$agegroup), function(x){
  
  tmp <- all_counts %>% 
    filter(agegroup == x) %>%
    compute_expected(exclude = exclude_dates, keep.components = TRUE)
  
  data.frame(agegroup = x, date = tmp$counts$date, trend= tmp$trend)
})

# -- Trend estimates
supp_fig1a <- res %>%
  mutate(agegroup = str_replace(agegroup, "-Inf", "+")) %>%
  mutate(agegroup = factor(agegroup, levels = c("0-4", "5-19", "20-39", 
                                                "40-59", "60-74", "75+"))) %>%
  ggplot(aes(date, trend)) +
  geom_line(size=1, show.legend = FALSE) +
  scale_x_date(limit = c(make_date(1985,1,1), make_date(2020,1,1))) +
  ylab("Estimated trend effect") + 
  xlab("") +
  facet_wrap(~agegroup, scale = "free_y") +
  theme(axis.title = element_text(face="bold", color="black"),
        axis.text  = element_text(face="bold", color="black"),
        strip.text   = element_text(face="bold", color="black"),
        legend.title = element_text(face="bold", color="black"),
        legend.text  = element_text(face="bold", color="black"))

# -- Save supp-figure-1a
ggsave("figs/supp-figure-1a.pdf",
       plot   = supp_fig1a,
       dpi    = 300, 
       height = 6,
       width  = 8)

# -- Getting seasonal estimates for each age group
res <- map_df(unique(all_counts$agegroup), function(x){
  
  all_counts %>% 
    filter(agegroup == x) %>%
    compute_expected(exclude = exclude_dates, keep.components = TRUE) %>%
    .$seasonal %>%
  mutate(agegroup = x)
})

# -- Seasonal estimates
supp_fig1b <- res %>%
  mutate(agegroup = str_replace(agegroup, "-Inf", "+")) %>%
  mutate(agegroup = factor(agegroup, levels = c("0-4", "5-19", "20-39", 
                                                "40-59", "60-74", "75+"))) %>%
  ggplot(aes(day, s)) +
  geom_line(size=1, show.legend = FALSE) +
  scale_x_continuous(limit = c(-5,365))+
  ylab("Estimated seasonal effect") + 
  xlab("Day of the year") +
  facet_wrap(~agegroup, scale = "free_y") +
  theme(axis.title = element_text(face="bold", color="black"),
        axis.text  = element_text(face="bold", color="black"),
        strip.text   = element_text(face="bold", color="black"),
        legend.title = element_text(face="bold", color="black"),
        legend.text  = element_text(face="bold", color="black"))

# -- Save supp-figure-1b
ggsave("figs/supp-figure-1b.pdf",
       plot   = supp_fig1b,
       dpi    = 300, 
       height = 6,
       width  = 8)

# -- Supp-figure-1
supp_fig1 <- ggarrange(supp_fig1a, supp_fig2a, 
                       labels = c("A", "B"), 
                       nrow = 2, ncol = 1)

# -- Save supp-figure-1
ggsave("figs/supp-figure-1.pdf",
       plot   = supp_fig1,
       dpi    = 300, 
       height = 6,
       width  = 8)
### -- ---------------------------------------------- -----------------------------------------------------
### -- END Supp Figure 1 (Trend & Seasonal Estimates) -----------------------------------------------------
### -- ---------------------------------------------- -----------------------------------------------------

### -- -------------------------------- -----------------------------------------------------
### -- Supp Figure 2 (Correlated errors) -----------------------------------------------------
### -- -------------------------------- -----------------------------------------------------
# -- Mortality data for individuals 75 years and older
# counts <- filter(all_counts, agegroup == "75-Inf")
counts <- puerto_rico_counts %>%
            group_by(date) %>%
            summarize(outcome    = sum(outcome),
                      population = sum(population)) %>%
            ungroup()

# -- Computing expected mortality counts
counts <- compute_expected(counts, exclude = exclude_dates)

# -- Dates for model check
example_dates <- control_dates[year(control_dates) <= 2005]

# -- Computing z scores for example dates (H0 true)
r <- tibble(date = counts$date, observed = counts$outcome, expected = counts$expected) %>%
  filter(date %in% example_dates) %>%
  mutate(r = (observed - expected)/sqrt(expected)) %>%
  pull(r)

# -- Supp Figure 2A: Show that poisson model doesn't quite fit
# qqnorm(r); abline(0,1): OLD VERSION
supp_fig2a <- tibble(r=r) %>%
  ggplot(aes(sample=r)) +
  stat_qq(alpha=0.50) + 
  geom_abline(intercept = 0, slope = 1, color="red", lty=2) +
  ylab("Sample quantiles") +
  xlab("Theoretical quantiles") +
  theme(axis.title = element_text(face="bold"),
        axis.text  = element_text(face="bold"))

# -- Save supp-figure-2a
ggsave("figs/supp-figure-2a.pdf",
       plot   = supp_fig2a,
       dpi    = 300, 
       height = 6,
       width  = 8)

# -- Supp Figure 2B: There is correlation
# acf(r): OLD VERSION
auto_cor <- acf(r, plot=FALSE)
supp_fig2b <- tibble(acf = auto_cor$acf, lag = auto_cor$lag) %>%
  ggplot(aes(lag, acf)) +
  geom_col(color="black", fill="#252525", width = 0.5) +
  ylab("ACF") +
  xlab("Lag") +
  geom_hline(yintercept = qnorm((1+0.95)/2)/sqrt(auto_cor$n.used), 
             colour = "#cb181d",
             linetype = 2) +
  geom_hline(yintercept = -qnorm((1+0.95)/2)/sqrt(auto_cor$n.used), 
             colour = "#cb181d",
             linetype = 2) +
  theme(axis.title = element_text(face="bold"),
        axis.text  = element_text(face="bold"))

# -- Save supp-figure-2b
ggsave("figs/supp-figure-2b.pdf",
       plot   = supp_fig2b,
       dpi    = 300, 
       height = 6,
       width  = 8)

# -- The SD is not 1
sd(r)

# -- Now show that we if "prewhiten" the errors things get better
tmp <- counts %>%
        filter(date %in% example_dates) %>%
        mutate(r = (outcome - expected)/expected)

# -- Extracting percent change
r     <- tmp$r

# -- Expected value
mu    <- tmp$expected

# -- Number of data points
n     <- length(r)

# -- Estimating AR(p) with control dates
arfit <- excessdeaths:::fit_ar(tmp,  control.dates = control_dates)

# -- Estimated rhos
rhos  <- ARMAacf(ar = arfit$ar, ma = 0, lag.max = n)

# -- Estimated sd
s     <- arfit$sigma

# -- Computing covariance matrix
Sigma <- apply(abs(outer(1:n, 1:n, "-")) + 1, 1, function(i) rhos[i]) * outer(sqrt(s^2 + 1/mu), sqrt(s^2 + 1/mu))

# -- Cholesky decomposition of Sigma & computing inverse
V     <- chol(Sigma)
V_inv <- backsolve(r = V, x = diag(ncol(V))) ## V is upper triangular so backsolve faster

# -- Supp Figure 2C: QQplot of adjusted residuals
# qqnorm(V_inv %*% r); abline(0,1): OLD VERSION
supp_fig2c <- tibble(r=V_inv %*% r) %>%
  ggplot(aes(sample=r)) +
  stat_qq(alpha=0.50) + 
  geom_abline(intercept = 0, slope = 1, color="red", lty=2) +
  ylab("Sample quantiles") +
  xlab("Theoretical quantiles") +
  theme(axis.title = element_text(face="bold"),
        axis.text  = element_text(face="bold"))

# -- Save supp-figure-2c
ggsave("figs/supp-figure-2c.pdf",
       plot   = supp_fig2c,
       dpi    = 300, 
       height = 6,
       width  = 8)

# -- Supp Figure 2D: No correlation in adjusted residuals
# acf(V_inv %*% r): OLD VERSION
auto_cor   <- acf(V_inv %*% r, plot = FALSE)
supp_fig2d <- tibble(acf = auto_cor$acf, lag = auto_cor$lag) %>%
  ggplot(aes(lag, acf)) +
  geom_col(color="black", fill="#252525", width = 0.5) +
  ylab("ACF") +
  xlab("Lag") +
  geom_hline(yintercept = qnorm((1+0.95)/2)/sqrt(auto_cor$n.used), 
             colour = "#cb181d",
             linetype = 2) +
  geom_hline(yintercept = -qnorm((1+0.95)/2)/sqrt(auto_cor$n.used), 
             colour = "#cb181d",
             linetype = 2) +
  theme(axis.title = element_text(face="bold"),
        axis.text  = element_text(face="bold"))

# -- Save supp-figure-2d
ggsave("figs/supp-figure-2d.pdf",
       plot   = supp_fig2d,
       dpi    = 300, 
       height = 6,
       width  = 8)

# -- Estimated sd
sd(V_inv %*% r)

# -- Supp-figure-2
supp_fig2 <- ggarrange(supp_fig2a, supp_fig2b, supp_fig2c, supp_fig2d, 
                       labels = c("A", "B", "C", "D"), 
                       nrow = 2, ncol = 2)

# -- Save supp-figure-2
ggsave("figs/supp-figure-2.pdf",
       plot   = supp_fig2,
       dpi    = 300, 
       height = 6,
       width  = 8)
### -- ------------------------------------- -----------------------------------------------------
### -- END Supp Figure 2 (Correlated errors) -----------------------------------------------------
### -- ------------------------------------- -----------------------------------------------------



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
  
  
  
  






