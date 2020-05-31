### -- ------ ------------------------------------------------------------------
### -- Set up ------------------------------------------------------------------
### -- ------ ------------------------------------------------------------------
# -- Libraries
library(scales)
library(ggpubr)
library(tidyverse)
library(lubridate)
library(excessmort)
dslabs::ds_theme_set()

# -- Hurricanes information
hurricane_dates        <- as.Date(c("1989-09-18","1998-09-21","2017-09-20"))
hurricane_effect_ends  <- as.Date(c("1990-03-18","1999-03-21","2018-03-20"))
names(hurricane_dates) <- c("Hugo", "Georges", "Maria")

# -- Control & exclude periods
control_dates <- seq(as.Date("2006-01-01"), as.Date("2013-12-31"), by = "day")
exclude_dates <- c(seq(hurricane_dates[1], hurricane_effect_ends[1], by = "day"),
                   seq(hurricane_dates[2], hurricane_effect_ends[2], by = "day"),
                   seq(hurricane_dates[3], hurricane_effect_ends[3], by = "day"),
                   seq(as.Date("2004-09-01"), as.Date("2005-12-31"), by = "day"),
                   seq(as.Date("2014-09-01"), as.Date("2015-03-21"), by = "day"),
                   seq(as.Date("2001-01-01"), as.Date("2001-01-15"), by = "day"),
                   seq(as.Date("2020-01-01"), lubridate::today(), by = "day"))


# -- Loading data
data("puerto_rico_counts")
counts <- puerto_rico_counts %>%
  group_by(date) %>%
  summarize(outcome    = sum(outcome),
            population = sum(population)) %>%
  ungroup()

# -- All possible dates
dates <- counts$date

# -- Number of knots per year
nknots <- 4
### -- ---------- ------------------------------------------------------------------
### -- END Set up ------------------------------------------------------------------
### -- ---------- ------------------------------------------------------------------

### -- -------- ------------------------------------------------------------------
### -- Figure 5 ------------------------------------------------------------------
### -- -------- ------------------------------------------------------------------
# -- Number of intervals
set.seed(1)
n <- 100

# -- Size of intervals in days
d <- c(10, 50, 100)

# -- Dates to choose from
idx <- dates[between(dates, first(control_dates), last(control_dates))]

# -- Fitting models
res <- map_df(d, function(i){
  
  # -- Randomly sampling intervals 
  intervals <- lapply(sample(idx, size=n), function(x){ seq(ymd(x), ymd(x) + i, by="days") })
  
  # -- Fitting poisson
  pois <- excess_model(counts         = counts,
                       intervals      = intervals, 
                       exclude        = exclude_dates,
                       control.dates  = control_dates,
                       knots.per.year = nknots,
                       weekday.effect = TRUE,
                       model          = "poisson") %>%
    mutate(model = "poisson")
  
  # -- Fitting overdispersed poisson
  q_pois <- excess_model(counts         = counts,
                         intervals      = intervals, 
                         exclude        = exclude_dates,
                         control.dates  = control_dates,
                         knots.per.year = nknots,
                         weekday.effect = TRUE,
                         model          = "quasipoisson") %>%
    mutate(model = "quasipoisson")
  
  # -- Fitting correlated
  correlated <- excess_model(counts         = counts,
                             intervals      = intervals, 
                             exclude        = exclude_dates,
                             control.dates  = control_dates,
                             knots.per.year = nknots,
                             weekday.effect = TRUE,
                             model          = "correlated") %>%
    mutate(model = "correlated")
  
  
  bind_rows(pois, q_pois, correlated) %>%
    mutate(sample_size = i)
})

# -- Figure 5A
fig5a <- res %>%
  as_tibble() %>%
  mutate(model = case_when(model == "poisson" ~ "Poisson",
                           model == "quasipoisson" ~ "Overdispersed",
                           model == "correlated" ~ "Correlated"),
         model = factor(model, levels = c("Poisson", "Overdispersed", "Correlated"))) %>%
  filter(sample_size == 10) %>%
  group_by(model) %>%
  mutate(a = excess / sd) %>%
  ungroup() %>%
  ggplot(aes(sample = a, color = model)) +
  stat_qq(alpha=0.50, size=1, show.legend = FALSE) +
  geom_abline(lty=2) +
  ylab("Sample quantiles") +
  xlab("Theoretical quantiles") +
  scale_y_continuous(limits = c(-6.7,6.7),
                     breaks = seq(-6, 6, by=3)) +
  scale_x_continuous(limits = c(-3,3),
                     breaks = seq(-3, 3, by=1)) +
  theme(axis.title = element_text(size=13),
        axis.text  = element_text(size=13))


ggsave("figs/figure-5a.pdf",
       plot   = fig5a,
       dpi    = 300, 
       height = 4,
       width  = 6)

# -- Figure 5B
fig5b <- res %>%
  as_tibble() %>%
  mutate(model = case_when(model == "poisson" ~ "Poisson",
                           model == "quasipoisson" ~ "Overdispersed",
                           model == "correlated" ~ "Correlated"),
         model = factor(model, levels = c("Poisson", "Overdispersed", "Correlated"))) %>%
  filter(sample_size == 50) %>%
  group_by(model) %>%
  mutate(a = excess / sd) %>%
  ungroup() %>%
  ggplot(aes(sample = a, color = model)) +
  stat_qq(alpha=0.50, size=1) +
  geom_abline(lty=2) +
  ylab("Sample quantiles") +
  xlab("Theoretical quantiles") +
  scale_y_continuous(limits = c(-6.7,6.7),
                     breaks = seq(-6, 6, by=3)) +
  scale_x_continuous(limits = c(-3,3),
                     breaks = seq(-3, 3, by=1)) +
  theme(axis.title = element_text(size=13),
        axis.text  = element_text(size=13),
        legend.title = element_blank(),
        legend.background = element_rect(fill="white",
                                         color="black"),
        legend.position = c(0.80, 0.40))

ggsave("figs/figure-5b.pdf",
       plot   = fig5b,
       dpi    = 300, 
       height = 4,
       width  = 6)

# -- Figure 5C
fig5c <- res %>%
  as_tibble() %>%
  mutate(model = case_when(model == "poisson" ~ "Poisson",
                           model == "quasipoisson" ~ "Overdispersed",
                           model == "correlated" ~ "Correlated"),
         model = factor(model, levels = c("Poisson", "Overdispersed", "Correlated"))) %>%
  filter(sample_size == 100) %>%
  group_by(model) %>%
  mutate(a = excess / sd) %>%
  ungroup() %>%
  ggplot(aes(sample = a, color = model)) +
  stat_qq(alpha=0.50, size=1, show.legend = FALSE) +
  geom_abline(lty=2) +
  ylab("Sample quantiles") +
  xlab("Theoretical quantiles") +
  scale_y_continuous(limits = c(-6.7,6.7),
                     breaks = seq(-6, 6, by=3)) +
  scale_x_continuous(limits = c(-3,3),
                     breaks = seq(-3, 3, by=1)) +
  theme(axis.title = element_text(size=13),
        axis.text  = element_text(size=13))

ggsave("figs/figure-5c.pdf",
       plot   = fig5c,
       dpi    = 300, 
       height = 4,
       width  = 6)
### -- ------------ ------------------------------------------------------------------
### -- END Figure 5 ------------------------------------------------------------------
### -- ------------ ------------------------------------------------------------------