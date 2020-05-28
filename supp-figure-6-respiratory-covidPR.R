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
control_dates <- seq(as.Date("2002-01-01"), as.Date("2013-12-31"), by = "day")
exclude_dates <- c(seq(hurricane_dates[1], hurricane_effect_ends[1], by = "day"),
                   seq(hurricane_dates[2], hurricane_effect_ends[2], by = "day"),
                   seq(hurricane_dates[3], hurricane_effect_ends[3], by = "day"),
                   seq(as.Date("2014-09-01"), as.Date("2015-03-21"), by = "day"),
                   seq(as.Date("2001-01-01"), as.Date("2001-01-15"), by = "day"),
                   seq(as.Date("2020-01-01"), lubridate::today(), by = "day"))

# -- Loading data
data("puerto_rico_counts")

# -- Number of knots per year
nknots <- 6

# -- Causes of interest
icds   <- c("[J00,J99]")
names  <- c("Respiratory System")
causes <- tibble(icds, names)
### -- ---------- ------------------------------------------------------------------
### -- END Set up ------------------------------------------------------------------
### -- ---------- ------------------------------------------------------------------

### -- ------------------------------------------------------------ ------------------------------------------------------------------
### -- Supp Figure 6: Respiratory mortality index in PR for Covid19 ------------------------------------------------------------------
### -- ------------------------------------------------------------ ------------------------------------------------------------------
# -- Fitting models
fit <- map_df(icds, function(x){
  
  # -- Subsetting and computing offset
  tmp_counts <- puerto_rico_icd %>%
    group_by(date) %>%
    mutate(population = sum(outcome)) %>%
    ungroup() %>%
    filter(icd == x)
  
  # -- Fitting excess model
  tmp_fit <- excess_model(tmp_counts, 
                          start = make_date(2019,07,01),
                          end   = make_date(2020,05,01),
                          exclude        = exclude_dates,
                          control.dates  = control_dates,
                          knots.per.year = nknots,
                          weekday.effect = FALSE,
                          model          = "correlated",
                          discontinuity  = FALSE)
  
  tibble(date = tmp_fit$date, expected = tmp_fit$expected, observed = tmp_fit$observed  ,fitted = tmp_fit$fitted, se = tmp_fit$se, icd=x)
})

# -- Supp figure 6
supp_fig6 <- fit %>%
  mutate(lwr=fitted-1.96*se,
         upr=fitted+1.96*se) %>%
  left_join(causes, by=c("icd"="icds")) %>%
  ggplot(aes(date, 100*fitted)) +
  geom_hline(yintercept = 0, color="red3", lty=2) +
  geom_ribbon(aes(ymin=100*lwr, ymax=100*upr), alpha=0.50) +
  geom_line() +
  xlab("") +
  ylab("Percent increase from expected mortality") +
  scale_y_continuous(limits = c(-70, 40),
                     breaks = seq(-60, 40, by=10)) +
  scale_x_date(date_breaks = "2 month", 
               date_labels = "%b %Y") +
  theme(axis.text  = element_text(size=12),
        axis.title = element_text(size=13))

# -- Save figure 5
ggsave("figs/supp-figure-6.pdf",
       plot   = supp_fig6,
       dpi    = 300, 
       height = 6,
       width  = 8)
### -- ------------------------------------------------------------ ------------------------------------------------------------------
### -- END Supp Figure 6: Respiratory mortality index in PR for Covid19 ------------------------------------------------------------------
### -- ------------------------------------------------------------ ------------------------------------------------------------------
