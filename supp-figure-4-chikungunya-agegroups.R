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

# -- Control & exclude periods
control_dates <- seq(as.Date("2002-01-01"), as.Date("2013-12-31"), by = "day")
exclude_dates  <- c(seq(puerto_rico_hurricane_dates[1], puerto_rico_hurricane_effect_ends[1], by = "day"),
                    seq(puerto_rico_hurricane_dates[2], puerto_rico_hurricane_effect_ends[2], by = "day"),
                    seq(puerto_rico_hurricane_dates[3], puerto_rico_hurricane_effect_ends[3], by = "day"),
                    seq(as.Date("2014-01-01"), as.Date("2015-12-31"), by = "day"),
                    seq(as.Date("2001-01-01"), as.Date("2001-01-15"), by = "day"),
                    seq(as.Date("2020-01-01"), lubridate::today(), by = "day"))

# -- Age groups
# the_breaks <- c(0, 5, 20, 40, 60, Inf)
the_breaks <- c(0, 5, 20, 40, 60, 75, Inf)
all_counts <- collapse_counts_by_age(puerto_rico_counts, the_breaks)

# -- Number of knots per year
nknots <- 6
### -- ---------- ------------------------------------------------------------------
### -- END Set up ------------------------------------------------------------------
### -- ---------- ------------------------------------------------------------------

### -- ---------------------------------------- ------------------------------------------------------------------
### -- Supp Figure 4: Chikungunya by age groups ------------------------------------------------------------------
### -- ---------------------------------------- ------------------------------------------------------------------
# -- Determining period of interest
start <- lubridate::make_date(2014, 01, 01)
end   <- lubridate::make_date(2015, 12, 31)

# -- Fitting model only to the groups of interest
res <- map_df(levels(all_counts$agegroup), function(x){
  
  tmp_counts <- filter(all_counts, agegroup == x)
  tmp_fit    <- excess_model(counts         = tmp_counts, 
                             start          = start,
                             end            = end,
                             exclude        = exclude_dates,
                             control.dates  = control_dates,
                             knots.per.year = nknots,
                             model          = "correlated",
                             discontinuity  = FALSE)
  
  tibble(date = tmp_fit$date, fhat = tmp_fit$fitted, se = tmp_fit$se, agegroup = x)
})

# -- Supplemental figure 4
supp_fig4 <- res %>%
  filter(date >= "2014-05-01", date <= "2015-05-01") %>%
  mutate(lwr = fhat-1.96*se,
         upr = fhat+1.96*se) %>%
  mutate(agegroup = factor(agegroup, levels = c("0-4", "5-19", "20-39", "40-59", "60-74", "75-Inf"))) %>%
  ggplot(aes(date, 100*fhat)) +
  geom_hline(yintercept = 0, lty=2, color="gray") +
  geom_ribbon(aes(ymin=100*lwr, ymax=100*upr), alpha=0.40, color="transparent") +
  geom_line(size=1) +
  xlab("") +
  ylab("Percent change in mortality") +
  scale_x_date(date_breaks = "3 month", date_labels = "%b %y") +
  facet_wrap(~agegroup) +
  theme(axis.title = element_text(face="bold", color="black"),
        axis.text  = element_text(face="bold", color="black"),
        strip.text = element_text(face="bold", color="black"))

# -- Save # -- Supplemental figure 4
ggsave("figs/supp-figure-4.pdf",
       plot   = supp_fig4,
       dpi    = 300, 
       height = 4,
       width  = 8)
### -- -------------------------------------------- ------------------------------------------------------------------
### -- END Supp Figure 4: Chikungunya by age groups ------------------------------------------------------------------
### -- -------------------------------------------- ------------------------------------------------------------------