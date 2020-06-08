### -- ------ ------------------------------------------------------------------
### -- Set up ------------------------------------------------------------------
### -- ------ ------------------------------------------------------------------
# -- Libraries
library(scales)
library(ggpubr)
library(tidyverse)
library(lubridate)
library(excessmort)
library(directlabels)
dslabs::ds_theme_set()

# -- Loading cook county data
data("cook_records")

# -- Wrangling mortality & population
the_breaks <- c(0, 5, 20, 40, 60, 75, Inf)
counts     <- compute_counts(cook_records, demo = cook_demographics, by = c("agegroup"), breaks = the_breaks) %>%
  mutate(agegroup = factor(agegroup, levels = c("0-4", "5-19", "20-39", "40-59", "60-74", "75-Inf")))

# -- Exclude & control dates
exclude_dates <- seq(ymd("2020-01-01"), max(counts$date), by="days")
control_dates <- seq(min(counts$date), ymd("2019-12-31"), by="days")
### -- ---------- ------------------------------------------------------------------
### -- END Set up ------------------------------------------------------------------
### -- ---------- ------------------------------------------------------------------

### -- ------------------------------------------------- ------------------------------------------------------------------
### -- Figure 3: Cook county Covid19 fhat white vs black ------------------------------------------------------------------
### -- ------------------------------------------------- ------------------------------------------------------------------
# -- Fitting models
cook <- map_df(unique(counts$agegroup), function(x){
  
  print(x)
  tmp_counts <- filter(counts, agegroup == x)
  f          <- excess_model(counts         = tmp_counts,
                             exclude        = exclude_dates,
                             control.dates = control_dates,
                             start          = ymd("2020-01-01"),
                             end            = max(counts$date),
                             weekday.effect = TRUE,
                             model          = "correlated",
                             verbose        = FALSE)
  
  tibble(date = f$date, obs = f$observed, mu = f$expected, fitted = f$fitted, se = f$se, agegroup = tmp_counts$agegroup[1])
}) %>%
  mutate(lwr=fitted-1.96*se,
         upr=fitted+1.96*se)

# -- Supp figure 
supp_fig6 <- cook %>%
  ggplot(aes(date, 100*fitted)) +
  geom_hline(yintercept = 0, lty=2, color="gray") +
  geom_ribbon(aes(ymin = 100*lwr, ymax=100*upr), alpha=0.50) + 
  geom_line() +
  ylab("Percent increase from expected mortality") +
  xlab("Date") +
  scale_x_date(date_breaks = "5 weeks", 
               date_labels = "%b %d") +
  scale_y_continuous(limits = c(-200, 2000),
                     breaks = seq(0, 2000, by=1000),
                     labels = scales::comma) +
  facet_wrap(~agegroup) +
  theme(axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15, angle=45, hjust=1),
        axis.title = element_text(size=17))

# -- Save figure 6
ggsave("figs/supp-figure-6.pdf",
       plot   = supp_fig6,
       dpi    = 300, 
       height = 5,
       width  = 7)
### -- ----------------------------------------------------- ------------------------------------------------------------------
### -- Figure 3: END Cook county Covid19 fhat white vs black ------------------------------------------------------------------
### -- ----------------------------------------------------- ------------------------------------------------------------------