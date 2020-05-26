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

# -- Age groups
the_breaks <- c(0, 5, 20, 40, 60, 75, Inf)
all_counts <- collapse_counts_by_age(puerto_rico_counts, the_breaks)
dates      <- unique(all_counts$date)

# -- Hurricanes information
hurricane_dates        <- as.Date(c("1989-09-18","1998-09-21","2017-09-20"))
hurricane_effect_ends  <- as.Date(c("1990-03-18","1999-03-21","2018-03-20"))
names(hurricane_dates) <- c("Hugo", "Georges", "Maria")

# -- Exclude periods
exclude_dates <- c(seq(hurricane_dates[1], hurricane_effect_ends[1], by = "day"),
                   seq(hurricane_dates[2], hurricane_effect_ends[2], by = "day"),
                   seq(hurricane_dates[3], hurricane_effect_ends[3], by = "day"),
                   seq(as.Date("2014-09-01"), as.Date("2015-03-21"), by = "day"),
                   seq(as.Date("2001-01-01"), as.Date("2001-01-15"), by = "day"),
                   seq(as.Date("2020-01-01"), lubridate::today(), by = "day"))
### -- ---------- ------------------------------------------------------------------
### -- END Set up ------------------------------------------------------------------
### -- ---------- ------------------------------------------------------------------

### -- ----------------------------------------- ------------------------------------------------------------------
### -- Supp Figure 1a: Trend effect by age group ------------------------------------------------------------------
### -- ----------------------------------------- ------------------------------------------------------------------
# -- Getting trend component
trend <- map_df(levels(all_counts$agegroup), function(x){
  
  # -- Subset by age group
  counts <- filter(all_counts, agegroup == x)
  
  # -- Fit mean model
  fit <- compute_expected(counts, 
                          exclude         = exclude_dates,
                          weekday.effect  = FALSE,
                          keep.components = TRUE, 
                          verbose         = FALSE)
  
  # -- Extracting trend component
  tibble(date = dates, trend = fit$trend, agegroup = x)
})

# -- Supplemental figure 1a
supp_fig1a <- trend %>%
  mutate(agegroup = factor(agegroup, levels = c("0-4", "5-19", "20-39", "40-59", "60-74", "75-Inf"))) %>%
  ggplot(aes(date, trend)) +
  geom_line() +
  ylab("Estimated trend effect") +
  xlab("Date") +
  facet_wrap(~agegroup, scales="free_y")

# -- Save supplemental figure 1
ggsave("figs/supp-figure-1a.pdf",
       plot   = supp_fig1a,
       dpi    = 300, 
       height = 4,
       width  = 8)
### -- --------------------------------------------- ------------------------------------------------------------------
### -- END Supp Figure 1a: Trend effect by age group ------------------------------------------------------------------
### -- --------------------------------------------- ------------------------------------------------------------------

### -- -------------------------------------------- ------------------------------------------------------------------
### -- Supp Figure 2a: Seasonal effect by age group ------------------------------------------------------------------
### -- -------------------------------------------- ------------------------------------------------------------------
# -- Getting seasonal component
seasonal <- map_df(levels(all_counts$agegroup), function(x){
  
  # -- Subset by age group
  counts <- filter(all_counts, agegroup == x)
  
  # -- Fit mean model
  fit <- compute_expected(counts, 
                          exclude         = exclude_dates,
                          weekday.effect  = FALSE,
                          keep.components = TRUE, 
                          verbose         = FALSE)
  
  # -- Extracting seasonal component
  tibble(fit$s) %>%
    mutate(agegroup = x)
})

# -- Supplemental figure 1b
supp_fig1b <- seasonal %>%
  mutate(agegroup = factor(agegroup, levels = c("0-4", "5-19", "20-39", "40-59", "60-74", "75-Inf"))) %>%
  ggplot(aes(day, s)) +
  geom_line() +
  ylab("Estimated seasonal effect") +
  xlab("Day of the year") +
  facet_wrap(~agegroup, scales="free_y")

# -- Save supplemental figure 2
ggsave("figs/supp-figure-1b.pdf",
       plot   = supp_fig1b,
       dpi    = 300, 
       height = 4,
       width  = 8)
### -- ------------------------------------------------ ------------------------------------------------------------------
### -- END Supp Figure 2a: Seasonal effect by age group ------------------------------------------------------------------
### -- ------------------------------------------------ ------------------------------------------------------------------
