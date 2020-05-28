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

# -- Daily mortality data
counts <- puerto_rico_counts %>%
  group_by(date) %>%
  summarize(outcome    = sum(outcome),
            population = sum(population)) %>%
  ungroup()

# -- No displacement population
no_disp <- counts %>%
  filter(month(date) == 7, day(date) == 1) %>%
  mutate(year = year(date)) %>%
  select(year, population) %>%
  filter(year != 2018) %>%
  approx_demographics(demo=.,
                      first_day = ymd("1985-01-01"),
                      last_day = ymd("2020-12-31"))

# -- Number of knots per year
nknots <- 6
### -- ---------- ------------------------------------------------------------------
### -- END Set up ------------------------------------------------------------------
### -- ---------- ------------------------------------------------------------------



### -- ------------------------------------------- ------------------------------------------------------------------
### -- Supp Figure 8: Counterfactual consideration ------------------------------------------------------------------
### -- ------------------------------------------- ------------------------------------------------------------------
# -- Determining period of interest
start <- lubridate::make_date(2017, 01, 01)
end   <- lubridate::make_date(2018, 12, 31)

# -- Fitting the model with population displacement
fit <- excess_model(counts         = counts,
                    start          = start,
                    end            = end,
                    event          = hurricane_dates[3],
                    exclude        = exclude_dates,
                    control.dates  = control_dates,
                    knots.per.year = nknots,
                    model          = "correlated",
                    discontinuity  = TRUE)
disp <- tibble(date=fit$date, expected=fit$expected, fitted=fit$fitted, se=fit$se, flag="Displacement")
disp_ex <- excess_cumulative(fit, start=hurricane_dates[3], end=hurricane_dates[3]+180) %>%
  mutate(flag = "Displacement")

# -- Changing population to no disp
counts_no_disp <- counts %>%
  select(-population) %>%
  left_join(no_disp, by = "date")

# -- Fitting the model with no population displacement
fit2 <- excess_model(counts         = counts_no_disp,
                    start          = start,
                    end            = end,
                    event          = hurricane_dates[3],
                    exclude        = exclude_dates,
                    control.dates  = control_dates,
                    knots.per.year = nknots,
                    model          = "correlated",
                    discontinuity  = TRUE)
disp2 <- tibble(date=fit2$date, expected=fit2$expected, fitted=fit2$fitted, se=fit2$se, flag="No displacement")
disp2_ex <- excess_cumulative(fit2, start=hurricane_dates[3], end=hurricane_dates[3]+180) %>%
  mutate(flag = "No displacement")


supp_fig8a <- bind_rows(disp, disp2) %>%
  mutate(lwr = fitted-1.96*se,
         upr = fitted+1.96*se) %>%
  filter(date >= "2017-07-01", date <= "2018-07-01") %>%
  ggplot(aes(date, 100*fitted, color=flag, fill=flag)) +
  geom_hline(yintercept = 0, lty=2, color="gray") +
  geom_ribbon(aes(ymin=100*lwr, ymax=100*upr), alpha=0.50, color="transparent", show.legend = F) +
  geom_line(show.legend = FALSE) +
  geom_dl(aes(color=flag, label=flag),
          method=list("smart.grid")) +
  ylab("Percent increase from expected mortality") +
  xlab("") +
  scale_x_date(date_breaks = "2 month", 
               date_labels = "%b %Y") +
  scale_y_continuous(limits = c(-12, 90),
                     breaks = seq(-5, 85, by=10)) +
  theme(axis.text  = element_text(size=12),
        axis.title = element_text(size=13))

# -- Save supp figure 8A
ggsave("figs/supp-figure-8a.pdf",
       plot   = supp_fig8a,
       dpi    = 300, 
       height = 6,
       width  = 8)

supp_fig8b <- bind_rows(disp_ex, disp2_ex) %>%
  as_tibble() %>%
  mutate(lwr = fitted-1.96*se,
         upr = fitted+1.96*se,
         day = as.numeric(date - hurricane_dates[3])) %>%
  ggplot(aes(day, fitted, color=flag, fill=flag)) +
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.50, color="transparent", show.legend = F) +
  geom_line(show.legend = FALSE) +
  geom_dl(aes(color=flag, label=flag),
          method=list("smart.grid")) +
  ylab("Cumulative excess deaths") +
  xlab("Days after the event") +
  scale_x_continuous(limits = c(0, 180),
                     breaks = seq(0, 180, by=20)) +
  scale_y_continuous(limits = c(0, 3850),
                     labels = scales::comma,
                     breaks = seq(0, 3500, by=500)) +
  theme(axis.text  = element_text(size=12),
        axis.title = element_text(size=13))

# -- Save supp figure 8A
ggsave("figs/supp-figure-8b.pdf",
       plot   = supp_fig8b,
       dpi    = 300, 
       height = 6,
       width  = 8)
### -- ----------------------------------------------- ------------------------------------------------------------------
### -- END Supp Figure 8: Counterfactual consideration ------------------------------------------------------------------
### -- ----------------------------------------------- ------------------------------------------------------------------