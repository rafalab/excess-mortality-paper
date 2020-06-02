### -- --------------------------------------------------------- ------------------------------------------------------------------
### -- Figure 4: Mortality index: F-hat for bacterial infections ------------------------------------------------------------------
### -- --------------------------------------------------------- ------------------------------------------------------------------
# -- Set up
source("pr-init.R")

# -- Number of knots per year
nknots <- 6

# -- Causes of interest
icds   <- c("[A00,A79]", "[J00,J99]")
names  <- c("Bacterial infections", "Respiratory diseases")
causes <- tibble(icds, names)

# -- Fitting models
fit <- map_df(icds, function(x){
  
  # -- Subsetting and computing offset
  tmp_counts <- puerto_rico_icd %>%
                  group_by(date) %>%
                  mutate(population = sum(outcome)) %>%
                  ungroup() %>%
                  filter(icd == x)
  
  if(x == "[A00,A79]"){
    # -- Fitting excess model
    tmp_fit <- excess_model(tmp_counts, 
                            event = hurricane_dates[3],
                            start = make_date(2017,07,01),
                            end   = make_date(2018,07,01),
                            exclude        = exclude_dates,
                            control.dates  = control_dates,
                            knots.per.year = nknots,
                            weekday.effect = TRUE,
                            model          = "correlated",
                            discontinuity  = TRUE)
  } else {
    # -- Fitting excess model
    tmp_fit <- excess_model(tmp_counts,
                            start = make_date(2019,04,01),
                            end   = make_date(2020,04,01),
                            exclude        = exclude_dates,
                            control.dates  = control_dates,
                            # knots.per.year = nknots,
                            weekday.effect = TRUE,
                            model          = "correlated",
                            discontinuity  = TRUE)
  }
  tibble(date = tmp_fit$date, expected = tmp_fit$expected, observed = tmp_fit$observed  ,fitted = tmp_fit$fitted, se = tmp_fit$se, icd=x)
})

# -- Figure 4a
fig4a <- fit %>%
  filter(icd == "[A00,A79]") %>%
  mutate(lwr=fitted-1.96*se,
         upr=fitted+1.96*se) %>%
  left_join(causes, by=c("icd"="icds")) %>%
  ggplot(aes(date, 100*fitted)) +
  geom_hline(yintercept = 0, color="red3", lty=2) +
  geom_ribbon(aes(ymin=100*lwr, ymax=100*upr), alpha=0.50) +
  geom_line() +
  xlab("") +
  ylab("Percent increase from expected mortality") +
  scale_y_continuous(limits = c(-70, 140),
                     breaks = seq(-60, 140, by=30)) +
  scale_x_date(date_breaks = "3 month", 
               date_labels = "%b %Y") +
  theme(axis.text  = element_text(size=18),
        axis.title = element_text(size=17))

# -- Save figure 4a
ggsave("figs/figure-4a.pdf",
       plot   = fig4a,
       dpi    = 300, 
       height = 5,
       width  = 7)

# -- Figure 4b
fig4b <- fit %>%
  filter(icd == "[J00,J99]") %>%
  mutate(lwr=fitted-1.96*se,
         upr=fitted+1.96*se) %>%
  left_join(causes, by=c("icd"="icds")) %>%
  ggplot(aes(date, 100*fitted)) +
  geom_hline(yintercept = 0, color="red3", lty=2) +
  geom_ribbon(aes(ymin=100*lwr, ymax=100*upr), alpha=0.50) +
  geom_line() +
  xlab("") +
  ylab("Percent increase from expected mortality") +
  scale_y_continuous(limits = c(-70, 140),
                     breaks = seq(-60, 140, by=30)) +
  scale_x_date(date_breaks = "3 month", 
               date_labels = "%b %Y") +
  theme(axis.text  = element_text(size=18),
        axis.title = element_text(size=17))

# -- Save figure 4b
ggsave("figs/figure-4b.pdf",
       plot   = fig4b,
       dpi    = 300, 
       height = 5,
       width  = 7)
### -- ------------------------------------------------------------- ------------------------------------------------------------------
### -- END Figure 4: Mortality index: F-hat for bacterial infections ------------------------------------------------------------------
### -- ------------------------------------------------------------- ------------------------------------------------------------------