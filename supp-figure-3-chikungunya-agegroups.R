### -- ---------------------------------------- ------------------------------------------------------------------
### -- Supp Figure 3: Chikungunya by age groups ------------------------------------------------------------------
### -- ---------------------------------------- ------------------------------------------------------------------
# -- Set up
source("pr-init.R")

# -- Number of knots per year
nknots <- 6

# -- Determining period of interest
event  <- ymd("2014-07-14")
before <- days(365)
after  <- days(365)

# -- Fitting model only to the groups of interest
res <- map_df(levels(all_counts$agegroup), function(x){
  
  tmp_counts <- filter(all_counts, agegroup == x)
  tmp_fit    <- excess_model(counts         = tmp_counts, 
                             start          = event - before,
                             end            = event + after,
                             exclude        = exclude_dates,
                             control.dates  = control_dates,
                             knots.per.year = nknots,
                             weekday.effect = TRUE,
                             model          = "correlated",
                             discontinuity  = FALSE)
  
  tibble(date = tmp_fit$date, observed = tmp_fit$observed, expected = tmp_fit$expected,fitted = tmp_fit$fitted, se = tmp_fit$se, agegroup = x)
})

# -- Supplemental figure 3
supp_fig3 <- res %>%
  filter(date >= "2014-05-01", date <= "2015-05-01") %>%
  mutate(lwr = fitted-1.96*se,
         upr = fitted+1.96*se) %>%
  mutate(agegroup = factor(agegroup, levels = c("0-4", "5-19", "20-39", "40-59", "60-74", "75-Inf"))) %>%
  ggplot(aes(date, 100*fitted)) +
  geom_hline(yintercept = 0, lty=2, color="red") +
  geom_point(aes(date, 100*(observed/expected - 1)), alpha=0.50, size=0.20) +
  geom_ribbon(aes(ymin=100*lwr, ymax=100*upr), alpha=0.40, color="transparent") +
  geom_line() +
  xlab("") +
  ylab("Percent increase from expected mortality") +
  scale_x_date(date_breaks = "3 month", date_labels = "%b") +
  facet_wrap(~agegroup, scales="free_y") +
  theme(axis.text  = element_text(size=15),
        axis.title = element_text(size=17))

# -- Supplemental figure 3a
ggsave("figs/supp-figure-3.pdf",
       plot   = supp_fig3,
       dpi    = 300, 
       height = 5,
       width  = 7)
### -- -------------------------------------------- ------------------------------------------------------------------
### -- END Supp Figure 3: Chikungunya by age groups ------------------------------------------------------------------
### -- -------------------------------------------- ------------------------------------------------------------------