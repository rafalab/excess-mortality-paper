### -- ----------------------------------------- ------------------------------------------------------------------
### -- Supp Figure 1a: Trend effect by age group ------------------------------------------------------------------
### -- ----------------------------------------- ------------------------------------------------------------------

# -- Set up
source("pr-init.R")
# dates<- unique(all_counts$date)

# -- Getting trend component
trend <- map_df(levels(all_counts$agegroup), function(x){
  
  # -- Subset by age group
  counts <- filter(all_counts, agegroup == x)
  
  # -- Fit mean model
  fit <- compute_expected(counts, 
                          exclude         = exclude_dates,
                          weekday.effect  = TRUE,
                          keep.components = TRUE, 
                          verbose         = FALSE)
  
  
  # -- Extracting trend component
  tibble(date = fit$counts$date, trend = fit$trend, agegroup = x)
})

# -- Supplemental figure 1a
supp_fig1a <- trend %>%
  mutate(agegroup = factor(agegroup, levels = c("0-4", "5-19", "20-39", "40-59", "60-74", "75-Inf"))) %>%
  ggplot(aes(date, trend)) +
  geom_line() +
  ylab("Estimated trend effect") +
  xlab("") +
  facet_wrap(~agegroup, scales="free_y") +
  # facet_wrap(~agegroup) +
  theme(axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15, angle=45, hjust=1),
        axis.title = element_text(size=17))
  
# -- Save supplemental figure 1
ggsave("figs/supp-figure-1a.pdf",
       plot   = supp_fig1a,
       dpi    = 300, 
       height = 5,
       width  = 7)
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
                          weekday.effect  = TRUE,
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
  facet_wrap(~agegroup) +
  theme(axis.text  = element_text(size=15),
        axis.title = element_text(size=17))

# -- Save supplemental figure 2
ggsave("figs/supp-figure-1b.pdf",
       plot   = supp_fig1b,
       dpi    = 300, 
       height = 5,
       width  = 7)
### -- ------------------------------------------------ ------------------------------------------------------------------
### -- END Supp Figure 2a: Seasonal effect by age group ------------------------------------------------------------------
### -- ------------------------------------------------ ------------------------------------------------------------------
