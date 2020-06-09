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

# -- Creating age groups and wrangle
the_breaks <- c(0, 5, 20, 40, 60, 75, Inf)
counts <- compute_counts(dat = cook_records, demo = cook_demographics, by=c("sex", "race", "agegroup"), breaks = the_breaks)

# -- Creating new agegroups and wrangle
counts <- counts %>%
  filter(race %in% c("white", "black")) %>%
  group_by(date, race, agegroup) %>%
  summarize(outcome    = sum(outcome, na.rm = TRUE),
            population = sum(population, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(agegroup = factor(agegroup, levels = c("0-4", "5-19", "20-39", "40-59", "60-74", "75-Inf"))) %>%
  mutate(group = paste0(race," | ",agegroup)) %>%
  filter(!agegroup %in% c("0-4", "5-19", "20-39"))

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
cook <- map_df(unique(counts$group), function(x){
  
  print(x)
  tmp_counts <- filter(counts, group == x)
  f          <- excess_model(counts         = tmp_counts,
                             exclude        = exclude_dates,
                             control.dates = control_dates,
                             start          = ymd("2020-01-01"),
                             end            = max(counts$date),
                             weekday.effect = TRUE,
                             model          = "correlated",
                             verbose        = FALSE)
  
  tibble(date = f$date, obs = f$observed, mu = f$expected, fitted = f$fitted, se = f$se, race = tmp_counts$race[1], agegroup = tmp_counts$agegroup[1])
}) %>%
  mutate(group = case_when(agegroup == "40-59" ~ "40 to 59 years",
                          agegroup == "60-74" ~ "60 to 74 years",
                          agegroup == "75-Inf" ~ "75 years and older"),
         race = str_to_sentence(race)) %>%
  mutate(lwr=fitted-1.96*se,
         upr=fitted+1.96*se)

# -- Figure 3a
fig3a <- cook %>%
  filter(agegroup == "40-59") %>%
  ggplot(aes(date, 100*fitted, color=race, fill=race)) +
  geom_hline(yintercept = 0, lty=2, color="gray") +
  geom_ribbon(aes(ymin=100*lwr, ymax=100*upr), alpha=0.50, color="transparent", show.legend = F) +
  geom_line(show.legend = F) +
  ylab("Percent increase from expected mortality") +
  xlab("Date") +
  scale_x_date(date_breaks = "1 months", 
               date_labels = "%b %d") +
  scale_y_continuous(limits = c(-115, 2700),
                     breaks = seq(0, 2500, by=500),
                     labels = scales::comma) +
  theme(legend.background = element_rect(color    = "black",
                                         fill     = "white",
                                         linetype = "solid"),
        legend.title     = element_blank(),
        legend.position  = c(0.50, 0.95),
        legend.direction = "horizontal",
        legend.text      = element_text(size=13),
        axis.text  = element_text(size=18),
        axis.title = element_text(size=17))

# -- Save figure 3a
ggsave("figs/figure-3a-cook-40-59.pdf",
       plot   = fig3a,
       dpi    = 300, 
       height = 5,
       width  = 7)

# -- Figure 3b
fig3b <- cook %>%
  filter(agegroup == "60-74") %>%
  ggplot(aes(date, 100*fitted, color=race, fill=race)) +
  geom_hline(yintercept = 0, lty=2, color="gray") +
  geom_ribbon(aes(ymin=100*lwr, ymax=100*upr), alpha=0.50, color="transparent", show.legend = F) +
  geom_line() +
  ylab("Percent increase from expected mortality") +
  xlab("Date") +
  scale_x_date(date_breaks = "1 months", 
               date_labels = "%b %d") +
  scale_y_continuous(limits = c(-115, 2700),
                     breaks = seq(0, 2500, by=500),
                     labels = scales::comma) +
  theme(legend.background = element_rect(color    = "black",
                                         fill     = "white",
                                         linetype = "solid"),
        legend.title     = element_blank(),
        legend.text      = element_text(size=13),
        legend.position  = c(0.50, 0.90),
        legend.direction = "horizontal",
        axis.text  = element_text(size=18),
        axis.title = element_text(size=17))

# -- Save figure 3b
ggsave("figs/figure-3b-cook-60-74.pdf",
       plot   = fig3b,
       dpi    = 300, 
       height = 5,
       width  = 7)

# -- Figure 3c
fig3c <- cook %>%
  filter(agegroup == "75-Inf") %>%
  ggplot(aes(date, 100*fitted, color=race, fill=race)) +
  geom_hline(yintercept = 0, lty=2, color="gray") +
  geom_ribbon(aes(ymin=100*lwr, ymax=100*upr), alpha=0.50, color="transparent", show.legend = F) +
  geom_line(show.legend = F) +
  ylab("Percent increase from expected mortality") +
  xlab("Date") +
  scale_x_date(date_breaks = "1 months", 
               date_labels = "%b %d") +
  scale_y_continuous(limits = c(-115, 2700),
                     breaks = seq(0, 2500, by=500),
                     labels = scales::comma) +
  theme(legend.background = element_rect(color    = "black",
                                         fill     = "white",
                                         linetype = "solid"),
        legend.title     = element_blank(),
        legend.position  = c(0.50, 0.90),
        legend.direction = "horizontal",
        legend.text      = element_text(size=13),
        axis.text  = element_text(size=18),
        axis.title = element_text(size=17))

# -- Save figure 3b
ggsave("figs/figure-3c-cook-75-above.pdf",
       plot   = fig3c,
       dpi    = 300, 
       height = 5,
       width  = 7)
### -- ----------------------------------------------------- ------------------------------------------------------------------
### -- Figure 3: END Cook county Covid19 fhat white vs black ------------------------------------------------------------------
### -- ----------------------------------------------------- ------------------------------------------------------------------