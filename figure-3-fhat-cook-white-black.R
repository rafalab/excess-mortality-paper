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
the_breaks <- c(seq(0, 85, by=5), Inf)
counts     <- compute_counts(cook_records, by=c("sex", "race", "agegroup"), breaks = the_breaks)
counts     <- left_join(counts, cook_demographics, by=c("date", "race", "sex", "agegroup"))

# -- Creating new agegroups and wrangle
the_breaks <- c(0, 5, 20, 40, 60, 75, Inf)
counts     <- collapse_counts_by_age(counts, the_breaks) %>%
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

# -- Number of knots
nknots <- 10
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
                             start          = ymd("2020-01-01"),
                             end            = max(counts$date),
                             weekday.effect = FALSE,
                             knots.per.year = nknots,
                             verbose        = FALSE)
  
  tibble(date = f$date, obs = f$observed, mu = f$expected, fitted = f$fitted, se = f$se, race = tmp_counts$race[1], agegroup = tmp_counts$agegroup[1])
})

# -- Figure 3
fig3 <- cook %>%
  mutate(agegroup = case_when(agegroup == "40-59" ~ "40 to 59 years",
                              agegroup == "60-74" ~ "60 to 74 years",
                              agegroup == "75-Inf" ~ "75 years and older"),
         race = str_to_sentence(race)) %>%
  mutate(lwr=fitted-1.96*se,
         upr=fitted+1.96*se) %>%
  ggplot(aes(date, 100*fitted, color=race, fill=race)) +
  geom_hline(yintercept = 0, lty=2, color="gray") +
  geom_ribbon(aes(ymin=100*lwr, ymax=100*upr), alpha=0.50, color="transparent", show.legend = F) +
  geom_line(size=1) +
  scale_color_manual(name="",
                     values = c("#D55E00","#0571b0","#009E73","#56B4E9","#CC79A7","#E69F00","#ca0020","gray")) +
  scale_fill_manual(name="",
                    values = c("#D55E00","#0571b0","#009E73","#56B4E9","#CC79A7","#E69F00","#ca0020","gray")) +
  ylab("Percent increase from expected mortality") +
  xlab("") +
  scale_x_date(date_breaks = "1 months", 
               date_labels = "%b %d") +
  scale_y_continuous(limits = c(-90, 2700),
                     breaks = seq(0, 2500, by=500)) +
  facet_wrap(~agegroup) +
  theme(axis.title = element_text(face="bold", color="black"),
        axis.text  = element_text(face="bold", color="black"),
        strip.text = element_text(face="bold", color="black"),
        legend.title    = element_blank(),
        legend.text     = element_text(face="bold", color="black", size=8),
        legend.background = element_rect(color    = "black",
                                         fill     = "white",
                                         linetype = "solid"),
        legend.position = c(0.50, 0.95),
        legend.key.size = unit(0.1, "cm"),
        legend.direction = "horizontal")

# -- Save figure 3
ggsave("figs/figure-3.pdf",
       plot   = fig3,
       dpi    = 300, 
       height = 4,
       width  = 8)
### -- ----------------------------------------------------- ------------------------------------------------------------------
### -- Figure 3: END Cook county Covid19 fhat white vs black ------------------------------------------------------------------
### -- ----------------------------------------------------- ------------------------------------------------------------------