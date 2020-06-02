# -- Set up
source("pr-init.R")

# -- Number of intervals
n <- 100

# -- Size of intervals in weeks
d <- c(1, 7, 14)

# -- Weekly mortality data for individuals 75 and over
counts <- filter(all_counts, agegroup == "75-Inf") %>%
  mutate(year = year(date),
         week = week(date)) %>%
  filter(week != 53) %>%
  group_by(year, week) %>%
  summarize(date       = date[1],
            outcome    = sum(outcome),
            population = mean(population)) %>%
  ungroup() %>%
  select(date, outcome, population)

# -- Dates to choose from
idx <- counts$date[between(counts$date, first(control_dates), last(control_dates))]

# -- Fitting models
set.seed(1) 
res <- map_df(d, function(i){
  
  # -- Randomly sampling intervals 
  intervals <- lapply(sample(idx, size=n), function(x){ seq(ymd(x), ymd(x) + i, by="days") })
  
  # -- Fitting poisson
  pois <- excess_model(counts         = counts,
                       intervals      = intervals, 
                       exclude        = exclude_dates,
                       control.dates  = control_dates,
                       model          = "poisson") %>%
    mutate(model = "poisson")
  
  # -- Fitting overdispersed poisson
  q_pois <- excess_model(counts         = counts,
                         intervals      = intervals, 
                         exclude        = exclude_dates,
                         control.dates  = control_dates,
                         model          = "quasipoisson") %>%
    mutate(model = "quasipoisson")
  
  # -- Fitting correlated
  correlated <- excess_model(counts         = counts,
                             intervals      = intervals, 
                             exclude        = exclude_dates,
                             control.dates  = control_dates,
                             model          = "correlated") %>%
    mutate(model = "correlated")
  
  
  bind_rows(pois, q_pois, correlated) %>%
    mutate(sample_size = i)
})

# -- Figure 5A
for(i in seq_along(d)){
  
  p = fig <- res %>%
    as_tibble() %>%
    mutate(model = case_when(model == "poisson" ~ "Poisson",
                             model == "quasipoisson" ~ "Overdispersed",
                             model == "correlated" ~ "Correlated"),
           model = factor(model, levels = c("Poisson", "Overdispersed", "Correlated"))) %>%
    filter(sample_size ==d[i]) %>%
    group_by(model) %>%
    mutate(a = excess / sd) %>%
    ungroup() %>%
    ggplot(aes(sample = a, color = model)) +
    geom_abline(lty=2) +
    stat_qq(alpha=0.50, size=3, show.legend = FALSE) +
    ylab("Sample quantiles") +
    xlab("Theoretical quantiles") +
    scale_y_continuous(limits = c(-8,8),
                       breaks = seq(-6, 6, by=3)) +
    scale_x_continuous(limits = c(-3,3),
                       breaks = seq(-3, 3, by=1)) +
    theme(axis.title = element_text(size=18),
          axis.text  = element_text(size=18))
  
  if(i == 2){
    p = fig <- res %>%
      as_tibble() %>%
      mutate(model = case_when(model == "poisson" ~ "Poisson",
                               model == "quasipoisson" ~ "Overdispersed",
                               model == "correlated" ~ "Correlated"),
             model = factor(model, levels = c("Poisson", "Overdispersed", "Correlated"))) %>%
      filter(sample_size ==d[i]) %>%
      group_by(model) %>%
      mutate(a = excess / sd) %>%
      ungroup() %>%
      ggplot(aes(sample = a, color = model)) +
      geom_abline(lty=2) +
      stat_qq(alpha=0.50, size=3, show.legend = T) +
      ylab("Sample quantiles") +
      xlab("Theoretical quantiles") +
      scale_y_continuous(limits = c(-8,8),
                         breaks = seq(-6, 6, by=3)) +
      scale_x_continuous(limits = c(-3,3),
                         breaks = seq(-3, 3, by=1)) +
      theme(axis.title  = element_text(size=18),
            axis.text   = element_text(size=18),
            legend.text = element_text(size=13),
            legend.title      = element_blank(),
            legend.direction  = "horizontal",
            legend.background = element_rect(color="black"),
            legend.position   = c(0.50, 0.90))
  }
  print(p)
  fn <- paste0("figs/supp-figure-9",letters[i], ".pdf")
  ggsave(fn, plot = fig, dpi  = 300, height = 5, width  = 7)
}
