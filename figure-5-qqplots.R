# -- Set up
source("pr-init.R")

# -- Number of intervals
n <- 100

# -- Size of intervals in days
d <- c(10, 50, 100)

# -- Dates to choose from
counts <- filter(all_counts, agegroup == "75-Inf")
idx    <- counts$date[between(counts$date, first(control_dates), last(control_dates))]

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
                       weekday.effect = TRUE,
                       model          = "poisson") %>%
    mutate(model = "poisson")
  
  # -- Fitting overdispersed poisson
  q_pois <- excess_model(counts         = counts,
                         intervals      = intervals, 
                         exclude        = exclude_dates,
                         control.dates  = control_dates,
                         weekday.effect = TRUE,
                         model          = "quasipoisson") %>%
    mutate(model = "quasipoisson")
  
  # -- Fitting correlated
  correlated <- excess_model(counts         = counts,
                             intervals      = intervals, 
                             exclude        = exclude_dates,
                             control.dates  = control_dates,
                             weekday.effect = TRUE,
                             model          = "correlated") %>%
    mutate(model = "correlated")
  
  
  bind_rows(pois, q_pois, correlated) %>%
    mutate(sample_size = i)
})

# -- Figure 5A
for(i in seq_along(d)){
  
  fig <- res %>%
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
    stat_qq(alpha=0.50, size=1, show.legend = FALSE) +
    ylab("Sample quantiles") +
    xlab("Theoretical quantiles") +
    scale_y_continuous(limits = c(-8,8),
                       breaks = seq(-6, 6, by=3)) +
    scale_x_continuous(limits = c(-3,3),
                       breaks = seq(-3, 3, by=1)) +
    theme(axis.title = element_text(size=18),
          axis.text  = element_text(size=18))
  
  if(i == 2){
    fig <- res %>%
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
      stat_qq(alpha=0.50, size=1, show.legend = T) +
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
  fn <- paste0("figs/figure-5",letters[i], ".pdf")
  ggsave(fn, plot = fig, dpi  = 300, height = 5, width  = 7)
}