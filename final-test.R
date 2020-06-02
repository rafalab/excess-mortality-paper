library(tidyverse)
library(lubridate)
library(excessmort)
dslabs::ds_theme_set()

# -- Hurricanes information
hurricane_dates  <- c(Hugo = make_date(1989, 9, 118),
                      Georges = make_date(1998, 9, 21),
                      Maria = make_date(2017, 9, 20))


# Excluding  1) years (starting in 7/1) that include events, 2) 2020 and 3) some outliers in 2001 -------------------

exclude_dates <- c(make_date(1989, 7, 1) + 0:364,
                   make_date(1998, 7, 1) + 0:365,
                   make_date(2017, 7, 1) + 0:365,
                   make_date(2014, 7, 1) + 0:365,
                   seq(make_date(2001, 1, 1), make_date(2001, 1, 15), by = "day"),
                   seq(make_date(2020, 1, 1), today(), by = "day"))
  
control_dates <- seq(as.Date("2002-01-01"), as.Date("2013-12-31"), by = "day")


data("puerto_rico_counts")

##Look at the all the fits

the_breaks <- c(75, Inf)
counts <- collapse_counts_by_age(puerto_rico_counts, the_breaks)



starts <- make_date(seq(1985, 2020, 3), 7, 1)
ends <- pmin(starts + years(3), max(counts$date))

e <- compute_expected(counts, exclude_dates)

 
for(i in seq_along(starts)){
  
  p1 <- expected_plot(e,  start = starts[i],
                      end = ends[i], title = starts[i])
  
  fit <- excess_model(counts, 
                      start = starts[i],
                      end = ends[i],
                      exclude = exclude_dates,
                      control.dates = control_dates,
                      model = "correlated",
                      weekday.effect = TRUE)
  
  p2 <- excess_plot(fit, show.data = FALSE)
  
  gridExtra::grid.arrange(p1, p2, nrow = 2)
}
                    
                    

### combined 
the_breaks <- c(0, 5, 20, 40, 60, 75, Inf)
all_counts <- collapse_counts_by_age(puerto_rico_counts, the_breaks)



starts <- make_date(seq(1985, 2020, 3), 7, 1)
ends <- pmin(starts + years(3), max(counts$date))


for(i in seq_along(starts)){
  tmp <- map_df(unique(all_counts$agegroup), function(x){
    
    f <- all_counts %>% 
      filter(agegroup == x) %>%
      excess_model(start = starts[i],
                   end = ends[i],
                   exclude = exclude_dates,
                   control.dates = control_dates,
                   model = "correlated",
                   weekday.effect = TRUE)
    tibble(date = f$date, expected = f$expected, fitted = f$fitted, se = f$se, agegroup = x)
  })
  p <- tmp %>% 
    group_by(date) %>% 
    summarize(fitted = sum(expected * fitted) / sum(expected), 
              se     = sqrt(sum(expected^2 * se^2)) / sum(expected)) %>%
    ungroup() %>%
    ggplot(aes(date, fitted)) + 
    geom_ribbon(aes(ymin = fitted - 2*se, ymax = fitted + 2*se), alpha = 0.25)+
    geom_line() +
    geom_hline(yintercept = 0, lty = 2)
    print(p)
}

      
  
  




