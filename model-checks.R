# !diagnostics off
library(tidyverse)
dslabs::ds_theme_set()
library(excessdeaths)
library(lubridate)
data("new_jersey_counts")
# Simulation study --------------------------------------------------------
counts <- new_jersey_counts
control_dates <-  seq(make_date(2007,1,1), make_date(2011, 10, 1), by = "day")
exclude_dates <- make_date(2012,0,29)+0:180
counts <- compute_expected(counts, exclude = exclude_dates)

### Simluation
dates <- counts$date
year <- 2012; el <- 30; M <- .5 ## el = event size
event_day<- lubridate::make_date(year, 9, 20)  
out_dates <- seq(event_day - 3.5*el, event_day + 3.5*el, by = "day")
f <- dnorm(as.numeric(dates - event_day), 0, el); f <- f*M/max(f)

## fit ar model
B <- 25
fits <- lapply(1:B, function(i){
  print(i)
  mu <- counts$expected
  sim <- counts
  
  e <- pmax(0, 1 + arima.sim(list(order = c(5,0,0), ar = seq(0.1, 0.01, len = 5)), sd = 0.05, n = length(dates)))
  
  sim$outcome <- rpois(length(dates), mu * (1+f) * e)
  
  tmp <- excess_model(sim, 
                     start = lubridate::make_date(year-1, 9, 20),
                     end =  lubridate::make_date(year+1, 9, 20),
                     exclude = exclude_dates,
                     control.dates = control_dates,
                     nknots = 24,
                     disc = FALSE, 
                     model  = "correlated")
  
  return(tmp)
})
rafalib::mypar(1,1)
index <- match(fits[[1]]$date, dates)
plot(dates[index], f[index], type="l") 
sapply(fits, function(f) lines(f$date, f$fitted, col=2))
lines(dates[index], f[index], type="l") 

##check fit
m <- sapply(fits, function(f) f$fitted)
s <- sapply(fits, function(f) f$se)

##Bias
plot(dates[index], rowMeans(m))
lines(dates[index], f[index], col = 2)
## Actual SE veruese estimated SEs
plot(dates[index], matrixStats::rowSds(m), ylim = range(c(matrixStats::rowSds(m)), range(rowMeans(s))))
matlines(dates[index], s)

## excess deaths

cumu <- lapply(fits, function(fit){
  excess_cumulative(fit, min(dates), max(dates)) 
})

true_f <- cumsum(counts$expected[index]*f[index])
observed <- sapply(cumu, function(r) r$observed)
fitted <- sapply(cumu, function(r) r$fitted)
fitted_se <- sapply(cumu, function(r) r$fitted_se)
se <- sapply(cumu, function(r) r$se)

mf <- rowMeans(fitted)
plot(true_f)
lines(mf, col = 2)

true_s <- matrixStats::rowSds(fitted)
est_s <- rowMeans(fitted_se)
plot(true_s)
matlines(fitted_se, col="grey", type="l", lty=1)
lines(est_s, col = 2)

true_s <- matrixStats::rowSds(observed)
est_s <- rowMeans(se)
plot(true_s)
matlines(se, col="grey", type="l")
lines(est_s, col = 2)



# Now check naive model ---------------------------------------------------

B <- 100
fits <- lapply(1:B, function(i){
  print(i)
  mu <- counts$expected
  sim <- counts
  
  e <- pmax(0, 1 - rnorm(length(dates), 0, 0.05))

  sim$outcome <- rpois(length(dates), mu * (1+f) * e)
  
  tmp <- excess_model(sim, 
                         start = lubridate::make_date(year-1, 9, 20),
                         end =  lubridate::make_date(year+1, 9, 20),
                         exclude = out_dates,
                         nknots = 24,
                         disc=FALSE, 
                      model = "quasipoisson")
  
  return(tmp)
})
rafalib::mypar(1,1)
index <- match(fits[[1]]$date, dates)
plot(dates[index], f[index], type="l") 
sapply(fits, function(f) lines(f$date, f$fitted, col=2))
lines(dates[index], f[index], type="l") 

##check fit
m <- sapply(fits, function(f) f$fitted)
s <- sapply(fits, function(f) f$se)

##Bias
plot(dates[index], rowMeans(m))
lines(dates[index], f[index], col = 2)
## Actual SE veruese estimated SEs
plot(dates[index], matrixStats::rowSds(m), ylim = range(c(matrixStats::rowSds(m)), range(rowMeans(s))))
matlines(dates[index], s)


## excess deaths

cumu <- lapply(fits, function(fit){
  excess_cumulative(fit, min(dates), max(dates)) 
})

true_f <- cumsum(res$expected[index]*f[index])
observed <- sapply(cumu, function(r) r$observed)
fitted <- sapply(cumu, function(r) r$fitted)
fitted_se <- sapply(cumu, function(r) r$fitted_se)
se <- sapply(cumu, function(r) r$se)

mf <- rowMeans(fitted)
plot(true_f)
lines(mf, col = 2)

true_s <- matrixStats::rowSds(fitted)
est_s <- rowMeans(fitted_se)
plot(true_s)
matlines(fitted_se, col="grey", type="l", lty=1)
lines(est_s, col = 2)


true_s <- matrixStats::rowSds(observed)
est_s <- rowMeans(se)
plot(true_s)
matlines(se, col="grey", type="l")
lines(est_s, col = 2)

# Diagnostics: empirical check for normality with correct mean and sd ---------------

intervals <- lapply(sample(control_dates, 100), function(s) seq(s, s+14, by="day"))
tab1 <- purrr::map_df(levels(all_counts$agegroup), function(g){
  counts <- dplyr::filter(all_counts, agegroup == g)
  tab <- excess_model(counts, 
                      intervals = intervals, 
                      exclude = out_dates,
                      model = "correlated",  
                      control.dates = control_dates)$excess
  dplyr::mutate(tab, agegroup = g)
})

for(g in levels(all_counts$agegroup)){
  with(dplyr::filter(tab1, agegroup == g), {
    z= excess/excess_se; 
    qqnorm(z, main = g);
    abline(0,1)})
}



