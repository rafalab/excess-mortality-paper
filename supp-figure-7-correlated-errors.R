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

# -- Age groups
the_breaks <- c(0, 5, 20, 40, 60, 75, Inf)
all_counts <- collapse_counts_by_age(puerto_rico_counts, the_breaks)
dates      <- unique(all_counts$date)

# -- Hurricanes information
hurricane_dates        <- as.Date(c("1989-09-18","1998-09-21","2017-09-20"))
hurricane_effect_ends  <- as.Date(c("1990-03-18","1999-03-21","2018-03-20"))
names(hurricane_dates) <- c("Hugo", "Georges", "Maria")

# -- Control & exclude periods
control_dates <- seq(as.Date("2006-01-01"), as.Date("2013-12-31"), by = "day")
exclude_dates <- c(seq(hurricane_dates[1], hurricane_effect_ends[1], by = "day"),
                   seq(hurricane_dates[2], hurricane_effect_ends[2], by = "day"),
                   seq(hurricane_dates[3], hurricane_effect_ends[3], by = "day"),
                   seq(as.Date("2004-09-01"), as.Date("2005-12-31"), by = "day"),
                   seq(as.Date("2014-09-01"), as.Date("2015-03-21"), by = "day"),
                   seq(as.Date("2001-01-01"), as.Date("2001-01-15"), by = "day"),
                   seq(as.Date("2020-01-01"), lubridate::today(), by = "day"))
### -- ---------- ------------------------------------------------------------------
### -- END Set up ------------------------------------------------------------------
### -- ---------- ------------------------------------------------------------------

### -- --------------------------------------------------------- -----------------------------------------------------
### -- Supp Figure 7a: QQ-plot without adjusting for correlation -----------------------------------------------------
### -- --------------------------------------------------------- -----------------------------------------------------
# -- Daily mortality data
counts <- puerto_rico_counts %>%
            group_by(date) %>%
            summarize(outcome    = sum(outcome),
                      population = sum(population)) %>%
            ungroup()

# -- Computing expected mortality counts
counts <- compute_expected(counts, exclude = exclude_dates, weekday.effect = TRUE)

# -- Dates for model check
example_dates <- control_dates[year(control_dates) >= 2005]
first(example_dates)
last(example_dates)

# -- Computing z scores for example dates (H0 true)
r <- tibble(date = counts$date, observed = counts$outcome, expected = counts$expected) %>%
        filter(date %in% example_dates) %>%
        mutate(r = (observed - expected)/sqrt(expected)) %>%
        pull(r)

# -- Supp Figure 7A: Poisson model doesn't fit the tails
supp_fig7a <- tibble(r=r) %>%
  ggplot(aes(sample=r)) +
  stat_qq(alpha=0.50) + 
  geom_abline(intercept = 0, slope = 1, color="red", lty=2) +
  xlab("Sample quantiles") +
  ylab("Theoretical quantiles") +
  theme(axis.text  = element_text(size=18),
        axis.title = element_text(size=18))

# -- Save supp-figure-7a
ggsave("figs/supp-figure-7a.pdf",
       plot   = supp_fig7a,
       dpi    = 300, 
       height = 5,
       width  = 7)
### -- ------------------------------------------------------------- -----------------------------------------------------
### -- END Supp Figure 7a: QQ-plot without adjusting for correlation -----------------------------------------------------
### -- ------------------------------------------------------------- -----------------------------------------------------

### -- -------------------------------------------- -----------------------------------------------------
### -- Supp Figure 7b: Show correlation in the data -----------------------------------------------------
### -- -------------------------------------------- -----------------------------------------------------
# -- Supp Figure 7b
auto_cor <- acf(r, plot=FALSE)
supp_fig7b <- tibble(acf = auto_cor$acf, lag = auto_cor$lag) %>%
  ggplot(aes(lag, acf)) +
  geom_col(color="black", fill="#252525", width = 0.5) +
  ylab("ACF") +
  xlab("Lag") +
  geom_hline(yintercept = qnorm((1+0.95)/2)/sqrt(auto_cor$n.used), 
             colour = "#cb181d",
             linetype = 2) +
  geom_hline(yintercept = -qnorm((1+0.95)/2)/sqrt(auto_cor$n.used), 
             colour = "#cb181d",
             linetype = 2) +
  theme(axis.text  = element_text(size=18),
        axis.title = element_text(size=18))

# -- Save supp-figure-7b
ggsave("figs/supp-figure-7b.pdf",
       plot   = supp_fig7b,
       dpi    = 300, 
       height = 5,
       width  = 7)

# -- The SD is not 1
sd(r)
### -- ------------------------------------------------ -----------------------------------------------------
### -- END Supp Figure 7b: Show correlation in the data -----------------------------------------------------
### -- ------------------------------------------------ -----------------------------------------------------

### -- ------------------------------------------------- -----------------------------------------------------
### -- Supp Figure 7c: QQ-plot adjusting for correlation -----------------------------------------------------
### -- ------------------------------------------------- -----------------------------------------------------
# -- Now show that we if adjust for correlation things get better
tmp <- counts %>%
  filter(date %in% example_dates) %>%
  mutate(r = (outcome - expected)/expected)

# -- Extracting percent change
r <- tmp$r

# -- Expected value
mu <- tmp$expected

# -- Number of data points
n <- length(r)

# -- Estimating AR(p) with control dates
arfit <- excessmort:::fit_ar(tmp,  control.dates = control_dates)

# -- Estimated rhos
rhos <- ARMAacf(ar = arfit$ar, ma = 0, lag.max = n)

# -- Estimated sd
s <- arfit$sigma

# -- Computing covariance matrix
Sigma <- apply(abs(outer(1:n, 1:n, "-")) + 1, 1, function(i) rhos[i]) * outer(sqrt(s^2 + 1/mu), sqrt(s^2 + 1/mu))

# -- Cholesky decomposition of Sigma & computing inverse
V     <- chol(Sigma)
V_inv <- backsolve(r = V, x = diag(ncol(V))) ## V is upper triangular so backsolve faster

# -- Supp Figure 7C
supp_fig7c <- tibble(r=V_inv %*% r) %>%
  ggplot(aes(sample=r)) +
  stat_qq(alpha=0.50) + 
  geom_abline(intercept = 0, slope = 1, color="red", lty=2) +
  xlab("Sample quantiles") +
  ylab("Theoretical quantiles") +
  theme(axis.text  = element_text(size=18),
        axis.title = element_text(size=18))

# -- Save supp-figure-7c
ggsave("figs/supp-figure-7c.pdf",
       plot   = supp_fig7c,
       dpi    = 300, 
       height = 5,
       width  = 7)
### -- ----------------------------------------------------- -----------------------------------------------------
### -- END Supp Figure 7c: QQ-plot adjusting for correlation -----------------------------------------------------
### -- ----------------------------------------------------- -----------------------------------------------------

### -- ---------------------------------------------------- -----------------------------------------------------
### -- Supp Figure 7d: No correlation in adjusted residuals -----------------------------------------------------
### -- ---------------------------------------------------- -----------------------------------------------------
# -- Supp Figure 7d
auto_cor   <- acf(V_inv %*% r, plot = FALSE)
supp_fig7d <- tibble(acf = auto_cor$acf, lag = auto_cor$lag) %>%
  ggplot(aes(lag, acf)) +
  geom_col(color="black", fill="#252525", width = 0.5) +
  ylab("ACF") +
  xlab("Lag") +
  geom_hline(yintercept = qnorm((1+0.95)/2)/sqrt(auto_cor$n.used), 
             colour = "#cb181d",
             linetype = 2) +
  geom_hline(yintercept = -qnorm((1+0.95)/2)/sqrt(auto_cor$n.used), 
             colour = "#cb181d",
             linetype = 2) +
  theme(axis.text  = element_text(size=18),
        axis.title = element_text(size=18))

# -- Save supp-figure-2d
ggsave("figs/supp-figure-7d.pdf",
       plot   = supp_fig7d,
       dpi    = 300, 
       height = 5,
       width  = 7)

# -- Estimated sd
sd(V_inv %*% r)
### -- -------------------------------------------------------- -----------------------------------------------------
### -- END Supp Figure 7d: No correlation in adjusted residuals -----------------------------------------------------
### -- -------------------------------------------------------- -----------------------------------------------------