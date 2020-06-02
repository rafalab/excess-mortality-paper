### -- --------------------------------------------------------- -----------------------------------------------------
### -- Supp Figure 7a: QQ-plot without adjusting for correlation -----------------------------------------------------
### -- --------------------------------------------------------- -----------------------------------------------------
# -- Set up 
source("pr-init.R")
counts <- filter(all_counts, agegroup == "75-Inf")

# -- Computing expected mortality counts
counts <- compute_expected(counts, exclude = exclude_dates, weekday.effect = TRUE)

# -- Dates for model check
example_dates <- control_dates[year(control_dates) > 2005]
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
  ylab("Sample quantiles") +
  xlab("Theoretical quantiles") +
  scale_y_continuous(limits = c(-4, 5),
                     breaks = seq(-4, 5, by=2)) +
  scale_x_continuous(limits = c(-4, 4),
                     breaks = seq(-4, 4, by=2)) +
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
  scale_y_continuous(limits = c(-4, 5),
                     breaks = seq(-4, 5, by=2)) +
  scale_x_continuous(limits = c(-4, 4),
                     breaks = seq(-4, 4, by=2)) +
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