### -- -------------------------------------------- ------------------------------------------------------------------
### -- Supp figure 1: Population displacement in PR ------------------------------------------------------------------
### -- -------------------------------------------- ------------------------------------------------------------------
# -- Set up
source("pr-init.R")

supp_fig <- puerto_rico_counts %>%
  group_by(date) %>%
  summarize(population = sum(population)) %>%
  ungroup() %>%
  filter(date >= "2017-07-01", date <= "2018-07-01") %>%
  ggplot(aes(date, population)) +
  geom_line() +
  scale_y_continuous(labels = scales::comma,
                     limits = c(2981110, 3330000),
                     breaks = seq(3000000, 3300000, by=50000)) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") +
  ylab("Population") +
  xlab("") +
  theme(axis.text  = element_text(size=15),
        axis.title = element_text(size=18))

# -- Supp figure 1
ggsave("figs/supp-figure-1-population-pr.pdf",
       plot   = supp_fig,
       dpi    = 300, 
       height = 5,
       width  = 7)
### -- -------------------------------------------- ------------------------------------------------------------------
### -- END Supp figure 1: Population displacement in PR ------------------------------------------------------------------
### -- -------------------------------------------- ------------------------------------------------------------------
