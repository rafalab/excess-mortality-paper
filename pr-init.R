library(lubridate)

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


# define control regions --------------------------------------------------

control_dates <- seq(as.Date("2002-01-01"), as.Date("2013-12-31"), by = "day")


# collapse age groups -----------------------------------------------------

data("puerto_rico_counts")

the_breaks <- c(0, 5, 20, 40, 60, 75, Inf)
all_counts <- collapse_counts_by_age(puerto_rico_counts, the_breaks)


