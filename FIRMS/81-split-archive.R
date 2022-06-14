library(sfarrow)
library(sf)
library(arrow)
library(dplyr)
library(lubridate)

parse_timestamp <- function(dat) {
  dat %>%
    mutate(
      h1 = cast(substring(acq_time, 1, 2), int64()),
      m1 = cast(substring(acq_time, 3, 4), int64()),
      t1 = cast(h1 * 3600 + m1 * 60, int64()),
      d1 = cast(acq_date, int32()),  # Days since 1970-01-01
      d2 = cast((d1 * 86400 + t1), int64()),  # Seconds since 1970-01-01
      acq_timestamp = cast(d2, timestamp("s", "UTC")),
      year = year(acq_date),
      month = month(acq_date),
      mday = mday(acq_date)
    ) %>%
    select(-h1, -m1, -t1, -d1, -d2)
}

process_firms <- function(dat) {
  dat %>%
    parse_timestamp() %>%
    group_by(year, month, mday, version) %>%
    write_dataset(
      "processed-output/by-date-version",
      existing_data_behavior = "delete_matching"
    )
}

write_firms <- function(infile) {
  open_dataset(infile) %>%
    process_firms()
}

for (year in 2019:2021) {
  message("Year: ", year)
  write_firms(paste0("processed-output/FIRMS-archive-", year, ".parquet"))
}

for (year in 2019:2022) {
  message("Year: ", year)
  write_firms(paste0("processed-output/FIRMS-NRT-", year, ".parquet"))
}

# dailies <- list.files(
#   "processed-output",
#   "FIRMS-daily-.*",
#   full.names = TRUE
# )

# for (f in dailies) {
#   # f <- dailies[1]
#   message("Processing daily file: ", f)
#   d <- read_parquet(f) %>%
#     as_tibble() %>%
#     mutate(acq_timestamp = )
#     process_firms()
# }

# alldat <- open_dataset("processed-output/by-date-version")

# alldat %>%
#   summarize(max(acq_date)) %>%
#   collect()

# alldat %>%
#   filter(date == "2022-06-03")

# library(ggplot2)
# library(ggspatial)

# ggplot(d2) +
#   aes(color = frp) +
#   geom_sf()

# d2s <- read_sf_dataset()

# vals <- d2 %>%
#   filter(acq_date == "2019-01-01") %>%
#   as_tibble() %>%
#   st_as_sf()


# ?write_parquet

# date <- seq.Date(as.Date("2019-01-01"), as.Date("2019-12-31"), by = "1 day")

# arch1
