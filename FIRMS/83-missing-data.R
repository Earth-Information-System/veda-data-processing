library(arrow)
library(sfarrow)
library(dplyr)
library(tidyr)
library(lubridate)
library(forcats)

dat <- open_dataset("processed-output/by-date-version")

date_instrument <- dat %>%
  select(acq_date, satellite, instrument, version) %>%
  distinct() %>%
  collect() %>%
  as_tibble() %>%
  arrange(acq_date)

date_instrument %>%
  count(satellite, instrument)

date_instrument2 <- date_instrument %>%
  mutate(
    satellite = fct_recode(
      satellite,
      "Aqua" = "A", "Terra" = "T",
      "NOAA-20" = "1", "SUOMI-NPP" = "N"
    ),
    instrument = case_when(
      instrument == "MODIS" ~ paste(instrument, satellite),
      instrument == "VIIRS" ~ paste(satellite, instrument, sep = "-"),
      TRUE ~ instrument
    )
  ) %>%
  select(acq_date, instrument)

all_date_instrument <- date_instrument2 %>%
  complete(acq_date, instrument)

missing_date_instrument <- all_date_instrument %>%
  anti_join(date_instrument2)

missing_date_instrument %>%
  filter(grepl("VIIRS", instrument)) %>%
  arrange(desc(acq_date))

missing_date_instrument %>%
  filter(acq_date >= "2022-01-01") %>%
  arrange(desc(acq_date)) %>%
  print(n = 40)

date_instrument %>%
  count(acq_date, sort = TRUE) %>%
  count(n)

date_instrument %>%
  summarize(max(acq_date))

date_instrument %>%
  group_by(acq_date) %>%
  count(sort = TRUE) %>%
  filter(n == 3) %>%
  pull(acq_date)

zdsub <- dat %>%
  filter(acq_date == as.Date("2019-12-04")) %>%
  collect() %>%
  as_tibble()

zdsub %>%
  count(satellite, instrument)

dsub <- dat %>%
  filter(
    acq_timestamp > as.Date("2022-01-01"),
    acq_timestamp < as.Date("2022-01-10")
  ) %>%
  collect() %>%
  as_tibble()

dsub %>%
  count(instrument, satellite)

dat %>%
  group_by(mday, version) %>%
  write_dataset(
    "processed-output/by-date-version/year=2022/month=5",
    existing_data_behavior = "delete_matching",
    coerce_timestamps = "ms"
  )

# dat <- open_dataset("processed-output/by-date-version/year=2022/month=4") %>%
#   collect() %>%
#   as_tibble()

# for (mday in seq(2, 30)) {
#   # mday <- 3
#   dat <- open_dataset(paste0("processed-output/by-date-version/year=2022/month=6/mday=", mday)) %>%
#     collect() %>%
#     as_tibble()
# }

# dat <- open_dataset("processed-output/by-date-version/year=2022") %>%
#   filter(acq_date >= as.Date("2022-06-02")) %>%
#   collect()

#   select(acq_date, acq_time) %>%
#   collect() %>%
#   as_tibble()
# head(dat)

##################################################

dtime <- open_dataset(
  "processed-output/by-date-stability"
)
dtime %>% 
  select(acq_date, acq_time) %>%
  collect() %>% as_tibble()

dat <- open_dataset(
  "processed-output/by-date-stability",
  schema(
    brightness = double(),
    scan = double(),
    track = double(),
    acq_date = date32(),
    acq_time = string(),
    satellite = string(),
    instrument = string(),
    confidence = string(),
    version = string(),
    bright_t31 = double(),
    frp = double(),
    daynight = string(),
    type = double(),
    geometry = binary()
  )
)

dat %>% head() %>% as_tibble() %>% glimpse()

dtime <- dat %>%
  select(acq_date, acq_time) %>%
  as_tibble()

d2 <- dat %>%
  mutate(
    h1 = cast(substring(acq_time, 1, 2), int64()),
    m1 = cast(substring(acq_time, 3, 4), int64()),
    t1 = cast(h1 * 3600 + m1 * 60, int64()),
    d1 = cast(acq_date, int32()),  # Days since 1970-01-01
    d2 = cast((d1 * 86400 + t1), int64()),  # Seconds since 1970-01-01
    acq_t = cast(d2, timestamp("s", "UTC")),
    acq_timestamp = cast(acq_t, string()),
    year = year(acq_date),
    month = month(acq_date),
    mday = mday(acq_date)
  ) %>%
  select(-h1, -m1, -t1, -d1, -d2, -acq_t) %>%
  select(acq_date, acq_timestamp, everything())

d2 %>% head() %>% as_tibble() %>% glimpse()

d2 %>%
  group_by(year, month, mday) %>%
  write_dataset(
    "processed-output/FIRMS-ymd",
    existing_data_behavior = "delete_matching"
  )

dtest <- open_dataset("processed-output/FIRMS-ymd")
