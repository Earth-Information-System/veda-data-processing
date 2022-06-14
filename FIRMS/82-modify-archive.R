library(arrow)
library(sfarrow)
library(dplyr)
library(lubridate)

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
