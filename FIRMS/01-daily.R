library(arrow)
library(sfarrow)
library(sf)
library(vroom)
library(lubridate)
library(dplyr)

firms_url <- "https://nrt3.modaps.eosdis.nasa.gov/api/v2/content/archives/FIRMS"

n20_dir <- "noaa-20-viirs-c2"
n20_pattern <- "J1_VIIRS_C2_Global_VJ114IMGTDL_NRT_%d%03d.txt"
modis_dir <- "modis-c6.1"
modis_pattern <- "MODIS_C6_1_Global_MCD14DL_NRT_%d%03d.txt"
snpp_dir <- "suomi-npp-viirs-c2"
snpp_pattern <- "SUOMI_VIIRS_C2_Global_VNP14IMGTDL_NRT_%d%03d.txt"

tryCatch(
  token <- readLines(".lance_token"),
  error = function(e) {
    stop(
      "File '.lance_token' not found in current working directory. ",
      "This script requires a LANCE token. ",
      "Create a LANCE token at https://nrt3.modaps.eosdis.nasa.gov/profile/app-keys ",
      "and put it in the '.lance_token' file to run this script."
    )
  }
)

get_date <- function(today_date) {
  # today_date <- d
  url_modis <- file.path(
    firms_url, modis_dir, "Global",
    sprintf(modis_pattern, year(today_date), yday(today_date))
  )
  url_snpp <- file.path(
    firms_url, snpp_dir, "Global",
    sprintf(snpp_pattern, year(today_date), yday(today_date))
  )
  url_n20 <- file.path(
    firms_url, n20_dir, "Global",
    sprintf(n20_pattern, year(today_date), yday(today_date))
  )

  datlist <- Map(
    get_data,
    c(url_modis, url_snpp, url_n20),
    c("MODIS", "SUOMI-NPP-VIIRS", "NOAA-20-VIIRS")
  )

  dat <- bind_rows(datlist) %>%
    mutate(
      across(starts_with("bright_"), as.numeric),
      acq_timestamp = acq_date + seconds(acq_time),
      acq_timestr = strftime(acq_timestamp, "%H%M%S", tz = "UTC"),
      year = year(acq_date),
      month = month(acq_date),
      mday = mday(acq_date)
    )

  tout <- tempfile(fileext = ".parquet")
  on.exit(file.remove(tout), add = TRUE)
  dat %>%
    st_write_parquet(tout)

  open_dataset(tout) %>%
    group_by(year, month, mday, version) %>%
    write_dataset(
      "processed-output/by-date-version",
      existing_data_behavior = "delete_matching"
    )
}

get_data <- function(url, instrument) {
  # url <- url_modis
  tfile <- tempfile(fileext=".txt", tmpdir = tempdir(check=TRUE))
  on.exit(file.remove(tfile), add = TRUE)
  download.file(
    url, tfile,
    headers = c("Authorization" = paste0("Bearer ", token))
  )
  cc <- cols(confidence = "c", satellite = "c", version = "c")
  raw <- vroom(tfile, col_types = cc)
  dsf <- st_as_sf(
    raw,
    coords = c("longitude", "latitude"),
    crs = st_crs(4326)
  )
  dsf[["instrument"]] <- instrument
  return(dsf)
}

alldat <- open_dataset("processed-output/by-date-version")
last_date <- alldat %>%
  summarize(max(acq_date)) %>%
  pull()
today_date <- Sys.Date()

if (last_date == today_date) {
  message("No dates to process!")
  if (!interactive()) quit(save = "no")
}

dateseq <- seq.Date(last_date, today_date, "1 day")

message("Retrieving dates ", dateseq[1], " -- ", tail(dateseq, 1),
        " (", length(dateseq), " total)")

for (i in seq_along(dateseq)) {
  # i <- 1
  d <- dateseq[i]
  message("Retrieving date: ", d)
  get_date(d)
}
