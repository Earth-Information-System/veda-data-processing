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

today_date <- Sys.Date()
argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) > 0) {
  # argv <- "2022-06-02"
  today_date <- as.Date(argv[1])
}

outfile <- strftime(today_date, "processed-output/FIRMS-daily-%Y-%m-%d.parquet")
message("Writing daily output to: ", outfile)

if (file.exists(outfile)) {
  stop("Output file already exists.")
}

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

get_data <- function(url, instrument) {
  # url <- url_modis
  tfile <- tempfile(fileext=".txt")
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

datlist <- Map(
  get_data,
  c(url_modis, url_snpp, url_n20),
  c("MODIS", "SUOMI-NPP-VIIRS", "NOAA-20-VIIRS")
)

dat <- bind_rows(datlist) %>%
  mutate(across(starts_with("bright_"), as.numeric))

st_write_parquet(dat, outfile)
