# pak::pak(c("sfarrow"))
library(sfarrow)
library(sf)
library(dplyr)

outfile <- strftime(Sys.time(), "processed-output/FIRMS-%Y-%m-%d-%H%M%S.parquet")
message("Writing daily output to: ", outfile)

url_modis <- "https://firms.modaps.eosdis.nasa.gov/data/active_fire/modis-c6.1/csv/MODIS_C6_1_Global_24h.csv"
url_vsnpp <- "https://firms.modaps.eosdis.nasa.gov/data/active_fire/suomi-npp-viirs-c2/csv/SUOMI_VIIRS_C2_Global_24h.csv"
url_vn20 <- "https://firms.modaps.eosdis.nasa.gov/data/active_fire/noaa-20-viirs-c2/csv/J1_VIIRS_C2_Global_24h.csv"

get_data <- function(url, instrument) {
  raw <- read_sf(url)
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
  c(url_modis, url_vsnpp, url_vn20),
  c("MODIS", "SUOMI-NPP-VIIRS", "NOAA-20-VIIRS")
)

dat <- bind_rows(datlist) %>%
  mutate(across(starts_with("bright_"), as.numeric))

st_write_parquet(dat, outfile)
