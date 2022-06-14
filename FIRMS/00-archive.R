library(sfarrow)
library(sf)
library(vroom)
library(dplyr)
library(purrr)
library(lubridate)

outfile <- "processed-output/FIRMS-archive.parquet"
csvs <- list.files("raw-input", "*.csv", full.names = TRUE)

cc <- cols(confidence="c", satellite="c", version="c")

# Process the archive and NRT files separately 

message("Reading raw arch CSVs...")
raw_arch <- vroom(grep("_archive_", csvs, value = TRUE), col_types = cc)

# Split up the files by year. This makes it easier to track, and we can treat
# them as a single dataset later with Arrow functions.
years <- sort(unique(year(raw_arch$acq_date)), decreasing = TRUE)

to_sf <- function(dat, ...) {
  st_as_sf(dat, coords = c("longitude", "latitude"), crs = st_crs(4326))
}

for (year in sort(years, decreasing = TRUE)) {
  message("Processing year: ", year)
  outfile <- sprintf("processed-output/FIRMS-archive-%d.parquet", year)
  if (file.exists(outfile)) {
    message("Skipping year because file exists: ", outfile)
    next
  }
  raw_arch_sub <- raw_arch %>%
    filter(year(acq_date) == year)
  sf_arch <- to_sf(raw_arch_sub)
  st_write_parquet(sf_arch, outfile)
}

rm(raw_arch, sf_arch)

message("Reading raw NRT CSVs...")
raw_nrt <- vroom(grep("_nrt_", csvs, value = TRUE), col_types = cc) %>%
  mutate(type = NA_real_)

years <- sort(unique(year(raw_nrt$acq_date)), decreasing = TRUE)

for (year in years) {
  message("Processing year: ", year)
  outfile <- sprintf("processed-output/FIRMS-NRT-%d.parquet", year)
  if (file.exists(outfile)) {
    message("Skipping year because file exists: ", outfile)
    next
  }
  raw_sub <- raw_nrt %>%
    filter(year(acq_date) == year)
  sf_arch <- to_sf(raw_sub)
  st_write_parquet(sf_arch, outfile)
}

message("Done!")
quit(save="no")
