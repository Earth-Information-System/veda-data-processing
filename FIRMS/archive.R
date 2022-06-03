library(sfarrow)
library(sf)
library(vroom)
library(dplyr)
library(purrr)

outfile <- "processed-output/FIRMS-archive.parquet"
csvs <- list.files("raw-input", "*.csv", full.names = TRUE)

cc <- cols(confidence="c", satellite="c", version="c")

# Process the archive and NRT files separately 
message("Reading raw arch CSVs...")
raw_arch <- vroom(grep("_archive_", csvs, value = TRUE), col_types = cc)

message("Converting archive dataset to SF...")
sf_arch <- st_as_sf(raw_arch, coords = c("longitude", "latitude"), crs = st_crs(4326))
message("Writing out archive parquet file...")
st_write_parquet(sf_arch, "processed-output/FIRMS-archive.parquet")

rm(raw_arch, sf_arch)

message("Reading raw NRT CSVs...")
raw_nrt <- vroom(grep("_nrt_", csvs, value = TRUE), col_types = cc) %>%
  mutate(type = NA_real_)
message("Converting NRT dataset to SF...")
sf_nrt <- st_as_sf(raw_nrt, coords = c("longitude", "latitude"), crs = st_crs(4326))
message("Writing out NRT file...")
st_write_parquet(sf_nrt, "processed-output/FIRMS-NRT.parquet")

message("Done!")
quit(save="no")
