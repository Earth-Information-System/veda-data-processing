# Convert archive NRT files (as downloaded from FIRMS "retrieve archive"
# service) to daily files (as downloaded directly from MODAPS, and as expected
# by FEDS code).

library(tidyverse)
library(vroom)
library(forcats)
library(lubridate)

reference_modis_path <- "reference-files/MODIS_C6_1_Global_MCD14DL_NRT_2022253.txt"
reference_snpp_path <- "reference-files/SUOMI_VIIRS_C2_Global_VNP14IMGTDL_NRT_2022195.txt"
reference_n20_path <- "reference-files/J1_VIIRS_C2_Global_VJ114IMGTDL_NRT_2022195.txt"
outdir <- "plain-text-daily"

write_daily <- function(dat, grp, prefix) {
  date <- grp[["acq_date"]]
  filename <- sprintf("%s%04d%03d.txt", prefix, year(date), yday(date))
  vroom_write(dat, filename, delim = ",")
}

ref_n20 <- vroom(reference_n20_path, col_types = cols(
  acq_time = "c"
))
archive_n20 <- vroom("raw-input/fire_nrt_J1V-C2_274164.csv")

head(ref_n20, 3)
head(archive_n20, 3)

ref_n20 %>% distinct(confidence)
archive_n20 %>% distinct(confidence)

archive_n20_reformat <- archive_n20 %>%
  select(-instrument) %>%
  mutate(
    confidence = factor(confidence) %>%
      fct_recode(nominal = "n", high= "h", low = "l"),
    acq_time = paste0(substr(acq_time, 1, 2), ":", substr(acq_time, 3, 4))
  )

head(ref_n20, 3)
head(archive_n20_reformat, 3)

prefix_n20 <- file.path(outdir, "VJ114IMGTDL", "J1_VIIRS_C2_Global_VJ114IMGTDL_NRT_")
dir.create(dirname(prefix_n20), recursive = TRUE, showWarnings = FALSE)
archive_n20_reformat %>%
  group_by(acq_date) %>%
  group_walk(write_daily, prefix = prefix_n20)

# ref_modis <- vroom(reference_modis_path, col_types = cols(
#   acq_time = "c"
# ))
# archive_mod <- vroom("raw-input/fire_nrt_M-C61_274163.csv")

# archive_modis_reformat <- archive_mod %>%
#   select(-instrument) %>%
#   mutate(
#     satellite = factor(satellite) %>%
#       fct_recode("A" = "Aqua", "T" = "Terra"),
#     acq_time = paste0(substr(acq_time, 1, 2), ":", substr(acq_time, 3, 4))
#   )

# prefix_modis <- file.path(outdir, "MODIS", "MODIS_C6_1_Global_MCD14DL_NRT_")
# dir.create(dirname(prefix_modis), recursive = TRUE)

# archive_modis_reformat %>%
#   group_by(acq_date) %>%
#   group_map(vroom_write, write_daily)

# archive_mod %>%
#   distinct(acq_time) %>%
#   print(n = 100)

# ref_modis %>%
#   distinct(acq_time) %>%
#   tail(100)

