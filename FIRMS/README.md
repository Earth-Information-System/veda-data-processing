# Active fire detections from FIRMS

[NASA FIRMS](https://firms.modaps.eosdis.nasa.gov) provides conveniently extracted (e.g., CSV) versions of active fire detection data from MODIS (Terra, Aqua) and VIIRS (Suomi-NPP and NOAA-20).
Recent (<7 days old) data can be downloaded directly and are updated at least daily.
Archival data (> 7 days old) need to be requested.
For performance reasons, archival data can only be requested up to 1 year at a time.

Here, we have two scripts for converting:

- `daily.R` processes the daily data. This should be run daily, and reads data directly from FIRMS. Output files follow the pattern: `FIRMS-YYYY-MM-DD-HHMMSS.parquet`, where `YYYY-MM-DD-HHMMSS` is the timestamp at which the file was run.
- `archive.R` processes the archive data. This should only be run once (or whenever the archival data are downloaded). This produces two output files:
  - `FIRMS-archive.parquet` -- These are re-processed archive versions of the data that will not be further altered. These data are typically older than the NRT data (see below).
  - `FIRMS-NRT.parquet` -- These use a more provisional algorithm intended for rapid processing for near-real-time (NRT) applications, but are less likely to be accurate. Over time, NRT data will transition to archive.

Since we are treating these as georeferenced point observations, we store them in geoparquet format, which is the preferred vector format for the VEDA team.

Raw input data (mostly, yearly archives for `archive.R`) are stored in `raw-input/`.
Processed parquet files produced by both scripts are stored in `processed-output/`.
