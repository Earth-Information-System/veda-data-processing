# Active fire detections from FIRMS

[NASA FIRMS](https://firms.modaps.eosdis.nasa.gov) provides conveniently extracted (e.g., CSV) versions of active fire detection data from MODIS (Terra, Aqua) and VIIRS (Suomi-NPP and NOAA-20).
Recent (<7 days old) data can be downloaded directly and are updated at least daily.
Archival data (> 7 days old) need to be requested.
For performance reasons, archival data can only be requested up to 1 year at a time.

This directory has several data processing scripts:

- `00-archive.R` -- This is a script that reads the raw archive data (in CSV format) and produces yearly files for the archive (`-archive-`) and near-real-time (`-NRT-`) algorithm versions. The resulting files are quite large and inefficient to process, so a subsequent script (`81-split-archive.R`) splits these up using the Parquet Arrow "partitioning" approach according to year, month, month day (`mday`), and algorithm version -- e.g., `processed-output/by-date-version/year=2020/month=6/mday=3/version=2.0NRT`.
  - Parquet readers can access any part of this directory hierarchy to read all of the corresponding data in a single pass. E.g., in R, all of the following will work:
  ```r
  # Read only 2020-06-03
  sfarrow::st_read_parquet("processed-output/by-date-version/year=2020/month=6/mday=3")
  # Read all data from 2021
  sfarrow::st_read_parquet("processed-output/by-date-version/year=2021")
  # Read all data
  sfarrow::st_read_parquet("processed-output/by-date-version")
  ```
- `01-daily.R` -- Process data on a daily basis. This script assumes the existence of `processed-output/by-date-version`. It will first read this dataset, identify the last available date, download all the data between that day and the present, and append it to the archive in the same format.

Since we are treating these as georeferenced point observations, we store them in geoparquet format, which is an efficient, cloud-optimized columnar data store.
See below for example code for working with `geoparquet` data.
See also the `pyarrow-test.py` script.

Raw input data (mostly, yearly archives for `archive.R`) are stored in `raw-input/`.
Processed parquet files produced by both scripts are stored in `processed-output/`.

## Working with Geoparquet data

```python
import pyarrow
import pyarrow.dataset as pds
import geopandas as gpd
```

**Option 1**:
Use `pyarrow` to read from the entire dataset; then "export" to GeoPandas.

```python
basepath = "s3://veda-data-store-staging/EIS/geoparquet/FIRMS-active-fires"
pdat = pds.dataset(basepath)
dat = pdat.head(10).to_pandas()
gs = gpd.GeoSeries.from_wkb(dat["geometry"])
gdat = gpd.GeoDataFrame(dat, geometry=gs, crs="EPSG:4326")
```

**Option 2**:
Read directly from parquet using the external partitioning to do temporal subsetting.

```python
d2 = gpd.read_parquet(
    "s3://veda-data-store-staging/EIS/geoparquet/FIRMS-active-fires/year=2022/month=1/mday=15"
)
```
