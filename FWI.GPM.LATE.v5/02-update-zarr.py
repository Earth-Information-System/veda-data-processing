import datetime
import re

import xarray as xr
import pandas as pd
import zarr

timevar = "time"
zarrpath = "FWI_GPM_LATE_v5_Daily-test.zarr"
infile = "raw-input/2022/FWI.GPM.LATE.v5.Daily.Default.20220121.nc"

lastday = xr.open_zarr(zarrpath)[timevar].max().values

def parse_date(path):
    ms = re.search(r"\.(\d{4})(\d{2})(\d{2})\.nc$", path)
    year = int(ms.group(1))
    month = int(ms.group(2))
    day = int(ms.group(3))
    date = datetime.date(year, month, day)
    return date

dnew_raw = xr.open_dataset(infile)
dnew = dnew_raw\
    .squeeze("time")\
    .drop("time")\
    .expand_dims(time=[pd.to_datetime(parse_date(infile))], axis=0)

dnew.to_zarr(zarrpath, append_dim="time")
zarr.consolidate_metadata(zarrpath)

dtest = xr.open_zarr(zarrpath, consolidated=False)
dtest.chunks

# Update algorithm
# (1) On creation, fill the last time dimension to the end with `nan`s
# (2) On update:
    # (a) Open existing store. Check last timestamp and time chunk size.
    # (b1) If timestamp in file, replace.
    # (b2) If not, create new chunk of same size and replace.
