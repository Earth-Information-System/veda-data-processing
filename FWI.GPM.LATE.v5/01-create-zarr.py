#!/usr/bin/env python3
import glob
import datetime
import re

import xarray as xr
import pandas as pd
import zarr

from dask.diagnostics.progress import ProgressBar

prefix = "FWI_GPM_LATE_v5_Daily-test"
files = sorted(list(glob.glob("raw-input/2022/*.nc")))
# files = files[0:20]

def parse_date(path):
    ms = re.search(r"\.(\d{4})(\d{2})(\d{2})\.nc$", path)
    year = int(ms.group(1))
    month = int(ms.group(2))
    day = int(ms.group(3))
    date = datetime.date(year, month, day)
    return date

dat_raw = xr.open_mfdataset(files, combine="nested", concat_dim="time")
# The time coordinate here is wrong, so we reassign it by parsing the file name
dates = pd.to_datetime([parse_date(f) for f in files])
dat = dat_raw.assign_coords(time=dates).chunk({
    "time": 100,
    "lat": 100,
    "lon": 100
})

target = f"{prefix}.zarr"
with ProgressBar():
    dat.to_zarr(target, mode="w")

print("testing zarr")
dtest = xr.open_zarr(target)
