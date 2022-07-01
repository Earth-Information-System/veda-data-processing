#!/usr/bin/env python3
import glob
import datetime
import re

import xarray as xr
import pandas as pd
import numpy as np

from dask.diagnostics.progress import ProgressBar

prefix = "FWI_GPM_LATE_v5_Daily"
files = sorted(list(glob.glob("raw-input/*/*.Daily.*.nc")))
# files = files[:-3]

chunking = {
    "time": 50,
    "lat": 100,
    "lon": 100
}

def parse_date(path):
    ms = re.search(r"\.(\d{4})(\d{2})(\d{2})\.nc$", path)
    year = int(ms.group(1))
    month = int(ms.group(2))
    day = int(ms.group(3))
    date = datetime.date(year, month, day)
    return date

print("Reading dataset...")
dat_raw = xr.open_mfdataset(files, combine="nested", concat_dim="time")
# The time coordinate here is wrong, so we reassign it by parsing the file name
print("Fixing time coordinate...")
dates = pd.to_datetime([parse_date(f) for f in files])
dat = dat_raw.assign_coords(time=dates)

ntime = len(dat.time)
residual_chunks = chunking["time"] - (ntime % chunking["time"])

print("Extending time series for even chunking...")
# dummy = xr.full_like(dat.isel(time = slice(residual_chunks)), np.nan)
# assert len(dummy.time) == residual_chunks
newdates = pd.date_range(
    dates.max() + pd.Timedelta(1, 'D'),
    dates.max() + pd.Timedelta(residual_chunks, 'D')
)
assert len(newdates) == residual_chunks
alldates = np.concatenate((dates, newdates))
assert len(alldates) % chunking["time"] == 0
dat_extended = dat.reindex({"time": alldates}, fill_value = np.nan)
assert len(dat_extended.time) % chunking["time"] == 0
dat_extended = dat_extended.chunk(chunking)

print("Writing output zarr...")
target = f"{prefix}.zarr"
with ProgressBar():
    dat_extended.to_zarr(target, mode="w")

with open("processed-files.txt", "w") as f:
    f.writelines(files)

print("Testing zarr")
dtest = xr.open_zarr(target)
print("Done!")
