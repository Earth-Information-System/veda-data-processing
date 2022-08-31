#!/usr/bin/env python3
import glob
import datetime
import re

import xarray as xr
import pandas as pd
import numpy as np

from dask.diagnostics.progress import ProgressBar

prefix = "OCO2_GEOS_L3CO2_day"
target = f"{prefix}-v2.zarr"
procfile = "processed-files.txt"

files = sorted(list(glob.glob("./raw-data/*/oco2_GEOS_L3CO2_day_*.nc4")))
# files = files[0:60]

chunking = {
    "time": 100,
    "lat": 100,
    "lon": 100
}

print("Reading dataset...")
dat = xr.open_mfdataset(files)

ntime = len(dat.time)
residual_chunks = chunking["time"] - (ntime % chunking["time"])

print("Extending time series for even chunking...")
dates = dat.time.values
date_max = dates.max()
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
with ProgressBar():
    dat_extended.to_zarr(target, mode="w")

with open(procfile, "w") as f:
    f.writelines(file + "\n" for file in files)

print("Testing zarr")
dtest = xr.open_zarr(target)
dtest.sel(time = slice("2015-02-01", "2015-02-10"))
print("Done!")
