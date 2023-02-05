# 01-create-zarr-hour
# Create initial zarr for GEOS-5 hourly Fire Weather Index (FWI)
# Author: Katrina Sharonin

import re
import os
import glob
from datetime import datetime

import xarray as xr
import numpy as np
import pandas as pd

import dask

from dask.diagnostics.progress import ProgressBar
from tqdm.auto import tqdm

import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)

zar_file = 'FWI.GEOS-5.Hourly.zarr' # flag - may be different

chunking = {
    "lat": 100,
    "lon": 100,
    "time": 120,
}

# Keep all the dataset in the data directory
files = []
years = range(2017, int(datetime.now().year) + 1)

# base directory
basedir = "/lovelace/brewer/rfield1/storage/observations/GFWED/Sipongi/fwiCalcs.GEOS-5/Default/GEOS-5/"

# Parse paths and add existing hourly files
assert os.path.exists(basedir)
for y in years:
    assert os.path.exists(f"{basedir}/{y}")
    files += sorted(glob.glob(f"{basedir}/{y}/FWI.GEOS-5.Hourly.*.nc"))

d0 = xr.open_dataset(files[0])
dvars = list(d0.keys())
blank_array = dask.array.empty(
    (d0.time.shape[0], d0.lat.shape[0], d0.lon.shape[0]),
    dtype=np.float32
)

dims = ["time", "lat", "lon"]
blank_ds = xr.Dataset(
    data_vars = {k: (dims, blank_array) for k in dvars},
    coords = {
        "lat": d0.lat,
        "lon": d0.lon,
    }
)

def make_blank(pd_date):
    return (blank_ds.
            expand_dims({"time": [pd_date.replace(hour=n) for n in [*range(0,24)]]}, axis=1).
            transpose("time", "lat", "lon"))

dlist = []
# Sort files using annual directories
for file in tqdm(files):
    date_match = re.match(r".*\.(\d{8})\.nc", os.path.basename(file))
    assert date_match
    date = date_match.group(1)
    pd_date = pd.to_datetime(date)
    assert os.path.exists(file)

    try:
        dat = xr.open_mfdataset(file)
        hours = [*range(0,24)]
        dnew = dat.assign_coords({
            "time": [pd_date.replace(hour=n) for n in hours],
        })

        dlist.append(dnew)
    except Exception as err:
        print(f"Failed to open date {date} with error: `{str(err)}`")
        print("Skipping this timestep... Trying blank insert")
        dlist.append(make_blank(pd_date))
        continue

# Concatenate all files along time dimension
dat = xr.concat(dlist, dim="time")

# Extend time series to complete chunk length
ntime = len(dat.time)
residual_chunks = (chunking["time"] - (ntime % chunking["time"]))

print("Extending time series for even chunking...")
dates = dat.time.values
date_max = dates.max()
# add +1 day to max date
range_start = dates.max() + pd.Timedelta(1, 'D')
range_start = range_start.replace(hour=0)
# add remaining days using residual / 24
residual_days = int(residual_chunks / 24)
range_end = dates.max() + pd.Timedelta(residual_days, 'D')
range_end = range_end.replace(hour=23)
newdates = pd.date_range(
    range_start,
    range_end,
    freq='1H'
)

assert len(newdates) == residual_chunks
alldates = np.concatenate((dates, newdates))
assert len(alldates) % chunking["time"] == 0
dat_extended = dat.reindex({"time": alldates}, fill_value = np.nan)
assert len(dat_extended.time) % chunking["time"] == 0
dat_extended = dat_extended.chunk(chunking)

# Write to file + progress bar
with ProgressBar():
    dat_extended.to_zarr(zar_file, mode="w")

# Test output
# See Katrina Sharonin's local testing module

