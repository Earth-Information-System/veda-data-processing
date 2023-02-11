# 02-update-zarr-hour.py
# If an GEOS-5 hourly FWI zarr exists, update content with recent files
# Author: Katrina Sharonin

import datetime
import re
import glob
import sys
from datetime import datetime
from os.path import exists

import dask
import os
from tqdm.auto import tqdm
import netCDF4 as nc

import xarray as xr
import pandas as pd
import numpy as np


timevar = "time"
zarrpath = 'FWI.GEOS-5.Hourly.zarr' # flag - may be different
procfile = "processed-files-bak-hourly.txt" # flag - may be different

allfiles = []
basedir = "/lovelace/brewer/rfield1/storage/observations/GFWED/Sipongi/fwiCalcs.GEOS-5/Default/GEOS-5/"
assert os.path.exists(basedir)
years = range(2017, int(datetime.now().year) + 1)
years = range(2022, 2023) # TODO - remove limit after testing

for y in years:
    assert os.path.exists(f"{basedir}/{y}")
    allfiles += sorted(glob.glob(f"{basedir}/{y}/FWI.GEOS-5.Hourly.*.nc"))

# Open txt file to id current files in zarr
with open(procfile, "r") as f:
    procfiles = [l for l in f.read().splitlines() if l != '']
newfiles = sorted(set(allfiles) - set(procfiles))

# Terminate program if no new files
if len(newfiles) == 0:
    sys.exit("No new files to process!")

# Function to extract date from file string
def parse_date(path):
    ms = re.search(r"\.(\d{4})(\d{2})(\d{2})\.nc$", path)
    year = int(ms.group(1))
    month = int(ms.group(2))
    day = int(ms.group(3))
    date = datetime.date(year, month, day)
    return date

# given list of new files, create DataArray with new data
d0 = xr.open_dataset(newfiles[0]) # read-only access
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


# base directory
dlist = []

for file in tqdm(newfiles):
    date_match = re.match(r".*\.(\d{8})\.nc", os.path.basename(file))
    assert date_match
    date = date_match.group(1)
    pd_date = pd.to_datetime(date)
    assert os.path.exists(file)

    try:
        dat = xr.open_mfdataset(file)
        hours = [*range(0,24)]
        dnew = dat.assign_coords({
            "time": [pd_date.replace(hour=n) for n in hours]
        })
        dlist.append(dnew)
    except Exception as err:
        print(f"Failed to open date {date} with error: `{str(err)}`")
        print("Skipping this timestep...Trying blank insert")
        dlist.append(make_blank(pd_date))
        continue

# Concatenate all new files via time dimension
dat = xr.concat(dlist, dim="time")
dnew_all = dat

chunking = {
    "lat": 100,
    "lon": 100,
    "time": 120
}

# chunk data
dnew_all = dnew_all.chunk(chunking)

# Given current zarr, look for timestamps of the new files
# if timestamp appears once -> no new extension needed
# if timestamp not present -> extend by chunk size to force presence
for idt, t in enumerate(dnew_all.time):

    # skip over any non-0 hour time, only iterate on every 24 hours of data
    if t.values not in dnew_all.time[0::24]:
        continue

    # select single time dim
    dnew = dnew_all.isel(time=slice(idt, idt + 24))
    assert dnew_all.time[idt].values == dnew.time.values[0], "Mismatching enumeration"
    print(dnew.time.values[0])
    # check within current zarr if time exists
    dtarget = xr.open_zarr(zarrpath)
    is_in = dtarget.time.isin(dnew.time[0])
    n_is_in = is_in.sum()

    if n_is_in > 1:
        raise Exception(f"Found {n_is_in} matching times! Something is wrong with this situation.")

    # timestamp with date exists -> no need for block extension
    elif n_is_in == 1:
        # identify location of the timestamp
        itime = int(np.where(is_in)[0])
        assert dtarget.isel(time=itime).time == dnew.time[0], "Mismatch in times"
        print(f"Inserting into time location {itime}")

        # slice and insert all hours
        dnew.to_zarr(zarrpath, region={
            "time": slice(itime, itime + 24),
            "lat": slice(0, dnew.sizes["lat"]),
            "lon": slice(0, dnew.sizes["lon"]),
        })

    # timestamp does NOT exist -> must extend our blocks to include time stamp
    elif n_is_in == 0:

        tchunk = dtarget.chunks["time"][-1]
        print(f"Creating new time chunk of size {tchunk}...")

        # add +1 day to max date
        range_start = dtarget.time.max().values + pd.Timedelta(1, 'D')
        range_start = range_start.replace(hour=0)
        # add remaining days using residual days / 24
        residual_days = int(tchunk/24)
        range_end = dtarget.time.max().values + pd.Timedelta(residual_days, 'D')
        range_end = range_end.replace(hour=23)

        newdates = pd.date_range(
            range_start,
            range_end,
            freq='1H'
        )

        print(f"...from {newdates.min()} to {newdates.max()}.")
        assert len(newdates) == tchunk
        # Extend the current timestep out to size `tchunk`, filling with `Nan`s
        dummy = dnew.reindex({"time": newdates}, fill_value=np.nan)
        assert len(dummy.time) == tchunk

        # Append to the existing Zarr data store
        dummy.to_zarr(zarrpath, append_dim="time")

print('Update zarr process complete!')

# print('Initiate testing...')
# Test output
# See Katrina Sharonin's local testing module

