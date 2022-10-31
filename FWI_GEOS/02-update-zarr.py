# Non-localized version of update zarr for forecast

import datetime
import re
import glob
import sys
import dask
import os
from tqdm.auto import tqdm
import netCDF4 as nc

import xarray as xr
import pandas as pd
import numpy as np

timevar = "time"
zarrpath = 'FWI.GEOS-5.zarr'

procfile = "processed-files-bak.txt" # flag - unsure if this would be different
allfiles = sorted(list(glob.glob("raw-input/*/*.Daily.*.nc")))  # flag - unsure if this would be different
with open(procfile, "r") as f:
    # Remove empty lines
    procfiles = [l for l in f.read().splitlines() if l != '']
newfiles = sorted(set(allfiles) - set(procfiles))

# if no new files, exit program
if len(newfiles) == 0:
    sys.exit("No new files to process!")

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
    (9, d0.lat.shape[0], d0.lon.shape[0]),
    dtype=np.float32
)
# include forecast dimension
dims = ["forecast", "lat", "lon"]
blank_ds = xr.Dataset(
    data_vars = {k: (dims, blank_array) for k in dvars},
    coords = {
        "lat": d0.lat,
        "lon": d0.lon,
        "forecast": list(range(0, 9))
    }
)

def make_blank(pd_date):
    return (blank_ds.
            expand_dims({"time": [pd_date]}, axis=1).
            transpose("forecast", "time", "lat", "lon"))

dlist = []
basedir = "/autofs/brewer/rfield1/storage/observations/GFWED/Sipongi/fwiCalcs.GEOS-5/Default/GEOS-5" # flag - unsure if this would be different

for file in tqdm(newfiles):
    date_match = re.match(r".*\.(\d{8})\.nc", os.path.basename(file))
    assert date_match
    date = date_match.group(1)
    pd_date = pd.to_datetime(date)
    year = date[0:4]
    datedir = f"{basedir}/{year}/{date}00"
    if not os.path.exists(datedir):
        print(f"Not found: {datedir}")
        print("Skipping this timestep.")
        dlist.append(make_blank(pd_date))
        continue
    # Get list of forecast files.
    pred_files = sorted(glob.glob(f"{datedir}/*.nc"))
    if not len(pred_files):
        print(f"No files found for date {date}. Skipping this timestep.")
        dlist.append(make_blank(pd_date))
        continue
    try:
        dat = xr.open_mfdataset([file] + pred_files, combine="nested", concat_dim="forecast")
        dnew = dat.assign_coords({
            "time": [pd_date],
            "forecast": range(0, 1+len(pred_files))
        })
        dlist.append(dnew)
    except Exception as err:
        print(f"Failed to open date {date} with error: `{str(err)}`")
        print("Skipping this timestep...")
        dlist.append(make_blank(pd_date))
        continue

# concatenate all new files via time dim
dat = xr.concat(dlist, dim="time")
dnew_all = dat

chunking = {
    "lat": 100,
    "lon": 100,
    "time": 100,
    "forecast": 9
}

# chunk data
dnew_all = dnew_all.chunk(chunking)

# now given "main" zarr, look for timestamps of the new files
# if timestamp appears once -> no new extension needed
# if timestamp not present -> extend to force presence
for idt, t in enumerate(dnew_all.time):

    # select single time dim
    dnew = dnew_all.isel(time=slice(idt, idt + 1))
    assert dnew_all.time[idt].values == dnew.time.values[0], "Mismatching enumeration"
    print(dnew.time.values[0])
    # check within current zarr if time exists
    dtarget = xr.open_zarr(zarrpath)
    is_in = dtarget.time.isin(dnew.time)
    n_is_in = is_in.sum()

    # there should only exist one unique time stamp
    if n_is_in > 1:
        raise Exception(f"Found {n_is_in} matching times! Something is wrong with this situation.")

    # timestamp exists -> no need for block extension
    elif n_is_in == 1:
        # identify location of the timestamp
        itime = int(np.where(is_in)[0])
        assert dtarget.isel(time=itime).time == dnew.time, "Mismatch in times"
        print(f"Inserting into time location {itime}")

        # slice and insert
        dnew.to_zarr(zarrpath, region={
            "time": slice(itime, itime + 1),
            "lat": slice(0, dnew.sizes["lat"]),
            "lon": slice(0, dnew.sizes["lon"]),
            # added slice - check if appropriately placed
            "forecast": slice(0, 9)
        })

    # timestamp DNE -> must extend our blocks to include time stamp
    elif n_is_in == 0:

        tchunk = dtarget.chunks["time"][-1]
        print(f"Creating new time chunk of size {tchunk}...")

        newdates = pd.date_range(
            dtarget.time.max().values + pd.Timedelta(1, 'D'),
            dtarget.time.max().values + pd.Timedelta(tchunk, 'D')
        )

        print(f"...from {newdates.min()} to {newdates.max()}.")
        assert len(newdates) == tchunk
        # Extend the current timestep out to size `tchunk`, filling with `Nan`s
        dummy = dnew.reindex({"time": newdates}, fill_value=np.nan)
        assert len(dummy.time) == tchunk
        # Now, append to the existing Zarr data store
        dummy.to_zarr(zarrpath, append_dim="time")

print('Update zarr process complete!')
# for testing, see katrina's local file set up 
