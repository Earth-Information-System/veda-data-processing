# Algorithm for continuously updating a Zarr data store
#
# (1) On creation, fill the last time dimension to the end of the chunk with
# `nan`s.
# 
# (2) On update:
#   (a) Open existing store. Check last timestamp and time chunk size.
#   (b1) If timestamp in file, replace.
#   (b2) If not, create new chunk of same size and replace.

import datetime
import re
import glob
import sys

import xarray as xr
import pandas as pd
import numpy as np

prefix = "OCO2_GEOS_L3CO2_day"
procfile = "processed-files.txt"

zarrpath = f"{prefix}-v2.zarr"

allfiles = sorted(list(glob.glob("./raw-data/*/oco2_GEOS_L3CO2_day_*.nc4")))
# allfiles = allfiles[0:110]
with open(procfile, "r") as f:
    # Remove empty lines
    procfiles = [l for l in f.read().splitlines() if l != '']
newfiles = sorted(set(allfiles) - set(procfiles))

if len(newfiles) == 0:
    sys.exit("No new files to process!")

dnew_all = xr.open_mfdataset(newfiles)
dates = dnew_all.time.values

# For now, do this individually for each timestep. It's very inefficient (lots
# of writes) but safe and conceptually simple, and shouldn't take too long for
# small numbers of files.
# TODO The much faster way to do this is to split up `dnew_all` into existing
# vs. new chunks and then do those separately.
for idt, t in enumerate(dnew_all.time):
    # idt = 1; t = dnew_all.time[idt]

    dnew = dnew_all.isel(time=slice(idt, idt+1))
    print(dnew.time.values[0])
    # Is the current timestamp in the Zarr file?
    # NOTE: Have to re-read `dtarget` each iteration because it may have been
    # expanded by the code below.
    dtarget = xr.open_zarr(zarrpath)
    is_in = dtarget.time.isin(dnew.time)
    n_is_in = is_in.sum()
    if n_is_in > 1:
        raise Exception(f"Found {n_is_in} matching times! Something is wrong with this situation.")
    elif n_is_in == 1:
        # Yes! Where exactly?
        itime = int(np.where(is_in)[0])
        assert dtarget.isel(time=itime).time == dnew.time, "Mismatch in times"
        print(f"Inserting into time location {itime}")
        # Write to that exact location (replace NaNs with new data)
        dnew.to_zarr(zarrpath, region={
            "time": slice(itime, itime+1),
            "lat": slice(0, dnew.sizes["lat"]),
            "lon": slice(0, dnew.sizes["lon"])
        })
    elif n_is_in == 0:
        # No, so we need to extend the time series.
        # For performance, we extend by one full time chunksize (`tchunk`)
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

    ifile = newfiles[idt]
    with open(procfile, "a") as f:
        # Note: Flipping the position of the \n means this adds a blank line
        # between "groups" of processed files. That's a minor feature -- it lets us
        # see when groups of files have been processed.
        f.write("\n" + ifile)

# Test the result
print("Testing new Zarr...")
dtest = xr.open_zarr(zarrpath)
dtest.sel(time = slice("2022-06-01", "2022-06-15"))

print("Done!")
