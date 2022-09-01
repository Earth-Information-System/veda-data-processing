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
import importlib
import os

import xarray as xr
import pandas as pd
import numpy as np

czmod = importlib.import_module("02-create-zarr")
importlib.reload(czmod)
mypreproc = czmod.mypreproc

short_name = "SPL3SMP"
procfile = "processed-files.txt"
zarrpath = f"{short_name}.zarr/"

datafiles_all = sorted(os.listdir(short_name))
datafiles = [short_name + '/' + f for f in datafiles_all if f.endswith("h5")]

with open(procfile, "r") as f:
    # Remove empty lines
    procfiles = [l for l in f.read().splitlines() if l != '']
newfiles = sorted(set(datafiles) - set(procfiles))
# newfiles = newfiles[0:10]

if len(newfiles) == 0:
    sys.exit("No new files to process!")

print(f"First processed file is: {procfiles[1]}")
print(f"Last processed file is: {procfiles[-1]}")
print(f"First new file is: {newfiles[0]}")
print(f"Last new file is: {newfiles[-1]}")
print(f"Files to process: {len(newfiles)}")

print("Reading AM files...")
dnew_am = xr.open_mfdataset(newfiles,
                           preprocess=mypreproc,
                           engine='h5netcdf',
                           phony_dims='sort',
                           combine='by_coords',
                           group="Soil_Moisture_Retrieval_Data_AM")

print("Reading PM files...")
dnew_pm = xr.open_mfdataset(newfiles,
                           preprocess=mypreproc,
                           engine='h5netcdf',
                           phony_dims='sort',
                           combine='by_coords',
                           group="Soil_Moisture_Retrieval_Data_PM")

print("Combining datasets...")
dnew_all = xr.merge((dnew_am, dnew_pm))

assert len(dnew_all.datetime) == len(newfiles), "Mismatch between number of files and number of timesteps"

# For now, do this individually for each timestep. It's very inefficient (lots
# of writes) but safe and conceptually simple, and shouldn't take too long for
# small numbers of files.
# TODO The much faster way to do this is to split up `dnew_all` into existing
# vs. new chunks and then do those separately.
for idt, t in enumerate(dnew_all.datetime):
    # idt = 0; t = dnew_all.datetime[idt]
    dnew = dnew_all.isel(datetime=slice(idt, idt+1))
    print(dnew.datetime.values[0])
    # Is the current timestamp in the Zarr file?
    dtarget = xr.open_zarr(zarrpath)
    is_in = dtarget.datetime.isin(dnew.datetime)
    n_is_in = is_in.sum()
    if n_is_in > 1:
        raise Exception(f"Found {n_is_in} matching times! Something is wrong with this situation.")
    elif n_is_in == 1:
        # Yes! Where exactly?
        itime = int(np.where(is_in)[0])
        print(f"Inserting into time location {itime}")
        # Write to that exact location (replace NaNs with new data)
        dnew.to_zarr(zarrpath, region={
            "datetime": slice(itime, itime+1),
            "northing_m": slice(0, dnew.sizes["northing_m"]),
            "easting_m": slice(0, dnew.sizes["easting_m"])
        })
    elif n_is_in == 0:
        # No, so we need to extend the time series.
        # For performance, we extend by one full time chunksize (`tchunk`)
        tchunk = dtarget.chunks["datetime"][1]
        print(f"Creating new time chunk of size {tchunk}...")
        newdates = pd.date_range(
            dtarget.datetime.max().values + pd.Timedelta(1, 'D'),
            dtarget.datetime.max().values + pd.Timedelta(tchunk, 'D')
        )
        print(f"...from {newdates.min()} to {newdates.max()}.")
        assert len(newdates) == tchunk
        # Extend the current timestep out to size `tchunk`, filling with `Nan`s
        dummy = dnew.reindex({"datetime": newdates}, fill_value=np.nan)
        assert len(dummy.datetime) == tchunk
        # Now, append to the existing Zarr data store
        dummy.to_zarr(zarrpath, append_dim="datetime")
    ifile = newfiles[idt]
    with open(procfile, "a") as f:
        # Note: Flipping the position of the \n means this adds a blank line
        # between "groups" of processed files. That's a minor feature -- it lets us
        # see when groups of files have been processed.
        f.write("\n" + ifile)

# Test the result
print("Testing new Zarr...")
dtest = xr.open_zarr(zarrpath, consolidated = True)
dtest.sel(datetime = slice("2018-01-01", "2018-02-15"))
print("Done!")
