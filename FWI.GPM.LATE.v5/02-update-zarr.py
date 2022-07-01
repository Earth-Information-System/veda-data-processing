import datetime
import re

import xarray as xr
import pandas as pd
import numpy as np

timevar = "time"
zarrpath = "FWI_GPM_LATE_v5_Daily-test.zarr"
# zarrpath = target
infile = "raw-input/2022/FWI.GPM.LATE.v5.Daily.Default.20220630.nc"

# lastday = xr.open_zarr(zarrpath)[timevar].max().values

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

# Is the current timestamp in the Zarr file?
dtarget = xr.open_zarr(zarrpath)
is_in = dtarget.time.isin(dnew.time)
n_is_in = is_in.sum()
if n_is_in > 1:
    raise Exception(f"Found {n_is_in} matching times! Something is wrong with this situation.")
elif n_is_in == 1:
    # Yes! Where exactly?
    itime = int(np.where(is_in)[0])
    # Write to that exact location (replace NaNs with new data)
    dnew.to_zarr(zarrpath, region={
        "time": slice(itime, itime+1),
        "lat": slice(0, dnew.sizes["lat"]),
        "lon": slice(0, dnew.sizes["lon"])
    })
elif n_is_in == 0:
    # No, so we need to extend the time series.
    # For performance, we extend by one full time chunksize -- `tchunk`
    tchunk = dtarget.chunks["time"][-1]
    newdates = pd.date_range(
        dtarget.time.max().values + pd.Timedelta(1, 'D'),
        dtarget.time.max().values + pd.Timedelta(tchunk, 'D')
    )
    assert len(newdates) == len(tchunk)
    # Extend the current timestep out to size `tchunk`, filling with `Nan`s
    dummy = dnew.reindex({"time": newdates}, fill_value=np.nan)
    assert len(dummy.time) == len(tchunk)
    # Now, append to the existing Zarr data store
    dummy.to_zarr(zarrpath, append_dim="time")

# Test the result
print("Testing new Zarr...")
dtest = xr.open_zarr(zarrpath)
print("Done!")

# Update algorithm
# (1) On creation, fill the last time dimension to the end with `nan`s
# (2) On update:
    # (a) Open existing store. Check last timestamp and time chunk size.
    # (b1) If timestamp in file, replace.
    # (b2) If not, create new chunk of same size and replace.
