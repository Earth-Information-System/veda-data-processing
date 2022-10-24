import xarray as xr
import numpy as np
import glob
import pandas as pd
import os
import re

import dask
from dask.diagnostics.progress import ProgressBar
from tqdm.auto import tqdm

import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)

# Keep all the dataset in the data directory
# files = sorted(list(glob.glob("~/Sipongi/fwiCalcs.GEOS-5/Default/GEOS-5/**/FWI.GEOS-5.Daily.Default.*.nc")))
files = []
# years = range(2017, 2023)
years = range(2021, 2023)
for y in years:
    files += sorted(glob.glob(f"/autofs/brewer/rfield1/storage/observations/GFWED/Sipongi/fwiCalcs.GEOS-5/Default/GEOS-5/{y}/FWI.GEOS-5.Daily.*.nc"))

zar_file = 'FWI.GEOS-5.zarr'
chunking = {
    "time": 100,
    "lat": 100,
    "lon": 100
}

def parse_date(fname):
    date_match = re.match(r".*\.(\d{8})\.nc", os.path.basename(fname))
    assert date_match
    date_str = date_match.group(1)
    date = pd.to_datetime(date_str)
    return date

dates = [parse_date(f) for f in files]

d0 = xr.open_dataset(files[0])
blank = dask.array.empty((9, len(dates), d0.lat.shape[0], d0.lon.shape[0]),
                         chunks=(9, 100, 100, 100),
                         dtype=np.float32)
dims = ["forecast", "time", "lat", "lon"]

# Pre-allocate dataset
complete_dat = xr.Dataset(
    data_vars = {
        "GEOS-5_DC" : (dims, blank),
        "GEOS-5_DMC" : (dims, blank),
        "GEOS-5_FFMC" : (dims, blank),
        "GEOS-5_ISI" : (dims, blank),
        "GEOS-5_BUI" : (dims, blank),
        "GEOS-5_FWI" : (dims, blank),
        "GEOS-5_DSR" : (dims, blank)
    },
    coords = {
        "lat": d0.lat,
        "lon": d0.lon,
        "time": dates,
        "forecast": list(range(0, 9))
    },
    attrs = {**d0.attrs}
)

# Extend
ntime = len(complete_dat.time) 
residual_chunks = chunking["time"] - (ntime % chunking["time"]) 
dates = complete_dat.time.values
date_max = dates.max()
newdates = pd.date_range(
    dates.max() + pd.Timedelta(1, 'D'), #added one day to max_date
    dates.max() + pd.Timedelta(residual_chunks, 'D') #added 100days to max_date
)
assert len(newdates) == residual_chunks ## check if length of new dates= residual chunks
alldates = np.concatenate((dates, newdates))
assert len(alldates) % chunking["time"] == 0
dat_extended = complete_dat.reindex({"time": alldates}, fill_value = np.nan)
assert len(dat_extended.time) % chunking["time"] == 0

with ProgressBar():
    dat_extended.to_zarr('FWI.GEOS-5.zarr', mode="w")

# Now, populate this Zarr one date at a time...

for file in tqdm(files):
    # date = file[-11:-3]
    date_match = re.match(r".*\.(\d{8})\.nc", os.path.basename(file))
    assert date_match
    date = date_match.group(1)
    pred_files = sorted(list(glob.glob(f"data/{date}00/*.nc")))   
    dat = xr.open_dataset(file)
    dat = dat.squeeze()
    dat = dat.assign_coords({"time": pd.to_datetime(date),
                             "forecast": 0})
    for i,f in enumerate(pred_files):
        tmp = xr.open_dataset(f)
        tmp = tmp.squeeze()
        tmp = tmp.assign_coords({"time": pd.to_datetime(date),
                                 "forecast": i+1})
        dat = xr.concat((dat,tmp), dim='forecast')
    
    dnew = dat.expand_dims("time")
    # dat.to_zarr("FWI.GEOS-5.zarr",)
    # datalist.append(xr.concat((data, dat), dim="time"))
