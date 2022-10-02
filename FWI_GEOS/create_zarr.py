import argparse
import xarray as xr
import numpy as np
import glob
import pandas as pd
import os
from dask.diagnostics.progress import ProgressBar

# Keep all the dataset in the data directory
files = sorted(list(glob.glob("data/*.nc")))
zar_file = 'FWI.zarr'
chunking = {
    "time": 100,
    "lat": 100,
    "lon": 100
}
data = xr.open_dataset(files[0])
c=0
for file in files:
    date = file[-11:-3]
    pred_files = sorted(list(glob.glob(f"data/{date}00/*.nc")))   
    dat = xr.open_dataset(file)
    dat=dat.squeeze()
    dat =dat.assign_coords({"time": pd.to_datetime(date),
                           "forecast": 0})
    for i,f in enumerate(pred_files):
        tmp= xr.open_dataset(f)
        tmp=tmp.squeeze()
        tmp=tmp.assign_coords({"time": pd.to_datetime(date),
                              "forecast": i+1})
        dat = xr.concat((dat,tmp), dim='forecast')
    
    dat=dat.expand_dims("time")
    data = dat if c == 0 else xr.concat((data,dat), dim="time")
    c=c+1
    
ntime = len(data.time) 
residual_chunks = chunking["time"] - (ntime % chunking["time"]) 
dates = data.time.values
date_max = dates.max()
newdates = pd.date_range(
    dates.max() + pd.Timedelta(1, 'D'), #added one day to max_date
    dates.max() + pd.Timedelta(residual_chunks, 'D') #added 100days to max_date
)
assert len(newdates) == residual_chunks ## check if length of new dates= residual chunks
alldates = np.concatenate((dates, newdates))
assert len(alldates) % chunking["time"] == 0
dat_extended = data.reindex({"time": alldates}, fill_value = np.nan)
assert len(dat_extended.time) % chunking["time"] == 0

dat_extended = dat_extended.chunk(chunking)
with ProgressBar():
    dat_extended.to_zarr(zar_file, mode="w")
