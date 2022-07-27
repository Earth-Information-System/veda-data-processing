import xarray as xr
import zarr
import pandas as pd
import numpy as np

import glob

zarrpath = "FWI_GPM_LATE_v5_Daily.zarr"
dtarget = xr.open_zarr(zarrpath)

len(dtarget["time"].values)
len(np.unique(dtarget["time"].values))

dt0 = dtarget.isel(time=slice(0, 1))
dt0_2 = dt0
dt0_2 = dt0_2.assign_coords(time=pd.to_datetime(["2014-05-01"]))
dtout = dt0_2[["time", "lat", "lon"]]
dtout.to_zarr(zarrpath, region={
    "time": slice(0, 1),
    "lat": slice(0, dt0_2.sizes["lat"]),
    "lon": slice(0, dt0_2.sizes["lon"])
})


zg = zarr.open(zarrpath)
# zarr.consolidate_metadata(zarrpath)

# NOTE: `NaT` was caused by different in fill_value.
# If fill_value = 0, then value of 0 is assumed to be nan

# Fix the time array
zg_time = zg["time"]
zg_time_vals = zg_time[:]
dict(zg_time.attrs)

nmax = zg_time_vals.max() + 1
nextra = len(zg_time_vals) - nmax
zg_time_fixed = zg_time[:-nextra]
# zg["time"] = zg_time_fixed

dt = xr.open_zarr("FWI_GPM_LATE_v5_Daily-test.zarr")
dt
zt = zarr.open("FWI_GPM_LATE_v5_Daily-test.zarr")
zt_time = zt["time"]
zt_time_vals = zt_time[:]
dict(zt_time.attrs)

# Data variables
zvars = [k for k in list(zg.keys()) if k.startswith("GPM")]

lastchunk = int(nmax / 50)

for zv in zvars:
    # zv = zvars[0]
    # 
    zblobs = sorted(glob.glob(f"{zarrpath}/{zv}/63.*.*"))
    zg_var = zg[zv]
    zg_var.

zg_time = zg["time"]
# Remove the last 150 entries -- they are duplicates
