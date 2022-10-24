import re
import os
import glob

import xarray as xr
import numpy as np
import pandas as pd

import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)

from dask.diagnostics.progress import ProgressBar
from tqdm.auto import tqdm

zar_file = 'FWI.GEOS-5.zarr'

# Change to true to force overwriting existing data.
overwrite = False

# Keep all the dataset in the data directory
files = []
# years = range(2017, 2023)
basedir = "/autofs/brewer/rfield1/storage/observations/GFWED/Sipongi/fwiCalcs.GEOS-5/Default/GEOS-5"
assert os.path.exists(basedir)
years = range(2021, 2023)
for y in years:
    assert os.path.exists(f"{basedir}/{y}")
    files += sorted(glob.glob(f"{basedir}/{y}/FWI.GEOS-5.Daily.*.nc"))

for file in tqdm(files):
    # file = files[0]
    date_match = re.match(r".*\.(\d{8})\.nc", os.path.basename(file))
    assert date_match
    date = date_match.group(1)
    pd_date = pd.to_datetime(date)
    dtarget = xr.open_zarr(zar_file)
    # NOTE: This doesn't work for some reason. Turns up true at every timestep,
    # even if they should be all nans. Comment out for now.
    # if not overwrite and dtarget.sel(time=pd_date)["GEOS-5_FWI"].notnull().any().values:
    #     print(f"Data exist for timestep {date}. Skipping this timestep...")
    #     continue
    year = date[0:4]
    datedir = f"{basedir}/{year}/{date}00"
    assert os.path.exists(datedir)
    # Get list of forecast files.
    pred_files = sorted(glob.glob(f"{datedir}/*.nc"))   
    if not len(pred_files):
        print(f"No files found for date {date}. Skipping this timestep.")
        continue
    # Open base date file (forecast = 0)
    dat = xr.open_mfdataset([file] + pred_files, combine="nested", concat_dim="forecast")
    dnew = dat.assign_coords({"forecast": range(0, 1+len(pred_files))})
    
    is_in = dtarget.time.isin(dnew.time)
    n_is_in = is_in.sum()
    if n_is_in > 1:
        raise Exception(f"Found {n_is_in} matching times! Something is wrong with this situation.")
    elif n_is_in == 1:
        # Yes! Where exactly?
        itime = int(np.where(is_in)[0])
        # print(f"Inserting into time location {itime}")
        # Write to that exact location (replace NaNs with new data)
        dnew.to_zarr(zar_file, region={
            "time": slice(itime, itime+1),
            "forecast": slice(0, dnew.sizes["forecast"]),
            "lat": slice(0, dnew.sizes["lat"]),
            "lon": slice(0, dnew.sizes["lon"])
        })
