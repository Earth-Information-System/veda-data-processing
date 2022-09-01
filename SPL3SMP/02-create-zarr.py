import os
import re
import datetime

import numpy as np
import pandas as pd
import xarray as xr

from dask.diagnostics.progress import ProgressBar

def main():

    short_name = "SPL3SMP"
    zarr_path = f"{short_name}-test.zarr"
    procfile = "processed-files.txt"

    datafiles_all = sorted(os.listdir(short_name))
    datafiles = [short_name + '/' + f for f in datafiles_all if f.endswith("h5")]

    # len(datafiles)
    # datafiles = datafiles[0:110]

    chunking = {
        "easting_m": 100,
        "northing_m": 100,
        "time": 100
    }

    print("Reading AM files...")
    dat_am = xr.open_mfdataset(datafiles,
                               preprocess=mypreproc,
                               engine='h5netcdf',
                               phony_dims='sort',
                               combine='by_coords',
                               group="Soil_Moisture_Retrieval_Data_AM")

    print("Reading PM files...")
    dat_pm = xr.open_mfdataset(datafiles,
                               preprocess=mypreproc,
                               engine='h5netcdf',
                               phony_dims='sort',
                               combine='by_coords',
                               group="Soil_Moisture_Retrieval_Data_PM")

    print("Combining datasets...")
    dat_both = xr.merge((dat_am, dat_pm))

    print("Extending time series for even chunking...")
    ntime = len(dat_both.time)
    residual_chunks = chunking["time"] - (ntime % chunking["time"])
    lastdate = dat_both.time.max().values
    newdates = pd.date_range(
        lastdate + pd.Timedelta(1, 'D'),
        lastdate + pd.Timedelta(residual_chunks, 'D')
    )
    assert len(newdates) == residual_chunks
    alldates = np.concatenate((dat_both.time.values, newdates))
    assert len(alldates) % chunking["time"] == 0
    dat_extended = dat_both.reindex({"time": alldates}, fill_value = np.nan)
    assert len(dat_extended.time) % chunking["time"] == 0
    dat_extended = dat_extended.chunk(chunking)

    print("Writing Zarr output...")
    dat_final = dat_extended.chunk(chunking)
    with ProgressBar():
        dat_final.to_zarr(zarr_path, consolidated=True, mode='w')

    print(f"Writing file list to {procfile} ...")
    with open(procfile, "w") as f:
        f.writelines(file + "\n" for file in datafiles)

    print("Testing ability to open Zarr...")
    dtest = xr.open_zarr(zarr_path)
    dtest.sel(time = slice("2018-01-01", "2018-02-15"))
    print("Done!")

# Preprocessing function -- takes dataset as an argument
def mypreproc(d_am_raw):
    fname = os.path.basename(d_am_raw.encoding["source"])
    fdate_str = re.search(r'(?<=_SM_P_)\d+(?=_)', fname).group()
    fdate = datetime.datetime.strptime(fdate_str, "%Y%m%d")
    dim_names = sorted(d_am_raw.dims.keys())
    is_am = dim_names[0] == "phony_dim_0"
    d_am = d_am_raw
    d_am = d_am.rename_dims({
        dim_names[0]: "northing_m",
        dim_names[1]: "easting_m",
        dim_names[2]: "ranked_land_cover"
    })

    # EASE grid information
    x_size = 36032.22  # meters
    y_size = 36032.22  # meters
    ulx = -17367530.45
    uly = 7314540.83
    # ncols = 964
    # nrows = 406
    # ease_epsg = "epsg:6933"

    # NOTE: Upper-left to lower-right; therefore, x increases but y decreases
    # NOTE: ulx,uly are the corners of the pixel. To get to the center, add half the size
    d_am = d_am.assign_coords({
        "easting_m": ulx + 0.5*x_size + x_size * d_am.easting_m,
        "northing_m": uly - 0.5*y_size - y_size * d_am.northing_m,
        "datetime": fdate
    })
    d_am = d_am.expand_dims('datetime', axis=2)
    use_vars = ["albedo",
                "bulk_density",
                "clay_fraction",
                "freeze_thaw_fraction",
                "grid_surface_status",
                "radar_water_body_fraction",
                "retrieval_qual_flag",
                "roughness_coefficient",
                "soil_moisture",
                "soil_moisture_error",
                "static_water_body_fraction",
                "surface_flag",
                "surface_temperature"]
    if not is_am:
        use_vars = [v + "_pm" for v in use_vars]
    return d_am[use_vars]

if __name__ == '__main__':
    main()
