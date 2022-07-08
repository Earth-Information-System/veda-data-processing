import xarray as xr
import rasterio
import rioxarray

from os import makedirs
from os.path import exists

outdir_root = "processed-output"

# Monthly output
outdir_m = f"{outdir_root}/monthly"
makedirs(outdir_m, exist_ok=True)
dmonth = xr.open_dataset("raw-input/Gridded Methane Data/GEPA_Monthly.nc")
dmonth = dmonth.rio.write_crs("epsg:4326")
dmonth = dmonth.rio.set_spatial_dims(x_dim="lon", y_dim="lat")

dmonth_vars = list(dmonth.keys())

for var in dmonth_vars:
    # var = dmonth_vars[0]
    outdir_mv = f"{outdir_m}/{var}"
    makedirs(outdir_mv, exist_ok=True)
    for itime, t in enumerate(dmonth.time):
        # itime = 0; t = dmonth.time[itime]
        tstring = f"{t.values:02.0f}"
        dtv = dmonth.isel(time = itime)[var]
        print(f"{var}: Month {tstring}".ljust(60), end="\r")
        outfile = f"{outdir_mv}/EPA-monthly-{var}-{tstring}.tif"
        if exists(outfile):
            continue
        dtv.rio.to_raster(outfile, driver="COG")

# Same thing for daily
outdir_d = f"{outdir_root}/daily"
makedirs(outdir_d, exist_ok=True)
dday = xr.open_dataset("raw-input/Gridded Methane Data/GEPA_Daily.nc")
dday = dday.rio.write_crs("epsg:4326")
dday = dday.rio.set_spatial_dims(x_dim="lon", y_dim="lat")

dday_vars = list(dday.keys())
for var in dday_vars:
    # var = dday_vars[0]
    outdir_dv = f"{outdir_d}/{var}"
    makedirs(outdir_dv, exist_ok=True)
    for itime, t in enumerate(dday.time):
        # itime = 0; t = dday.time[itime]
        tstring = f"{t.values:03.0f}"
        dtv = dday.isel(time = itime)[var]
        print(f"{var}: Day {tstring}".ljust(60), end="\r")
        outfile = f"{outdir_dv}/EPA-daily-{var}-{tstring}.tif"
        if exists(outfile):
            continue
        dtv.rio.to_raster(outfile, driver="COG")

# Annual emissions are just one timestep
outdir_a = f"{outdir_root}/annual"
makedirs(outdir_a, exist_ok=True)
dann = xr.open_dataset("raw-input/Gridded Methane Data/GEPA_Annual.nc")
dann = dann.rio.write_crs("epsg:4326")
dann = dann.rio.set_spatial_dims(x_dim="lon", y_dim="lat")

dann_vars = list(dann.keys())
for var in dann_vars:
    # var = dann_vars[0]
    dtv = dann[var]
    print(f"{var}".ljust(60), end="\r")
    outfile = f"{outdir_a}/EPA-annual-{var}.tif"
    if exists(outfile):
        continue
    dtv.rio.to_raster(outfile, driver="COG")
print("")
print("Done!")
