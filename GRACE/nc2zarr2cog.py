import argparse
from os import makedirs
from os.path import exists
import xarray as xr
import numpy as np
import rioxarray
import pandas as pd


parser = argparse.ArgumentParser(description='Convert netcdf dataset to zarr and to COG')
parser.add_argument('--dataset', type=str, help='Input Zarr dataset to convert')
parser.add_argument('--nc_dataset', type=str, help='Input netcdf dataset to convert')
parser.add_argument('--outdir', type=str, help='Output directory for COGs')
parser.add_argument('--prefix', type=str, help='Prefix for output file names')
parser.add_argument('--crs', type=str, help='Dataset CRS string (e.g., "epsg:4326")')
parser.add_argument('--x_dim', type=str, help='Name of X dimension')
parser.add_argument('--y_dim', type=str, help='Name of Y dimension')
parser.add_argument("--timevar", type=str, help='Name of time dimension')
parser.add_argument("--timeunit", type=str, help='Time unit. See `numpy.datetime_as_string` argument `unit`. (D=daily)')


args = parser.parse_args()
# args = parser.parse_args(['--outdir=GSFC-cog', '--nc_dataset=gsfc.glb_.200204_202112_rl06v2.0_obp-ice6gd_halfdegree.nc', '--dataset=GSFC.zarr', '--crs=EPSG:4326', '--x_dim=lon', '--y_dim=lat', '--timevar=time', '--timeunit=D', '--prefix=GSFC'])
# python nc2zarr2cog.py --outdir GSFC-cog --nc_dataset gsfc.glb_.200204_202112_rl06v2.0_obp-ice6gd_halfdegree.nc --dataset GSFC.zarr --crs EPSG:4326 --x_dim lon --y_dim lat --timevar time --timeunit D --prefix GSFC

outdir = args.outdir
dataset = args.dataset
nc_dataset= args.nc_dataset
crs = args.crs
y_dim = args.y_dim
x_dim = args.x_dim
timevar = args.timevar
prefix = args.prefix
timeunit = args.timeunit

ds = xr.open_dataset(nc_dataset)
ds_chunked=ds.chunk({
    "time": 100,
    "lat": 100,
    "lon": 100
})
ds_chunked.to_zarr(f"{prefix}.zarr", mode="w")
dat = xr.open_zarr(dataset, consolidated=True).transpose(y_dim, x_dim, ...)
dat.rio.set_spatial_dims(x_dim=x_dim, y_dim=y_dim, inplace=True)

makedirs(outdir, exist_ok=True)

for var in dat.keys():
    dat_v = dat[var]
    print(f"for variable {var}")
    if timevar not in dat_v.dims:
        print(f"Skipping variable `{var}` without time dimension.")
        continue
    for t in dat_v[timevar]:
        datestring = np.datetime_as_string(t, unit="M")
        print(f"{var}: {datestring}")
        x = dat_v.sel({timevar: t})
        if x_dim not in dat_v.coords:
            continue
        x['lon']= np.where(x.lon.values>180,x.lon.values-360,x.lon.values)
        x=x.sortby(['lon'])
        x=x.sortby("lat",ascending=False)
        outfile = f"{outdir}/{prefix}-{var}-{datestring}.tif"
        if exists(outfile):
            continue
        x.rename({y_dim: 'y', x_dim: 'x'}).rio.to_raster(outfile, driver='COG')