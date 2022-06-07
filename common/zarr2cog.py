import argparse
from os import makedirs
from os.path import exists

import xarray as xr
import numpy as np
import rioxarray

parser = argparse.ArgumentParser(description='Convert Zarr dataset to COG')
parser.add_argument('--dataset', type=str, help='Input Zarr dataset to convert')
parser.add_argument('--outdir', type=str, help='Output directory for COGs')
parser.add_argument('--crs', type=str, help='Dataset CRS string (e.g., "epsg:4326")')
parser.add_argument('--x_dim', type=str, help='Name of X dimension')
parser.add_argument('--y_dim', type=str, help='Name of Y dimension')
parser.add_argument("--timevar", type=str, help='Name of time dimension')

args = parser.parse_args()
# args = parser.parse_args(['--outdir=../SPL3SMP/SPL3SMP-cog/', '--dataset=../SPL3SMP/SPL3SMP.zarr/', '--crs=EPSG:6933', '--x_dim=easting_m', '--y_dim=northing_m', '--timevar=datetime'])

outdir = args.outdir
dataset = args.dataset
crs = args.crs
y_dim = args.y_dim
x_dim = args.x_dim
timevar = args.timevar

dat = xr.open_zarr(dataset, consolidated=True).transpose(y_dim, x_dim, ...)
dat.rio.write_crs(crs, inplace=True)
dat.rio.set_spatial_dims(x_dim=x_dim, y_dim=y_dim, inplace=True)

makedirs(outdir, exist_ok=True)

cogvars = sorted(dat.keys())

for var in cogvars:
    # var = list(cogvars)[0]
    dat_v = dat[var]
    if timevar not in dat_v:
        print(f"Skipping variable `{var}` without time dimension.")
        continue
    vardir = f"{outdir}/{var}"
    makedirs(vardir, exist_ok=True)
    for t in dat_v[timevar]:
        # t = dat[timevar][0]
        datestring = np.datetime_as_string(t, unit="D")
        print(f"{var}: {datestring}" + " "*40, end="\r")
        outfile = f"{vardir}/SPL3SMP-{var}-{datestring}.tif"
        if exists(outfile):
            continue
        dat_vt = dat_v.sel({timevar: t})
        dat_vt.\
            rename({y_dim: 'y', x_dim: 'x'}).\
            rio.to_raster(outfile)
