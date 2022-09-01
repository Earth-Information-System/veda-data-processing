import zarr
import numpy as np
import xarray as xr

dataset = "./SPL3SMP.zarr"

# Fix SPL3SMP coordinates
dz = zarr.open(dataset)
east = dz["/easting_m"]
north = dz["/northing_m"]

x_size = 36032.22  # meters
y_size = 36032.22  # meters
ulx = -17367530.45
uly = 7314540.83

north2 = uly - 0.5*y_size - y_size * np.arange(len(east))
north2 = zarr.array(north2, chunks=east.chunks)
for key,value in east.attrs.items():
    north2.attrs[key] = value
dz["/northing_m_correct"] = north2

east2 = ulx + 0.5*x_size + x_size * np.arange(len(north))
east2 = zarr.array(east2, chunks=north.chunks)
for key,value in north.attrs.items():
    east2.attrs[key] = value
dz["/easting_m_correct"] = east2

zarr.consolidate_metadata(dataset)

dat = xr.open_zarr(dataset, consolidated=True)
dat.easting_m.max().values
dat.isel(northing_m = slice(0, 100))

# ulx_gdal = -17385546.560
# ulx = dat.easting_m.min().values
# ulx_gdal - ulx
