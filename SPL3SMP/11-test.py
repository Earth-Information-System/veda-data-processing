import xarray as xr

dspath = "s3://veda-data-store-staging/EIS/zarr/SPL3SMP.zarr"
dat = xr.open_zarr(dspath, consolidated=True)
dat.isel(easting_m=300, northing_m=300).soil_moisture.values
dat.sel(datetime="2020-06-15").soil_moisture.values
