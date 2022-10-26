import xarray as xr

dtest = xr.open_zarr("FWI.GEOS-5.zarr/")

dtestc = dtest.chunk("auto")

dtest["GEOS-5_FWI"].max().values

tseries = dtest.sel(lat=slice(38,40), lon=slice(-119, 121)).mean(["lat", "lon"])

tseries["GEOS-5_FWI"].sel(forecast=0).values

d0 = xr.open_dataset(files[0])

d0["GEOS-5_FWI"].max().values
