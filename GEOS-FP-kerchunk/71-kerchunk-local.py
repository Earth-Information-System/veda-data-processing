import kerchunk.hdf
import kerchunk.combine
import fsspec
import ujson
import xarray as xr

fs = fsspec.filesystem('file')

ncfiles = sorted(fs.glob("raw-data/*.nc4"))

xr.open_mfdataset(ncfiles)
# <xarray.Dataset>
# Dimensions:   (lon: 1152, lat: 721, time: 7)
# Coordinates:
#   * lon       (lon) float64 -180.0 -179.7 -179.4 -179.1 ... 179.1 179.4 179.7
#   * lat       (lat) float64 -90.0 -89.75 -89.5 -89.25 ... 89.25 89.5 89.75 90.0
#   * time      (time) datetime64[ns] 2022-09-01T00:30:00 ... 2022-09-01T06:30:00
# Data variables: (12/43)
#     BSTAR     (time, lat, lon) float32 dask.array<chunksize=(1, 721, 1152), meta=np.ndarray>
#     CDH       (time, lat, lon) float32 dask.array<chunksize=(1, 721, 1152), meta=np.ndarray>
#     CDM       (time, lat, lon) float32 dask.array<chunksize=(1, 721, 1152), meta=np.ndarray>
# ...

for ncf in ncfiles:
    with fs.open(ncf, "rb") as f:
        chunks = kerchunk.hdf.SingleHdf5ToZarr(f, ncf)
        outfile = ncf.replace(".nc4", ".json")
        with fs.open(outfile, "wb") as of:
            of.write(ujson.dumps(chunks.translate()).encode())

jfiles = sorted(fs.glob("raw-data/*.json"))

xr.open_dataset(
    "reference://", engine = "zarr",
    backend_kwargs = {
        "storage_options": {"fo": jfiles[0]},
        "consolidated": False
    }
)
# <xarray.Dataset>
# Dimensions:   (time: 1, lat: 721, lon: 1152)
# Coordinates:
#   * lat       (lat) float64 -90.0 -89.75 -89.5 -89.25 ... 89.25 89.5 89.75 90.0
#   * lon       (lon) float64 -180.0 -179.7 -179.4 -179.1 ... 179.1 179.4 179.7
#   * time      (time) datetime64[ns] 2022-09-01T00:30:00
# Data variables: (12/43)
#     BSTAR     (time, lat, lon) float32 ...
#     CDH       (time, lat, lon) float32 ...
#     CDM       (time, lat, lon) float32 ...
# ...

multi = kerchunk.combine.MultiZarrToZarr(
    jfiles,
    concat_dims = ["time"],
    identical_dims = ["lat", "lon"],
    coo_map = {"time": "cf:time"}
)
combined = "test.json"
with fs.open(combined, "wb") as of:
    of.write(ujson.dumps(multi.translate()).encode())

# Try to read the result
dmulti = xr.open_dataset(
    "reference://", engine = "zarr",
    backend_kwargs = {
        "storage_options": { 
            "fo": combined
        },
        "consolidated": False
    }
)

dmulti

# <xarray.Dataset>
# Dimensions:   (time: 1, lat: 721, lon: 1152)
# Coordinates:
#   * lat       (lat) float64 -90.0 -89.75 -89.5 -89.25 ... 89.25 89.5 89.75 90.0
#   * lon       (lon) float64 -180.0 -179.7 -179.4 -179.1 ... 179.1 179.4 179.7
#   * time      (time) datetime64[ns] NaT
# Data variables: (12/43)
#     BSTAR     (time, lat, lon) float32 ...
#     CDH       (time, lat, lon) float32 ...
#     CDM       (time, lat, lon) float32 ...
# ...
