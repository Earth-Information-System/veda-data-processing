import kerchunk.hdf
import kerchunk.combine
import fsspec
import ujson
import xarray as xr

fs = fsspec.filesystem('https')

ncfiles = sorted(fs.glob("raw-data/*.nc4"))

# url_dir = "https://portal.nccs.nasa.gov/datashare/gmao/geos-fp/das/Y2022/M09/D01/"
# fname = "GEOS.fp.asm.tavg1_2d_flx_Nx.20220901_0030.V01.nc4"
# url = f"{url_dir}/{fname}"

url = "https://portal.nccs.nasa.gov/datashare/gmao/geos-fp/das/Y2022/M09/D03/GEOS.fp.asm.inst1_2d_lfo_Nx.20220903_0000.V01.nc4"

fs.exists(url)

with fs.open(url, "rb") as f:
    chunks = kerchunk.hdf.SingleHdf5ToZarr(f, url)
    outfile = "http-01.json"
    with fsspec.open(outfile, "wb") as of:
        of.write(ujson.dumps(chunks.translate()).encode())

################################################################################
dtest = xr.open_dataset("reference://", engine="zarr", backend_kwargs={
    "consolidated": False,
    "storage_options": {"fo": "http-01.json"}
})

################################################################################

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
