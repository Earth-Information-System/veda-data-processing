import xarray as xr
import pandas as pd

from dask.diagnostics.progress import ProgressBar 

from fsspec.implementations.local import LocalFileSystem

from pangeo_forge_recipes.patterns import ConcatDim, FilePattern
from pangeo_forge_recipes.recipes import XarrayZarrRecipe
from pangeo_forge_recipes.storage import FSSpecTarget, StorageConfig

prefix = "FWI_GPM_LATE_v5_Daily-test"
chunking = {"lon": 100, "lat": 100, "time": 100}
start_date = "2022-01-01"
end_date = "2022-02-01"
frequency = "1D"

def format_function(time):
    # "raw-input/2022/FWI.GPM.LATE.v5.Daily.Default.20220101.nc"
    base_pattern = "raw-input/{yyyy}/FWI.GPM.LATE.v5.Daily.Default.{yyyy}{mm}{dd}.nc"
    return base_pattern.format(
        yyyy=time.strftime("%Y"),
        mm=time.strftime("%m"),
        dd=time.strftime("%d")
    )

outfile = f"{prefix}.zarr"
outmeta = f"{prefix}.metadata"

times = pd.date_range(start_date, end_date, freq=frequency)
concat_dim = ConcatDim(name="time", keys=times, nitems_per_file=1)
pattern = FilePattern(format_function, concat_dim, file_type="netcdf3")

# id0, f0 = next(pattern.items())
# d0 = xr.open_dataset(f0)
# for idx, fx in pattern.items():
#     dx = xr.open_dataset(fx)

storage = StorageConfig(
    FSSpecTarget(LocalFileSystem(), outfile),
    FSSpecTarget(LocalFileSystem(), outmeta)
)

recipe = XarrayZarrRecipe(
    pattern, 
    storage_config = storage,
    target_chunks = chunking,
    cache_inputs = False
)

delayed = recipe.to_dask()
print('Executing recipe...')
with ProgressBar():
    try:
        delayed.compute()
    except RuntimeError as err:
        print(f"Hit error {type(err)=}, {err=}")
        print("Trying again")
        delayed.compute()

print("Testing ability to open dataset.")
testdat = xr.open_zarr(outfile, consolidated=True)
