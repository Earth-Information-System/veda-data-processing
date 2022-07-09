import xarray as xr
import pandas as pd

# from tqdm.auto import tqdm
from dask.diagnostics.progress import ProgressBar 

from fsspec.implementations.local import LocalFileSystem

from pangeo_forge_recipes.patterns import ConcatDim, FilePattern
from pangeo_forge_recipes.recipes import XarrayZarrRecipe
from pangeo_forge_recipes.storage import FSSpecTarget, StorageConfig

def format_function(time):
    # "./raw-data/2021/oco2_GEOS_L3CO2_day_20211030_B10206Ar.nc4.xml"
    base_pattern = "./raw-data/{yyyy}/oco2_GEOS_L3CO2_day_{yyyy}{mm}{dd}_B10206Ar.nc4"
    return base_pattern.format(
        yyyy=time.strftime("%Y"),
        mm=time.strftime("%m"),
        dd=time.strftime("%d")
    )

prefix = "OCO2_GEOS_L3CO2_day"
outfile = f"{prefix}.zarr"
outmeta = f"{prefix}.metadata"

start_date = "2015-01-01"
# start_date = "2021-01-01"
end_date = "2021-08-01"
chunking = {
    "lon": 100,
    "lat": 100,
    "time": 100
}

times = pd.date_range(start_date, end_date, freq="1D")
concat_dim = ConcatDim(name="time", keys=times, nitems_per_file=1)
pattern = FilePattern(format_function, concat_dim)

# id0, url0 = next(pattern.items())
# f0 = xr.open_dataset(url0)

storage = StorageConfig(FSSpecTarget(LocalFileSystem(), outfile),
                        FSSpecTarget(LocalFileSystem(), outmeta))

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
annual_mean = testdat["XCO2"].groupby("time.year").mean(["lat", "lon"]).values
print(annual_mean)
print("Done!")
