import xarray as xr
import pandas as pd

# from tqdm.auto import tqdm
from dask.diagnostics import ProgressBar 

from fsspec.implementations.local import LocalFileSystem

from pangeo_forge_recipes.patterns import ConcatDim, FilePattern
from pangeo_forge_recipes.recipes import XarrayZarrRecipe
from pangeo_forge_recipes.storage import FSSpecTarget, StorageConfig

def format_function(time):
    base_pattern = "/css/merra2/MERRA2_all/Y{yyyy}/M{mm}/MERRA2.inst1_2d_asm_Nx.{yyyy}{mm}{dd}.nc4"
    return base_pattern.format(
        yyyy=time.strftime("%Y"),
        mm=time.strftime("%m"),
        dd=time.strftime("%d")
    )

times = pd.date_range("2022-01-01", "2022-02-01", freq="1H")
concat_dim = ConcatDim(name="time", keys=times, nitems_per_file=24)
pattern = FilePattern(format_function, concat_dim)

# # Check the pattern
# i = 0
# for idx, url in pattern.items():
#     print(idx)
#     print(url)
#     i += 1
#     if i > 50:
#         break

storage = StorageConfig(FSSpecTarget(LocalFileSystem(), "test.zarr"),
                        FSSpecTarget(LocalFileSystem(), "test.metadata"))

recipe = XarrayZarrRecipe(
    pattern, 
    storage_config = storage,
    target_chunks = {"lon": 100, "lat": 100, "time": 96},
    cache_inputs = False
)

print('Preparing target...')
recipe.prepare_target()

print('Executing recipe...')
delayed = recipe.to_dask()
with ProgressBar():
    delayed.compute()

# all_chunks = list(recipe.iter_chunks())
# for chunk in tqdm(recipe.iter_chunks(), total=len(all_chunks)):
#     recipe.store_chunk(chunk)

testdat = xr.open_zarr("test.zarr", consolidated=True)
