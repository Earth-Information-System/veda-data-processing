import zarr
import xarray as xr
from rechunker import rechunk
from dask.diagnostics import ProgressBar

def main():
    source = "SPL3SMP-unchunked.zarr"
    source_array = xr.open_zarr(source, consolidated=True)
    strategy = {
        "easting_m": 100,
        "northing_m": 100,
        "datetime": 100
    }
    target = "SPL3SMP.zarr"
    temp_store = "rechunk-temp.zarr"
    max_mem = '16GB'

    print("Creating rechunking plan...")
    plan = rechunk(source_array, strategy, max_mem, target, temp_store=temp_store)
    print("Executing rechunking...")
    with ProgressBar():
        plan.execute()

    print("Consolidating metadata...")
    zarr.consolidate_metadata(target)

    print("Done!")

if __name__ == '__main__':
    main()
