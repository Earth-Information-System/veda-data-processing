# LIS Pangeo Forge Recipe

This repo contains a workflow used to convert LIS outputs from NetCDF to Zarr stores. This workflow assumes an AWS-based running environment. Data is read from and written to S3 object storage and inputs are cached temporarily in a local filesystem (an EFS mount in this case) during the conversion process.

The key files are `lis_netcdf_to_zarr.py` and YAML configuration files in `configs/`. The Jupyter Notebooks in `notebooks/` contain the initial prototype for the CLI script `lis_netcdf_to_zarr.py`.

## Environment

The required libraries to install are: `pangeo_forge_recipes`, `fsspec`, `s3fs`, `tqdm`, `numpy`, `xarray`, and `yaml`.

Two `conda` environment files are provided under `env/`:

1. `environment.yml` contains the exact dependencies used on SMCE during the EIS Freshwater pilot project. The environment created by this file may be excessive and may not work properly on a different system. It is provided here to document the exact dependencies used on SMCE because I ran into dependency conflicts between `fsspec` and `s3fs` after letting `conda` resolve the dependencies.
2. `minimal-environment.yml` *should* create a minimal working environment. This file installs only the libraries required in the script and pins `python`, `fsspec`, and `s3fs` to known working versions.

## Setup

1. Create a YAML config file based on examples under `config/`, for example:

```yaml
---

enable_logging:   False
bucket:           'eis-dh-hydro'
input_path:       'LIS_NETCDF/DA_10km/DA_LAI_PMW_HyMAP/SURFACEMODEL/**/LIS_HIST*.nc'
target_path:      'LIS/DA_10km/RECHUNKED/DA_LAI_PMW_HyMAP/SURFACEMODEL/LIS_HIST.d01.zarr'
target_chunks:    {'time': 365}
inputs_per_chunk: 365
nitems_per_file:  1
temp_dir:         '/home/jovyan/efs/tmp'
```

2. Set the values of `input_path` to a pattern matching the input filenames under hosted in the S3 bucket indicated by `bucket`. Set `target_path` to the desired name of the resulting Zarr store.

3. Set `target_chunks` and `inputs_per_chunk` as appropriate (aim for chunks of ~100mb).

4. Create the cache directory indicated by the `temp_dir` path.

## Execution

1. Activate the `lis_pangeo_forge_env` conda environment.

2. Run the script and pass the appropriate config file:

```
python src/lis_netcdf_to_zarr.py path/to/dataset_config.yaml
```

3. Let it run. A progress bar should appear after a few setup messages (and some warnings that can be ignored).

## TODOs

This script currently processes files serially. It should be possible to parallelize this with Dask, but we did not have `pangeo_forge_recipes` installed in the conda environment of our base image on SMCE so it was not available to the Dask workers and would crash.