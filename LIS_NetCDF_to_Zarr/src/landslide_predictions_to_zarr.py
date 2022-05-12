import argparse, yaml
import s3fs, tempfile, logging

from tqdm import tqdm

from fsspec.implementations.local import LocalFileSystem

from pangeo_forge_recipes.patterns import pattern_from_file_sequence
from pangeo_forge_recipes.recipes import XarrayZarrRecipe
from pangeo_forge_recipes.storage import FSSpecTarget, CacheFSSpecTarget, MetadataTarget

import numpy as np
import xarray as xr

if __name__ == '__main__':

    s3 = s3fs.S3FileSystem(anon=False)

    protocol = 's3://'

    #### parse CLI arguments ####

    # create argument parser
    parser = argparse.ArgumentParser(description="Convert LIS output from NetCDF to Zarr on S3.")

    # add YAML config argument
    parser.add_argument('config', metavar='CFG', type=str,
                       help='Path to a YAML configuration file')

    # parse args
    args = parser.parse_args()

    # TODO: add validation for file and contents

    # parse YAML config file
    with open(args.config) as f:
        config_dict = yaml.safe_load(f)

    # convert YAML key,value pairs to variables
    globals().update(config_dict)

    # turn on pangeo_forge_recipes logging?
    if enable_logging:
        logger = logging.getLogger('pangeo_forge_recipes')
        formatter = logging.Formatter('%(name)s:%(levelname)s - %(message)s')
        handler = logging.StreamHandler()
        handler.setLevel(logging.INFO)
        handler.setFormatter(formatter)
        logger.setLevel(logging.INFO)
        logger.addHandler(handler)

        
    def build_url(path):
        
        return protocol + '/'.join([bucket, path])

    input_url = build_url(input_path)
    
    # create input URLs by globbing input URL pattern
    input_urls = [protocol + s for s in s3.glob(input_url)]

    # define recipe file pattern
    pattern = pattern_from_file_sequence(input_urls,                      # source paths
                                         'time',                          # concat dimension
                                         nitems_per_file=nitems_per_file) # items per file


    # define local FS object
    fs_local = LocalFileSystem()

    # create target FS object
    target_url = build_url(target_path)
    target = FSSpecTarget(fs=s3, root_path=target_url)

    # create metadata target path and FS object
    meta_dir = tempfile.TemporaryDirectory(dir=temp_dir)
    meta_store = MetadataTarget(fs_local, meta_dir.name)

    print('Creating recipe...')
    
    # define the Pangeo Forge recipe
    recipe = XarrayZarrRecipe(pattern,                           # file URL pattern
                              inputs_per_chunk=inputs_per_chunk, # input files per chunk
                              cache_inputs=False,                # read inputs directly from S3
                              target_chunks=target_chunks)       # set chunking scheme for output

    # add metadata and target paths to recipe
    recipe.metadata_cache = meta_store
    recipe.target = target

    # get list of all chunks
    all_chunks = list(recipe.iter_chunks())

    # prepare the target (create empty Zarr store)
    print('Preparing target...')
    recipe.prepare_target()

    # execute the recipe (conversion)
    print('Executing recipe...')
    for chunk in tqdm(recipe.iter_chunks(), total=len(all_chunks)):
        recipe.store_chunk(chunk)

    # consolidate metadata
    print('Finalizing recipe...')
    recipe.finalize_target()

    print('Recipe executed!')
