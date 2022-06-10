
##################################################

merradir = Path("/css/merra2/MERRA2_all/")
merrafiles = sorted(glob("/css/merra2/MERRA2_all/Y2022/M01/MERRA2.inst1_2d_asm_Nx.????????.nc4"))

merra_pattern = pattern_from_file_sequence(merrafiles, "time", nitems_per_file=24)

merra_recipe = XarrayZarrRecipe(
    file_pattern = merra_pattern,
    storage_config = ..., # TODO
    cache_inputs = False,
    target_chunks = {"lon": 100, "lat": 100, "time": 96}
)

len(list(merra_recipe.iter_inputs()))
len(list(merra_recipe.iter_chunks()))

merrafiles = sorted(glob("/css/merra2/MERRA2_all/Y2022/M01/MERRA2.inst1_2d_asm_Nx.????????.nc4"))
d1 = xr.open_dataset(merrafiles[0])
d1.time

dat = xr.open_mfdataset(merrafiles[0:5])
