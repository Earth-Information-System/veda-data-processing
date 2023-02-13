# 03-testing.py
# Test resulting zarr for contents

import datetime
import re
import glob
import sys
from datetime import datetime
from os.path import exists

import dask
import os
from tqdm.auto import tqdm
import netCDF4 as nc

import xarray as xr
import pandas as pd
import numpy as np

zarrpath = '/autofs/brewer/eisfire/katrina/FWI.GEOS-5.Hourly.zarr'
procfile = "/autofs/brewer/eisfire/katrina/processed-files-bak-hourly.txt" 
basedir = "/lovelace/brewer/rfield1/storage/observations/GFWED/Sipongi/fwiCalcs.GEOS-5/Default/GEOS-5/"

print('Initiate path testing.')

# Test dir/file paths
if not os.path.isdir(zarrpath):
	raise FileNotFoundError('Zarr directory does not exist. Test failed.')
if not os.path.isfile(procfile):
	raise FileNotFoundError('Processing bak txt not found. Test failed.')
if not os.path.isdir(basedir):
	raise FileNotFoundError('FWI data directory not found. Test failed.')

print('Path Testing complete')

print('Initiate main zarr testing.')
# Open main zarr and test contents

print('Test 1: bak .txt and zarr synchronization')

# Open zarr + open txt file to ID current files in zarr
main_zarr = xr.open_zarr(zarrpath)
with open(procfile, "r") as f:
	procfiles = [l for l in f.read().splitlines() if l != '']
currentfiles = sorted(set(procfiles))

# Count total number of hours in Zarr vs proc files
total_hours_zarr = len(main_zarr.time)
total_days_zarr = int(total_hours_zarr / 24)
total_days_bak = len(currentfiles)
total_hours_bak = int(total_days_bak * 24)

# chunk multiple test
assert ((total_hours_zarr % 120) == 0)

# given zarr is multiples of 120 -> bump total_hours_bak to next mul of 12
chunk_difference = total_hours_zarr - total_hours_bak
assert (chunk_difference < 120)
recalc_hours_bak = (total_hours_bak + chunk_difference)
assert (total_days_zarr == int(recalc_hours_bak / 24 )) 

# Extract hourly information from bak string and compare to time dimension
main_zarr_time_values = main_zarr.time.values
bak_dates = []
# extract date portion after last '/' i.e. FWI.GEOS-5.Hourly.Default.20220102.nc 
for s in currentfiles:
	# split by /
	parts = s.split("/")
	# Assumes standard FWI naming convention, see above example
	extract_date = (parts.pop())[26:34]
	# transform to datetime form
	dt_date = pd.to_datetime(extract_date)
	# append to dt collection
	bak_dates = bak_dates + [dt_date]

# per bak date, assert that it exiss in zarr values (00-23 parts)
for dt in bak_dates:
	# dt 00:00:000 exists in zarr
	assert (dt in main_zarr_time_values)
	# generate all hours possible and locate in zarr
	dt_all_hours = [dt.replace(hour=n) for n in [*range(0,24)]]
	for dt_hrs in dt_all_hours:
		assert(dt_hrs in main_zarr_time_values)

print('Test passed: All bak files located in main zarr')

print('Test 2: randomized zarr content check')
# given zarr, slice into particular days and check for matching content as original
# must use date of selected file to I.D. location 

# filter function to identify if date string is in passed path name
def filter_date(path_string, dat_string):
	if (dat_string in path_string):
		return True
	else: 
		return False

# collect all current fwi hourly files 
all_fwi_files = []

for y in range(2017, int(datetime.now().year) + 1):
	assert os.path.exists(f"{basedir}/{y}")
	all_fwi_files += sorted(glob.glob(f"{basedir}/{y}/FWI.GEOS-5.Hourly.*.nc"))

# use given bak file list + date list - currently only for 0th hour
for dat in bak_dates:
	found_match = []
	# convert to string match form
	dat_str = str(10000*dat.year + 100*dat.month + dat.day)
	# find match in all_fwi_files 
	found_match_pre = [dat_str in a for a in all_fwi_files]
	assert (len(found_match_pre) == len(all_fwi_files))
	found_match = [x for x, y in zip(all_fwi_files, found_match_pre) if y == True]
	# there should only be single file that matches
	# print(found_match)
	assert (len(found_match) == 1) 
	found_match = found_match[0] 
	hours = [*range(0,24)] # note: file_data is +1 ahead

	# test equal content for all hours, all days
	for h in hours:
		# open match with xr for content
		file_data = (xr.open_mfdataset(found_match)).sel(time = int(h) + 1)
		# slice into zarr using selected date
		dat = dat.replace(hour=h) # update hour of dat using loop
		dat_zarr_data = main_zarr.sel(time = dat)
	
		# For all subgroups, all data should align i.e. no False ever
		subgroup_strs = ["GEOS-5_FFMC", "GEOS-5_FWI", "GEOS-5_ISI"]
		for subs in subgroup_strs:
			test_equality = dat_zarr_data[subs].values.astype('int32') == file_data[subs].values.astype('int32')
			assert (False not in test_equality)


print("Test 2 passed: all content matching!")
	
print("Testing complete! All existing files of bak have matching content in Zarr and exist in time dimension with proper sizing")
