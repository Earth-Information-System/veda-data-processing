#!/usr/bin/env python
'''
© 2021. California Institute of Technology. Government sponsorship acknowledge.
Script to open the MERRA-2 hourly (and 3 hourly for runoff) files
Resample to daily
Then save it in a zarr
'''

import numpy as np
import pandas as pd
import h5py as h5
import os
import sys
import pickle
import netCDF4 as nc
import glob
import xarray as xr
import re
import dateutil.parser as dparser

year1=1980
year2=2020
y = np.r_[year1:(year2+1)]
m = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
#m = ['01']

merra_path = '/discover/nobackup/projects/gmao/merra2/data/products/MERRA2_all/'
merra_out_path = '/discover/nobackup/projects/eis_sealevel/MERRA-2/'
# merra_path = merra_path = '/Volumes/Samsung_T1/MERRA/' #Max local
bbox_lon = [-80, -10] # bounding box
bbox_lat = [55, 90]

for yi in y:
	for mi in m:
		print('Year ' + str(yi) + ', Month ' +  mi)

		mdirSMB = merra_path + 'Y' + (str(yi)) + '/M' + mi + '/' # this is the directory that contains the tavg1_2d_flx_Nx files (SMB fluxes)
		# mdirSMB = merra_path + 'Greenland/Hourly/SMB/' # Max local
		fnsSMB = glob.glob(mdirSMB + '*tavg1_2d_flx_Nx*.nc*') # all of the netCDF files
		# fnsSMB = [ff[ii] for ii in ff]
		
		mdirTS = merra_path + 'Y' + (str(yi)) + '/M' + mi + '/' # this is the directory that contains the tavg1_2d_slv_Nx files (Temperatures)
		# mdirTS = merra_path + 'Greenland/Hourly/TS/' # Max local
		fnsTS = glob.glob(mdirTS + '*tavg1_2d_slv_Nx*.nc*')
		
		#I don't have the runoff files, so not sure what the file name type is.
		mdirRU = merra_path + 'Y' + (str(yi)) + '/M' + mi + '/' # this is the directory that contains the runoff files (not sure if I have the code correct; I do not have these on my computer)
		fnsRU = glob.glob(mdirRU + '*tavg3_2d_glc_Nx*.nc*')
		
		datasetsSMB = [xr.open_dataset(fname, chunks={}) for fname in fnsSMB]
		datasetsTS = [xr.open_dataset(fname, chunks={}) for fname in fnsTS]
		datasetsRU = [xr.open_dataset(fname, chunks={}) for fname in fnsRU]
		
		dsSMB = (xr.concat(datasetsSMB,dim='time')).sortby('time').sel(lon=slice(bbox_lon[0],bbox_lon[1]),lat=slice(bbox_lat[0],bbox_lat[1]))
		dsTS = (xr.concat(datasetsTS,dim='time')).sortby('time').sel(lon=slice(bbox_lon[0],bbox_lon[1]),lat=slice(bbox_lat[0],bbox_lat[1]))
		dsRU = (xr.concat(datasetsRU,dim='time')).sortby('time').sel(lon=slice(bbox_lon[0],bbox_lon[1]),lat=slice(bbox_lat[0],bbox_lat[1]))
		
		dsSMBr = (dsSMB.PRECTOT.resample(time = '1D').sum()*3600).to_dataset() 
		# I multiply by 3600 because MERRA-2 fluxes are in kg m^-2 s^-1, which I think is annoying during resampling - my preference is to have units of kg m^-2 per time resultion (e.g. per hour or per day). That way the data represents how much precip actually fell during that time interval, rather than how much fell per second over that time interval. But, delete the 3600 if you want (just tell me that you did!)
		dsSMBr['EVAP'] = dsSMB.EVAP.resample(time = '1D').sum()*3600
		
		# Not sure if this next line will work! (But I think it should)
		dsSMBr['RUNOFF'] = dsRU.RUNOFF.resample(time = '1D').sum()*3*3600 
		# Again, the *3*3600 give the runoff during that time interval (MERRA2 are 3 hourly, so values in native 3 hourly M2 data need to be multiplied by 3*3600 to get amount of runoff during those 3 hours, then summing for the resample to daily will give the total mass of runoff during that day.) 
		
		#dsTSr = dsTS.resample(time='1D').mean()
		dsTSr = (dsTS.T2M.resample(time = '1D').mean()).to_dataset()
		dsTSr['TS']=dsTS['TS'].resample(time='1D').mean()

		dsALL = dsTSr.merge(dsSMBr)
		
		dsALL.to_zarr(store = merra_out_path + 'MERRA-2_Greenland_' + (str(yi)) + '_' + mi  + '.zarr', mode = 'w')
