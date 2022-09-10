#!/usr/bin/env python

# Author: C. Max Stevens (@maximusjstevens)
# christopher.d.stevens-1@nasa.gov
#
# This script is effectively a MERRA-2 subsetter. It opens all files of a
# specific MERRA-2 file type (e.g., time-averaged Surface Flux Diagnostics; see
# https://gmao.gsfc.nasa.gov/pubs/docs/Bosilovich785.pdf), does a geographic subsetting
# and pulls just the variables you specify. It then resamples to daily resolution, and
# saves a .nc file for each year, containing all of the variables (i.e. it concatentates 
# from the numerous file types into a single file).
#
# Then, it takes each of those yearly files and combines them in a zarr.
# 
# to run this on discover, first load python module:
# module load python/GEOSpyD/Min4.11.0_py3.9
#
# MERRA-2 on discover are here:
# /discover/nobackup/projects/gmao/merra2/data/products
#
# Takes advantage of the "MERRA2_all" directory, as described in the
# MERRA-2 readme on discover.

import math, cmath
import numpy as np
from numba import jit
import netCDF4 as nc
import h5py as h5
import datetime
import s3fs
import xarray as xr
import pandas as pd
import calendar
import os
import sys
import glob
import time
import pickle
import re
import dateutil.parser as dparser

class FQS:
    '''
    Fast Quartic Solver: analytically solves quartic equations (needed to calculate melt)
    Takes methods from fqs package (@author: NKrvavica)
    full documentation: https://github.com/NKrvavica/fqs/blob/master/fqs.py
    '''
    def __init__(self):
        pass

    @jit(nopython=True)
    def single_quadratic(self, a0, b0, c0):
        ''' 
        Analytical solver for a single quadratic equation
        '''
        a, b = b0 / a0, c0 / a0

        # Some repating variables
        a0 = -0.5*a
        delta = a0*a0 - b
        sqrt_delta = cmath.sqrt(delta)

        # Roots
        r1 = a0 - sqrt_delta
        r2 = a0 + sqrt_delta

        return r1, r2


    @jit(nopython=True)
    def single_cubic(self, a0, b0, c0, d0):
        ''' 
        Analytical closed-form solver for a single cubic equation
        '''
        a, b, c = b0 / a0, c0 / a0, d0 / a0

        # Some repeating constants and variables
        third = 1./3.
        a13 = a*third
        a2 = a13*a13
        sqr3 = math.sqrt(3)

        # Additional intermediate variables
        f = third*b - a2
        g = a13 * (2*a2 - b) + c
        h = 0.25*g*g + f*f*f

        def cubic_root(x):
            ''' Compute cubic root of a number while maintaining its sign'''
            if x.real >= 0:
                return x**third
            else:
                return -(-x)**third

        if f == g == h == 0:
            r1 = -cubic_root(c)
            return r1, r1, r1

        elif h <= 0:
            j = math.sqrt(-f)
            k = math.acos(-0.5*g / (j*j*j))
            m = math.cos(third*k)
            n = sqr3 * math.sin(third*k)
            r1 = 2*j*m - a13
            r2 = -j * (m + n) - a13
            r3 = -j * (m - n) - a13
            return r1, r2, r3

        else:
            sqrt_h = cmath.sqrt(h)
            S = cubic_root(-0.5*g + sqrt_h)
            U = cubic_root(-0.5*g - sqrt_h)
            S_plus_U = S + U
            S_minus_U = S - U
            r1 = S_plus_U - a13
            r2 = -0.5*S_plus_U - a13 + S_minus_U*sqr3*0.5j
            r3 = -0.5*S_plus_U - a13 - S_minus_U*sqr3*0.5j
            return r1, r2, r3


    @jit(nopython=True)
    def single_cubic_one(self, a0, b0, c0, d0):
        ''' 
        Analytical closed-form solver for a single cubic equation
        '''
        a, b, c = b0 / a0, c0 / a0, d0 / a0

        # Some repeating constants and variables
        third = 1./3.
        a13 = a*third
        a2 = a13*a13

        # Additional intermediate variables
        f = third*b - a2
        g = a13 * (2*a2 - b) + c
        h = 0.25*g*g + f*f*f

        def cubic_root(x):
            ''' Compute cubic root of a number while maintaining its sign
            '''
            if x.real >= 0:
                return x**third
            else:
                return -(-x)**third

        if f == g == h == 0:
            return -cubic_root(c)

        elif h <= 0:
            j = math.sqrt(-f)
            k = math.acos(-0.5*g / (j*j*j))
            m = math.cos(third*k)
            return 2*j*m - a13

        else:
            sqrt_h = cmath.sqrt(h)
            S = cubic_root(-0.5*g + sqrt_h)
            U = cubic_root(-0.5*g - sqrt_h)
            S_plus_U = S + U
            return S_plus_U - a13


    @jit(nopython=True)
    def single_quartic(self, a0, b0, c0, d0, e0):
        '''
        Analytical closed-form solver for a single quartic equation
        '''
        a, b, c, d = b0/a0, c0/a0, d0/a0, e0/a0

        # Some repeating variables
        a0 = 0.25*a
        a02 = a0*a0

        # Coefficients of subsidiary cubic euqtion
        p = 3*a02 - 0.5*b
        q = a*a02 - b*a0 + 0.5*c
        r = 3*a02*a02 - b*a02 + c*a0 - d

        # One root of the cubic equation
        z0 = self.single_cubic_one(1, p, r, p*r - 0.5*q*q)

        # Additional variables
        s = cmath.sqrt(2*p + 2*z0.real + 0j)
        if s == 0:
            t = z0*z0 + r
        else:
            t = -q / s

        # Compute roots by quadratic equations
        r0, r1 = self.single_quadratic(1, s, z0 + t)
        r2, r3 = self.single_quadratic(1, -s, z0 - t)

        return r0 - a0, r1 - a0, r2 - a0, r3 - a0


    def multi_quadratic(self, a0, b0, c0):
        ''' 
        Analytical solver for multiple quadratic equations
        '''
        a, b = b0 / a0, c0 / a0

        # Some repating variables
        a0 = -0.5*a
        delta = a0*a0 - b
        sqrt_delta = np.sqrt(delta + 0j)

        # Roots
        r1 = a0 - sqrt_delta
        r2 = a0 + sqrt_delta

        return r1, r2


    def multi_cubic(self, a0, b0, c0, d0, all_roots=True):
        '''
        Analytical closed-form solver for multiple cubic equations
        '''
        a, b, c = b0 / a0, c0 / a0, d0 / a0

        # Some repeating constants and variables
        third = 1./3.
        a13 = a*third
        a2 = a13*a13
        sqr3 = math.sqrt(3)

        # Additional intermediate variables
        f = third*b - a2
        g = a13 * (2*a2 - b) + c
        h = 0.25*g*g + f*f*f

        # Masks for different combinations of roots
        m1 = (f == 0) & (g == 0) & (h == 0)     # roots are real and equal
        m2 = (~m1) & (h <= 0)                   # roots are real and distinct
        m3 = (~m1) & (~m2)                      # one real root and two complex

        def cubic_root(x):
            ''' Compute cubic root of a number while maintaining its sign
            '''
            root = np.zeros_like(x)
            positive = (x >= 0)
            negative = ~positive
            root[positive] = x[positive]**third
            root[negative] = -(-x[negative])**third
            return root

        def roots_all_real_equal(c):
            ''' Compute cubic roots if all roots are real and equal
            '''
            r1 = -cubic_root(c)
            if all_roots:
                return r1, r1, r1
            else:
                return r1

        def roots_all_real_distinct(a13, f, g, h):
            ''' Compute cubic roots if all roots are real and distinct
            '''
            j = np.sqrt(-f)
            k = np.arccos(-0.5*g / (j*j*j))
            m = np.cos(third*k)
            r1 = 2*j*m - a13
            if all_roots:
                n = sqr3 * np.sin(third*k)
                r2 = -j * (m + n) - a13
                r3 = -j * (m - n) - a13
                return r1, r2, r3
            else:
                return r1

        def roots_one_real(a13, g, h):
            ''' Compute cubic roots if one root is real and other two are complex
            '''
            sqrt_h = np.sqrt(h)
            S = cubic_root(-0.5*g + sqrt_h)
            U = cubic_root(-0.5*g - sqrt_h)
            S_plus_U = S + U
            r1 = S_plus_U - a13
            if all_roots:
                S_minus_U = S - U
                r2 = -0.5*S_plus_U - a13 + S_minus_U*sqr3*0.5j
                r3 = -0.5*S_plus_U - a13 - S_minus_U*sqr3*0.5j
                return r1, r2, r3
            else:
                return r1

        # Compute roots
        if all_roots:
            roots = np.zeros((3, len(a))).astype(complex)
            roots[:, m1] = roots_all_real_equal(c[m1])
            roots[:, m2] = roots_all_real_distinct(a13[m2], f[m2], g[m2], h[m2])
            roots[:, m3] = roots_one_real(a13[m3], g[m3], h[m3])
        else:
            roots = np.zeros(len(a))  # .astype(complex)
            roots[m1] = roots_all_real_equal(c[m1])
            roots[m2] = roots_all_real_distinct(a13[m2], f[m2], g[m2], h[m2])
            roots[m3] = roots_one_real(a13[m3], g[m3], h[m3])

        return roots


    def multi_quartic(self, a0, b0, c0, d0, e0):
        ''' 
        Analytical closed-form solver for multiple quartic equations
        '''
        a, b, c, d = b0/a0, c0/a0, d0/a0, e0/a0

        # Some repeating variables
        a0 = 0.25*a
        a02 = a0*a0

        # Coefficients of subsidiary cubic euqtion
        p = 3*a02 - 0.5*b
        q = a*a02 - b*a0 + 0.5*c
        r = 3*a02*a02 - b*a02 + c*a0 - d

        # One root of the cubic equation
        z0 = self.multi_cubic(1, p, r, p*r - 0.5*q*q, all_roots=False)

        # Additional variables
        s = np.sqrt(2*p + 2*z0.real + 0j)
        t = np.zeros_like(s)
        mask = (s == 0)
        t[mask] = z0[mask]*z0[mask] + r[mask]
        t[~mask] = -q[~mask] / s[~mask]

        # Compute roots by quadratic equations
        r0, r1 = self.multi_quadratic(1, s, z0 + t) - a0
        r2, r3 = self.multi_quadratic(1, -s, z0 - t) - a0

        return r0, r1, r2, r3


    def cubic_roots(self, p):
        '''
        A caller function for a fast cubic root solver (3rd order polynomial).
        '''
        # Convert input to array (if input is a list or tuple)
        p = np.asarray(p)

        # If only one set of coefficients is given, add axis
        if p.ndim < 2:
            p = p[np.newaxis, :]

        # Check if four coefficients are given
        if p.shape[1] != 4:
            raise ValueError('Expected 3rd order polynomial with 4 '
                             'coefficients, got {:d}.'.format(p.shape[1]))

        if p.shape[0] < 100:
            roots = [self.single_cubic(*pi) for pi in p]
            return np.array(roots)
        else:
            roots = self.multi_cubic(*p.T)
            return np.array(roots).T


    def quartic_roots(self, p):
        '''
        A caller function for a fast quartic root solver (4th order polynomial).
        '''
        # Convert input to an array (if input is a list or tuple)
        p = np.asarray(p)

        # If only one set of coefficients is given, add axis
        if p.ndim < 2:
            p = p[np.newaxis, :]

        # Check if all five coefficients are given
        if p.shape[1] != 5:
            raise ValueError('Expected 4th order polynomial with 5 '
                             'coefficients, got {:d}.'.format(p.shape[1]))

        if p.shape[0] < 100:
            roots = [self.single_quartic(*pi) for pi in p]
            return np.array(roots)
        else:
            roots = self.multi_quartic(*p.T)
            return np.array(roots).T

class MERRA_concat:
    '''
    Class to open the native MERRA-2 files, make a geographic subset, and pull 
    just the variables you want.

    Then resamples to daily resolution. All resampling is done using mean(), 
    which preserves units. 

    Could easily be altered to get different variables from different MERRA-2 products
    e.g., you could alter to get data from tavg1_2d_lnd_Nx (M2T1NXLND): Land Surface Diagnostics.
    To do so, follow the code format for the rad, int, slv, and flx. 
    '''
    def __init__(self):
        pass

    def read_netcdfs(self, files, lat_min, lat_max, lon_min, lon_max, var_list):
        '''
        read the netcdfs and subset. Appears faster than xr.open_mfdataset.
        '''
        def process_one_path(path):
            with xr.open_dataset(path) as ds:
                dss = ds.sel({'lat': slice(lat_min,lat_max),'lon': slice(lon_min,lon_max)})[var_list]
                dss.load()  # load data to ensure we can use it after closing each original file
                return dss

        paths = sorted(files)
        datasets = [process_one_path(p) for p in paths]
        combined = xr.concat(datasets, dim='time').sortby('time')
        return combined

    def M2combine(self,icesheet,out_path,allyears,calc_melt=True):
        '''
        main function to open all files, optionally calculate melt, resample to daily, 
        and merge all into a single netCDF for each year.
        '''

        if icesheet=='GrIS':
            lat_min=55
            lat_max=90
            lon_min=-80
            lon_max=-10
        elif icesheet=='AIS':
            lat_min=-90
            lat_max=-60
            lon_min=-180
            lon_max=180

        tnow = time.time() #time it

        ### lists of variables to get from each MERRA-2 file type
        rad_list = ['ALBEDO','LWGAB','LWGEM','LWGNT','SWGDN','SWGNT','TS']
        int_list = ['PRECCU','PRECLS','PRECSN']
        flx_list = ['EVAP','HFLUX','EFLUX','GHTSKIN']
        slv_list = ['PS','T2M','U2M','V2M','QV2M']

        for YY in allyears:

            out_filename = f'MERRA2_{icesheet}_Daily_{YY}.nc'

            if os.path.exists(out_path+out_filename):
                print(out_filename, 'found')
                print('moving to next year')
                continue

            print(YY)
            merra_path = f'/discover/nobackup/projects/gmao/merra2/data/products/MERRA2_all/Y{YY}/'
            
            ### lists of filenames for each M2 type
            fns_rad = sorted(glob.glob(merra_path+f'**/MERRA*tavg1_2d_rad_Nx.{YY}*.nc4'))
            fns_int = sorted(glob.glob(merra_path+f'**/MERRA*tavg1_2d_int_Nx.{YY}*.nc4'))
            fns_flx = sorted(glob.glob(merra_path+f'**/MERRA*tavg1_2d_flx_Nx.{YY}*.nc4'))
            fns_slv = sorted(glob.glob(merra_path+f'**/MERRA*tavg1_2d_slv_Nx.{YY}*.nc4'))

            ### concatenated datasets (i.e. all data for one year); these are hourly, geographically
            ### subsetted, and only contain the variables in the variable lists above
            ds_rad = self.read_netcdfs(fns_rad,lat_min, lat_max, lon_min, lon_max, rad_list)
            ds_int = self.read_netcdfs(fns_int,lat_min, lat_max, lon_min, lon_max, int_list)
            ds_flx = self.read_netcdfs(fns_flx,lat_min, lat_max, lon_min, lon_max, flx_list)
            ds_slv = self.read_netcdfs(fns_slv,lat_min, lat_max, lon_min, lon_max, slv_list)

            ### resample to daily resolution, keep attributes from original files
            for kk,rvar in enumerate(rad_list):
                if kk==0:
                    ds_rad_r = ds_rad[rvar].resample(time = '1D').mean().to_dataset()
                else:
                    ds_rad_r[rvar] = ds_rad[rvar].resample(time = '1D').mean()
                    ds_rad_r[rvar].attrs = ds_rad[rvar].attrs

            for kk,ivar in enumerate(int_list):
                if kk==0:
                    ds_int_r = ds_int[ivar].resample(time = '1D').mean().to_dataset()
                else:
                    ds_int_r[ivar] = ds_int[ivar].resample(time = '1D').mean()
                    ds_int_r[ivar].attrs = ds_int[ivar].attrs

            for kk,fvar in enumerate(flx_list):
                if kk==0:
                    ds_flx_r = ds_flx[fvar].resample(time = '1D').mean().to_dataset()
                else:
                    ds_flx_r[fvar] = ds_flx[fvar].resample(time = '1D').mean()
                    ds_flx_r[fvar].attrs = ds_flx[fvar].attrs

            for kk,svar in enumerate(slv_list):
                if kk==0:
                    ds_slv_r = ds_slv[svar].resample(time = '1D').mean().to_dataset()
                else:
                    ds_slv_r[svar] = ds_slv[svar].resample(time = '1D').mean()
                    ds_slv_r[svar].attrs = ds_slv[svar].attrs
            ###

            ### Below code calculates surface melt and Tsurf using radiation fluxes ###
            ### assumes a surface layer of 0.08 m with density=400 kg/m3
            if calc_melt:
                SBC = 5.67e-8
                CP_I = 2097.0 
                m = 400*0.08
                dt = 3600
                LF_I = 333500.0 #[J kg^-1]
                flux_df1 = (ds_rad['SWGNT'] + ds_rad['LWGAB'] - ds_flx['HFLUX'] - ds_flx['EFLUX'] + ds_flx['GHTSKIN'])
                flux_df1_r = flux_df1.values.reshape(flux_df1.shape[0],-1)
                oshape = np.shape(flux_df1.values)
                Tcalc = np.zeros_like(flux_df1_r)
                meltmass = np.zeros_like(flux_df1_r)

                dts = ds_rad.TS.values
                dts_r=dts.reshape(dts.shape[0],-1)
                dsha = dts_r.shape[-1]

                tindex = ds_rad.time.data
                fqs = FQS()
                for kk,mdate in enumerate(tindex):
                    pmat = np.zeros((dsha,5))
                    if kk==0:
                        T_0 = dts_r[kk,:]
                    else:
                        T_0 = dts_r[kk-1,:]

                    a = SBC*dt/(CP_I*m) * np.ones(dsha)
                    b = 0
                    c = 0
                    d = 1 * np.ones(dsha)
                    e = -1 * (flux_df1_r[kk,:]*dt/(CP_I*m)+T_0)

                    pmat[:,0] = a
                    pmat[:,3] = d
                    pmat[:,4] = e
                    pmat[np.isnan(pmat)] = 0

                    r = fqs.quartic_roots(pmat)
                    Tnew = (r[((np.isreal(r)) & (r>0))].real)
                    Tnew[np.isnan(e)] = np.nan
                    
                    imelt = np.where(Tnew>=273.15)[0]
                    Tnew[imelt] = 273.15
                    meltmass[kk,imelt] = (flux_df1_r[kk,imelt] - SBC*273.15**4) / LF_I #*dt #multiply by dt to put in units per day
                    ### meltmass has units [kg/m2/s]            
                    Tcalc[kk,:] = Tnew

                Tcalc_out = Tcalc.reshape(oshape)
                meltmass_out = meltmass.reshape(oshape)

                melt_array =  xr.DataArray(meltmass_out,dims=ds_rad.dims)
                tcalc_array = xr.DataArray(Tcalc_out,dims=ds_rad.dims)
                tcalc_array.attrs = ds_rad['TS'].attrs
                melt_array.attrs =  {'long_name': 'meltwater_production',
                                     'units': 'kg m-2 s-1',
                                     'fmissing_value':  1000000000000000.0,
                                     'standard_name':   'meltwater',
                                     'vmax':  1000000000000000.0,
                                     'vmin':  -1000000000000000.0,             
                                     'fullnamepath':    '/SMELT'
                                     }

                ds_rad['SMELT'] = melt_array
                ds_rad['TScalc'] = tcalc_array

                ds_rad_r['TScalc']     = ds_rad.TScalc.resample(time = '1D').mean()
                ds_rad_r['SMELT']     = ds_rad.SMELT.resample(time = '1D').mean()
            ### end melt calculation

            for vv in ds_rad_r.variables:
                ds_rad_r[vv].attrs = ds_rad[vv].attrs

            ### merge all datasets into one
            dsALL = ((ds_int_r.merge(ds_rad_r)).merge(ds_flx_r)).merge(ds_slv_r)

            dsALL.attrs = {
                'SouthernmostLatitude': str(dsALL.lat.min().values),
                'NorthernmostLatitude': str(dsALL.lat.max().values),
                'WesternmostLongitude': str(dsALL.lon.min().values),
                'EasternmostLongitude': str(dsALL.lon.max().values),
                'RangeBeginningDate': str(dsALL.time[0].values),
                'RangeEndingDate': str(dsALL.time[-1].values),
                'comment1': 'This is a custom resample of MERRA2 data used to force the Community Firn Model',
                'comment2': 'Data come from: tavg1_2d_rad_Nx, tavg1_2d_int_Nx, tavg1_2d_slv_Nx, tavg1_2d_flx_Nx',
                'comment3': 'Files are produced using MERRA_concat.py',
                'email': 'christopher.d.stevens-1@nasa.gov',
                'filename': out_filename
            }

            dsALL.to_netcdf(out_path + out_filename) #save as netcdf
            telap = (time.time()-tnow)/60
            print('iteration time:', telap)
            tnow=time.time()
        ### end loop through years

class make_zarr:
    '''
    class to take the daily resolution files (one for each year), 
    combine them and save as zarr
    optionally can save a big netCDF will all years in one file.
    '''
    def __init__(self):
        pass

    def read_daily_netcdfs(self,files):
        def process_one_path(path):
            print(path)
            with xr.open_dataset(path) as ds:
                ds.load()
                return ds
        paths = sorted(files)
        datasets = [process_one_path(p) for p in paths]
        combined = xr.concat(datasets, dim='time').sortby('time')
        return combined

    def make_zarr(self,icesheet,out_path,saveZarr=True,saveBigNC=False):
        allYearlyFiles = glob.glob(out_path + f'MERRA2_{icesheet}_Daily_*.nc')
        dsALL = self.read_daily_netcdfs(allYearlyFiles)        
        if saveBigNC: #optionally save a single netCDF with all year (big file)
            dsALL.to_netcdf(out_path + '/MERRA_Greenland_Daily_ALL.nc')   
        if saveZarr: # Put the .nc file into a zarr
            dsALL.to_zarr(store=os.getenv('NOBACKUP')+'/M2_GrIS_daily_melt.zarr',mode='w',consolidated=True) # Be aware of where you are saving this; does not go to merra_path

if __name__ == '__main__':
    icesheet='GrIS'
    out_path = os.getenv('NOBACKUP')+f'/MERRA2/Daily/{icesheet}/'
    if os.path.exists(out_path):
        pass
    else:
        os.makedirs(out_path)
    allyears = np.arange(1980,2021) # years to resample
    Mc = MERRA_concat()
    Mc.M2combine(icesheet,out_path,allyears)
    print('Done with making yearly files with daily resolution.')
    print('Now making zarr.')
    MZ = make_zarr()
    MZ.make_zarr(icesheet,out_path)


