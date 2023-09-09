#!/usr/bin/env python
# coding: utf-8

# ## Code to aggregate ERA5 fluxes over 6 hour intervals
# 
# 
# Notes
# 
# - Fuxes are (1) Forecasts made at 6Z and 18Z; and (2) are accumulated over the past hour.
# 
#     - For any day (D), 6 Z forecast_initial_time, Take forecast hours 1-6: Sum the fluxes and assign to 12 Z  for that day D
#     - For any day (D), 6 Z forecast_initial_time, Take forecast hours 7-12: Sum the fluxes and assign to 18 Z for that day D
#     - For any day (D), 18 Z forecast_initial_time, Take forecast hours 1-6: Sum the fluxes and assign to 0 Z  for D+1
#     - For any day (D), 18 Z forecast_initial_time, Take forecast hours 7-12: Sum the fluxes and assign to 6 Z for D+1
#   
# 
# ## Note
# 
# - Fluxes from NCAR RDA are actually coded as 'Watts per m-sq s' which is same as 'Joules per m-sq'
# 
# 
# ### NCSU Large Scale and Tropical Dynamics
# - A. Aiyyer (Jul 2023)
# - Sep 6: Updated to divide the 6 hourly accumulated flux by (6*3600) and implemented metpy


import numpy as np
import xarray as xr
import pandas as pd
from datetime import date, datetime
from numpy import absolute, exp, log

# Any import of metpy will activate the accessors
from metpy.units import units
import os
import glob

# for regridding
import xesmf as xe

import metpy

#------------------------------------------------------------------------

# SELECT SUBSET OF DATA

# grid edges
latS = -15.
latN =  35.
lonW = -160.
lonE =  50.

# grid spacing
dlat = 0.5
dlon = 0.5

# dates
year_start = 1998
year_end = 2022

#------------------------------------------------------------------------

# SELECT ERA5 FLUXES
# The ECMWF convention for vertical fluxes is positive downwards

# Surface sensible heat flux
#https://codes.ecmwf.int/grib/param-db/?id=146
#varId  = '146'
#varNam = 'sshf'
#variab = 'SSHF'

# Surface latent heat flux
#https://codes.ecmwf.int/grib/param-db/?id=146
#varId  = '147'
#varNam = 'slhf'
#variab = 'SLHF'    

# Surface net short-wave (solar) radiation 
# https://codes.ecmwf.int/grib/param-db/?id=176
#varId  = '176'
#varNam = 'ssr'
#variab = 'SSR'
   
# Surface net long-wave (thermal) radiation
# https://codes.ecmwf.int/grib/param-db/?id=177
#varId  = '177'
#varNam = 'str'
#variab = 'STR'

# Top net short-wave (solar) radiation
# https://codes.ecmwf.int/grib/param-db/?id=178
varId  = '178'
varNam = 'tsr'
variab = 'TSR'

# Top net long-wave (thermal) radiation
# https://codes.ecmwf.int/grib/param-db/?id=179
#varId  = '179'
#varNam = 'ttr'
#variab = 'TTR'

#------------------------------------------------------------------------

# daily era5
er5_sfc_accu = '/glade/collections/rda/data/ds633.0/e5.oper.fc.sfc.accumu/'

# output path for regridded data
path_out = '/glade/scratch/athornton/era5_processed_data/2d/2d_updated/'


# function for getting files
def fileNames (year):
    start = date(year-1,12,1)
    end   = date(year,12,1)
    dates_list = pd.date_range(start, end, freq='MS')
    
    fpath = er5_sfc_accu + dates_list.strftime('%Y%m') + '/e5.oper.fc.sfc.accumu.128_'+varId+'_'+varNam  
    fils = []
    for f in fpath:
        fils.append(glob.glob(f + '*nc'))

    infiles = [x for l in fils for x in l]
    
    return infiles
    
# function for selecting subset of data
def preprocess(ds):
    ds = ds.assign_coords(longitude=(((ds.longitude + 180) % 360) - 180))
    ds = ds.roll(longitude=int(len(ds['longitude']) / 2), roll_coords=True)
    ds = ds.sel(latitude=slice(latN,latS), longitude=slice(lonW,lonE))    
    return ds


firstPass = True
 
for year in range (year_start, year_end + 1):
    infiles = fileNames (year)
    #print(infiles)

    ds = xr.open_mfdataset(infiles,  preprocess=preprocess)

    # create regriddger only once and reuse it afterward
    if (firstPass):
        ds_out = xr.Dataset( 
            {
                "latitude": (["latitude"], np.arange(latN,latS, -dlat),  {"units": "degrees_north"}),
                "longitude": (["longitude"], np.arange(lonW, lonE, dlon), {"units": "degrees_east"}),

            }
        )
        ds_out.attrs = ds.attrs
        regridder = xe.Regridder(ds, ds_out, "conservative")
        firstPass = False
         
    dat_out = regridder(ds[variab], keep_attrs=True)  
    
    
    # now sum over 6 hours for 0,6,12 and 18 Z accumulations
    datA=dat_out.sel(forecast_hour=slice(1,6)).sum(dim='forecast_hour', keep_attrs=True)
    datB=dat_out.sel(forecast_hour=slice(7,12)).sum(dim='forecast_hour', keep_attrs=True)

    
    # divide the 6-hourly accumulation by 6*3600 seconds
    # note this division removes the attributes from datA and datB. we will add them later
    datA = datA.metpy.quantify()/(6*3600*units('s'))  
    datB = datB.metpy.quantify()/(6*3600*units('s')) 
   
    
    # adjust the time stamps to match the accumulations and rename forecast time to time
    datA['forecast_initial_time'] = datA.forecast_initial_time + pd.Timedelta(6, "h")
    datA = datA.rename({'forecast_initial_time': 'time'})
    datB['forecast_initial_time'] = datB.forecast_initial_time + pd.Timedelta(12, "h")
    datB = datB.rename({'forecast_initial_time': 'time'})

    dat_combined = datA.combine_first(datB)
    dat_combined = dat_combined.metpy.dequantify()

    
    # add some attributes
    dat_combined.attrs['info'] = '6 hour accumu from 6Z and 18Z forecasts. By NCSU Tropical Dynamics'
    dat_combined.attrs['long_name'] = dat_out.attrs['long_name']
    dat_combined.attrs['short_name'] = dat_out.attrs['short_name']
         
    # now loop over each day and write data to netcdf files (4x daily)
    for date1 in pd.date_range(str(year) + '-01-01-00', str(year) + '-12-31-00' , freq='D'):
        date2 = date1 + pd.Timedelta(18, "h")
        #print (date1, date2)          
        file_out = path_out + varNam+ '_' + date1.strftime("%Y%m%d") + '.nc'
        dat_combined.sel(time = slice(date1,date2)).to_netcdf(path=file_out, format='NETCDF4', mode='w')
        