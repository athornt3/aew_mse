#!/usr/bin/env python
# coding: utf-8

# Code to extract era5 data for a specified sub-region and regrid to a coarser grid and also coarser time step
# 
# 
# - NCSU Large Scale and Tropical Dynamics
# 
# Versions
# - A. Aiyyer, Jul 23, 2023
# - A. Thornton, Sep 9, 2023

# In[1]:



import numpy as np
import xarray as xr
import pandas as pd
from datetime import date
from numpy import absolute, exp, log

# Any import of metpy will activate the accessors
from metpy.units import units
import os
import glob

# for regridding
import xesmf as xe

# daily era5
era5_sfc_dir = '/glade/collections/rda/data/ds633.0/e5.oper.an.sfc/'
era5_pl_dir  = '/glade/collections/rda/data/ds633.0/e5.oper.an.pl/'


# output path to save regridded data
path_out = '/glade/scratch/athornton/era5_processed_data/3d/'

#----------------------------------------------------

# SELECT A VARIABLE FROM ERA5

#varId  = '129'
#varNam = 'z'
#variab = 'Z'  # the variable name in the data file

#varId  = '130'
#varNam = 't'
#variab = 'T'

#varId  = '060'
#varNam = 'pv'
#variab = 'PV'

#varId  = '138'
#varNam = 'vo'
#variab = 'VO'

varId  = '133'
varNam = 'q'
variab = 'Q'

#varId  = '135'
#varNam = 'w'
#variab = 'W'

#varId  = '131'
#varNam = 'u'
#variab = 'U'

#varId  = '132'
#varNam = 'v'
#variab = 'V'

#----------------------------------------------------

# SELECT SUBSET OF DATA

# lat/lons
latS = -15.
latN =  35.
lonW = -160.
lonE =  50.

# grid spacing
dlat = 0.5
dlon = 0.5

# levels
level_b = 1000
level_t = 100

# range of dates
year_start  = 1998
month_start = 1
day_start   = 1

year_end  = 2022
month_end = 12
day_end   = 31

date_series = [pd.date_range(date(i,month_start,day_start),date(i,month_end,day_end), freq ='D') for i in range(year_start,year_end+1)]
# date_series is a list of lists. Lets unpack it now
dates_list = [element for sublist in date_series for element in sublist]

print(dates_list[0].strftime("%Y%m%d"))
print(dates_list[-1].strftime("%Y%m%d"))

#----------------------------------------------------


firstPass = True

for a_date in dates_list:
    b_date = a_date + pd.DateOffset(hours=23)
    times = [a_date + pd.DateOffset(hour=h) for h in np.arange(0,24,6)]
    
    # input file name
    fname =   a_date.strftime('%Y%m') + '/e5.oper.an.pl.128_'+varId+'_'+varNam+'.ll025sc.' + a_date.strftime('%Y%m%d')  + '00_' +  a_date.strftime('%Y%m%d') + '23.nc'
    
    infile = era5_pl_dir + fname    
    ds  = xr.open_dataset(infile)  
    
    # prepare to roll the longitude from 0 to 360 --> -180 to 180
    ds = ds.assign_coords(longitude=(((ds.longitude + 180) % 360) - 180))
    
    dat = ds[variab].sel(time=times, latitude=slice(latN,latS), level = slice(level_t, level_b))
    dat = dat.roll(longitude=int(len(dat['longitude']) / 2), roll_coords=True)
    
    
    # create regriddger only once and reuse it afterward
    if (firstPass):
        ds_out = xr.Dataset( 
            {
                "latitude": (["latitude"], np.arange(latN,latS, -dlat),  {"units": "degrees_north"}),
                "longitude": (["longitude"], np.arange(lonW, lonE, dlon), {"units": "degrees_east"}),

            }
        )
        ds_in = xr.Dataset(
            {
                "latitude": (["latitude"], dat.latitude.values,  {"units": "degrees_north"}),
                "longitude": (["longitude"], dat.longitude.values, {"units": "degrees_east"}),
            }
        )
        regridder = xe.Regridder(ds_in, ds_out, "conservative")
        firstPass = False
        
        
    # regrid and write to file (4x daily)
    dat_out = regridder(dat, keep_attrs=True)    
    file_out = path_out + varNam+ '_' + a_date.strftime("%Y%m%d") + '.nc'
    dat_out.to_netcdf(path=file_out)