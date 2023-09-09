#!/usr/bin/env python
# coding: utf-8

# ## Long Term Average Fields
# This python script downloads the long term JAS average for era5 variables:
# 
# A. Thornton
# 
# Jul 2023

# Imports
import numpy as np
import xarray as xr
import pandas as pd
from datetime import date
from numpy import absolute, exp, log

from metpy.units import units
import os
import glob


#---------------------------------------------------------------------------
## SELECT VARIABLE

#varNam = 'mse'
#variab = 'MSE'

# Potential Vorticity
#varId  = '060'
#varNam = 'pv'
#variab = 'PV'

# Vorticity
#varId  = '138'
#varNam = 'vo'
#variab = 'VO'

# Specific humidity
#varId  = '133'
#varNam = 'q'
#variab = 'Q'

# Geopotential
#varId  = '129'
#varNam = 'z'
#variab = 'Z'

# Temperature
#varId  = '130'
#varNam = 't'
#variab = 'T'

# U-wind
#varId  = '131'
#varNam = 'u'
#variab = 'U'

# V-wind
#varId  = '132'
#varNam = 'v'
#variab = 'V'

# W-wind
#varId  = '135'
#varNam = 'w'
#variab = 'W'

# Surface sensible heat flux
varId  = '146'
varNam = 'sshf'
variab = 'SSHF'

# Surface latent heat flux
#varId  = '147'
#varNam = 'slhf'
#variab = 'SLHF'    

# Surface net short-wave (solar) radiation 
#varId  = '176'
#varNam = 'ssr'
#variab = 'SSR'
    
# Surface net long-wave (thermal) radiation
#varId  = '177'
#varNam = 'str'
#variab = 'STR'

# Top net short-wave (solar) radiation
#varId  = '178'
#varNam = 'tsr'
#variab = 'TSR'

# Top net long-wave (thermal) radiation
#varId  = '179'
#varNam = 'ttr'
#variab = 'TTR'

#---------------------------------------------------------------------------

# SELECT SUBSET OF DATES
year_start  = 1998
month_start = 7
day_start   = 1

year_end  = 2022
month_end = 9
day_end   = 30

date_series = [pd.date_range(date(i,month_start,day_start),date(i,month_end,day_end), freq ='D') for i in range(year_start,year_end+1)]
# date_series is a list of lists. Lets unpack it now
dates_list = [element for sublist in date_series for element in sublist]


path_to_files = '/glade/scratch/athornton/era5_processed_data/2d/2d_updated/'
os.chdir(path_to_files)
infiles = [varNam+ '_unfiltered_' + date.strftime("%Y%m%d") + '.nc' for date in dates_list]
ds = xr.open_mfdataset(infiles)

# CALCULATE MONTHLY MEAN
# keep_attrs retains units
ds = ds.mean(dim='time', keep_attrs=True)
ds = ds[variab].compute()

# PATH TO SAVE FILE
path_out = '/glade/scratch/athornton/era5_processed_data/jas_means/'

# SAVE FILE
file_out = path_out + varNam+ '_' + 'jas_'+ str(year_start) + '_' + str(year_end) + '.nc'
ds.to_netcdf(path=file_out)




