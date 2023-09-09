#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Reads the subset data from era5 and band-pass filter the meridional wind at selected level


# ### NCSU Tropical and Large Scale Dynamics
# 
# - A. Aiyyer
# - A. Thornton
# 

# In[1]:


import numpy as np
import xarray as xr
import pandas as pd
from datetime import date, timedelta
from numpy import absolute, exp, log

# Any import of metpy will activate the accessors
from metpy.units import units

#from metpy.calc import dewpoint_from_relative_humidity
from metpy.calc import first_derivative, geospatial_gradient, advection

import metpy.constants as constants

import matplotlib.pyplot as plt


# In[16]:


# regridded era5 4x daily data for variables
path_full = '/glade/scratch/athornton/era5_processed_data/3d/' 

# where to store the daily files. 
path_out = '/glade/scratch/athornton/era5_processed_data/filtered_era5_full_fields/600/' 

## Select Variable

############################################################# 

#ERA5 FULL VARIABLES

# Temp
#varNam = 't'
#variab = 'T'
#level = 600
#description = ' T at ' + str(level) + 'hpa'

#U-wind
#varNam = 'u'
#variab = 'U'
#level = 600
#description = ' U-wind ' + str(level) + 'hpa'

# V-wind
#varNam = 'v'
#variab = 'V'
#level = 600
#description = ' V-wind at ' + str(level) + 'hpa'

# omega
#varNam = 'w'
#variab = 'W'
#level = 600
#description = ' Omega at ' + str(level) + 'hpa'

###########################################################

## STATIC ENERGY VARIABLES (FULL)

# MSE
#varNam = 'mse'
#variab = 'MSE'
#level = 600
#description = ' MSE at ' + str(level) + 'hpa'

# DSE
varNam = 'dse'
variab = 'DSE'
level = 600
description = ' DSE at ' + str(level) + 'hpa'

###########################################################

## BUDGET TERMS (FLUXES & RADIATION)

# Surface sensible heat flux
#varNam = 'sshf'
#variab = 'SSHF'
#description = ' Surface sensible heat flux'

# Surface latent heat flux
#varNam = 'slhf'
#variab = 'SLHF'    
#description = ' Surface latent heat flux'

# Q_total
#varNam = 'q'
#variab = 'Q'
#description = ' Net radiative heating'

###########################################################

## BUDGET TERMS (ADVECTION)

# ubar dh'/dx
#varNam = 'ubar_hp_dx_vint'
#variab = 'ubar_hp_dx_vint'
#description = ' background u-advection of anomalous h'

# vbar dh'/dy
#varNam = 'vbar_hp_dy_vint'
#variab = 'vbar_hp_dy_vint'    
#description = ' background v-advection of anomalous h'

# wbar dh'/dp
#varNam = 'wbar_hp_dp_vint'
#variab = 'wbar_hp_dp_vint'    
#description = ' background w-advection of anomalous h'

# u' d/dx (hbar)
#varNam = 'up_hbar_dx_vint'
#variab = 'up_hbar_dx_vint'
#description = ' anomalous u-advection of background h'

# v' d/dy (hbar)
#varNam = 'vp_hbar_dy_vint'
#variab = 'vp_hbar_dy_vint'
#description = ' anomalous v-advection of background h'

# w' d/dp (hbar)
#varNam = 'wp_hbar_dp_vint'
#variab = 'wp_hbar_dp_vint'
#description = ' anomalous w-advection of background h'

# u' dh'/dx
#varNam = 'up_hp_dx_vint'
#variab = 'up_hp_dx_vint'    
#description = ' anomalous u-advection of anomalous h'

# v' dh'/dy
#varNam = 'vp_hp_dy_vint'
#variab = 'vp_hp_dy_vint'
#description = ' anomalous v-advection of anomalous h'

# w' dh'/dp
#varNam = 'wp_hp_dp_vint'
#variab = 'wp_hp_dp_vint'
#description = ' anomalous w-advection of anomalous h'


# define Lancozs filter function
def low_pass_weights(window, cutoff):
    """Calculate weights for a low pass Lanczos filter.

    Args:

    window: int
        The length of the filter window.

    cutoff: float
        The cutoff frequency in inverse time steps.

    """
    order = ((window - 1) // 2) + 1
    nwts = 2 * order + 1
    w = np.zeros([nwts])
    n = nwts // 2
    w[n] = 2 * cutoff
    k = np.arange(1.0, n)
    sigma = np.sin(np.pi * k / n) * n / (np.pi * k)
    firstfactor = np.sin(2.0 * np.pi * cutoff * k) / (np.pi * k)
    w[n - 1 : 0 : -1] = firstfactor * sigma
    w[n + 1 : -1] = firstfactor * sigma
    return w[1:-1]



High_cut_off = 2    # in days
Low_cut_off  = 10   # in days    
data_dt      = 6    # in hours. Here 6 hour interval for era5   
                    # that is, we have data every 6 hours = 4 per day
N_per_day    = int(24/data_dt)

# Window length for filters.
window_in_days = 90 # in days 
window = window_in_days*N_per_day+1
                  
# Construct 2-day low pass filter  
period_cut_off = High_cut_off*24    # same units as data_dt (here hours)
                                    # so 48  = 2 day cut-off in hours
print ('High Period cut-off = ' , period_cut_off)
    
freq = data_dt/period_cut_off
wgts_2d = xr.DataArray(low_pass_weights(window, freq ), dims = ['window'])
   
# Construct 8-day low pass filter  
period_cut_off = Low_cut_off*24    # in hours    
print ('Low Period cut-off = ' , period_cut_off)
    
freq = data_dt/period_cut_off
wgts_10d = xr.DataArray(low_pass_weights(window, freq ), dims = ['window'])

# now get the bandpass filter 2 to 8 days    
wgts_bp = wgts_2d - wgts_10d


# Create list of years to filter
year_start  = 1998
year_end  = 2023
dates_list = np.arange(year_start,year_end, 1)



for a_date in dates_list:
    # filter data temporal extent
    # we want to focus on Jul-Aug_Sep months. To account for some future lag-analysis
    # let us get additional days in the beginning and end
    year_start  = a_date
    month_start = 6
    day_start   = 20

    year_end  = a_date
    month_end   = 10
    day_end     = 10

    year = year_start
    date_filter_start = date(year,month_start,day_start)    
    date_filter_end   = date(year,month_end,day_end)    

    print("User specified range for filtered data = ", date_filter_start, date_filter_end)

    # we need additional data on either side for the filtering

    date_start = date_filter_start - timedelta(days = window_in_days/2)
    date_end = date_filter_end + timedelta(days = window_in_days/2)

    print("User specified buffer days = ", window_in_days)
    print("Need raw data for the time range ", date_start, date_end)
    # list of dates for era5 daily files 
    date_series = [pd.date_range(date_start,date_end, freq ='D') for i in range(year,year+1)]
    # date_series is a list of lists. Lets unpack it now
    dates_list = [element for sublist in date_series for element in sublist]
    data_files = [path_full + varNam + '_' + d.strftime("%Y%m%d") + '.nc' for d in dates_list] 
    ds = xr.open_mfdataset(data_files)
    # read the variable at a particular level
    Var=ds[variab].sel(level=level)
    # now filter it
    Varbpf = Var.rolling(time = len(wgts_bp), center = True).construct('window').dot(wgts_bp).compute()
    Varbpf=Varbpf.assign_attrs(units=Var.units, 
                         #short_name= Var.short_name, 
                         description="2-10 Day filtered ERA5" + description)
    Varbpf=Varbpf.dropna(dim='time')
    Varbpf.name = variab+'_prime'
    file_out = path_out + varNam + '_' + str(level) + '_filtered_' + str(a_date) + '.nc'
    Varbpf.to_netcdf(path=file_out,  format='NETCDF4', mode='w')



