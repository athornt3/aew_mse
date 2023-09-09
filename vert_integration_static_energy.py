#!/usr/bin/env python
# coding: utf-8

# Code to calculate vertical integration of both MSE and DSE
# Will be used to compare to terms in MSE budget
# 
# - NCSU Large Scale and Tropical Dynamics
# - Based on MSE Budget vertical integration provided by A. Aiyyer
#
# Versions
# - A. Thornton, Sep 5, 2023

import numpy as np
import xarray as xr
import pandas as pd
from datetime import date
from numpy import absolute, exp, log

# Any import of metpy will activate the accessors
from metpy.units import units

#from metpy.calc import dewpoint_from_relative_humidity
from metpy.calc import first_derivative, geospatial_gradient, advection

import metpy.constants as constants


# PATH TO DATA 
path_full = '/glade/scratch/athornton/era5_processed_data/3d/' 

# SELECT SUBSET OF DATES
year_start  = 1998
month_start = 3
day_start   = 1

year_end  = 2022
month_end = 12
day_end   = 31

date_series = [pd.date_range(date(i,month_start,day_start),date(i,month_end,day_end), freq ='D') for i in range(year_start,year_end+1)]

# date_series is a list of lists. Lets unpack it now
dates_list = [element for sublist in date_series for element in sublist]
print(dates_list[0].strftime("%Y%m%d"))
print(dates_list[-1].strftime("%Y%m%d"))


# mass weighted vertical integral of a quantity
def mass_weighted_vert_integral(data):
    # data is expected to be on pressure levels
    levels = (data.level*units(data.level.units)).metpy.convert_units('Pa')
    deltaP = (levels - levels.shift(level=1)).metpy.convert_units('Pa')
    vert_int_data = ((data.shift(level=1)+data)*.5*deltaP).sum(dim='level') / constants.earth_gravity   
    return  vert_int_data


def print_minmax(var_str, data):
    print( var_str , ' min, max = ', data.min().values, data.max().values, data.metpy.units )
    return

for a_date in dates_list:
    print (a_date)
    # Term 1
    
    # get the 4x daily mse
    path_data = path_full + 'mse_' + a_date.strftime("%Y%m%d") + '.nc'
    ds = xr.open_dataset(path_data)
    h = ds.MSE.metpy.convert_units('joule/kilogram')
    
    # get the 4x daily dse
    path_data = path_full + 'dse_' + a_date.strftime("%Y%m%d") + '.nc'
    ds = xr.open_dataset(path_data)
    s = ds.DSE.metpy.convert_units('joule/kilogram')
    
    #----------------------------------------------------------------------------------------   
    # vertically integrated 
    path = '/glade/scratch/athornton/era5_processed_data/budget_terms/unfiltered_terms/'
    
    mse_vint = mass_weighted_vert_integral(h)
    mse_vint.name = 'mse_vint'
    file_out = path + 'mse_vint' + '_unfiltered_' + a_date.strftime("%Y%m%d") + '.nc'
    mse_vint.metpy.dequantify().to_netcdf(path=file_out, format='NETCDF4', mode='w')
    
    dse_vint = mass_weighted_vert_integral(s)
    dse_vint.name = 'dse_vint'    
    file_out = path + 'dse_vint' + '_unfiltered_' + a_date.strftime("%Y%m%d") + '.nc'
    dse_vint.metpy.dequantify().to_netcdf(path=file_out, format='NETCDF4', mode='w')
    
    print_minmax('dse_vint', dse_vint)    
  






