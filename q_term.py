#!/usr/bin/env python
# coding: utf-8

# Code to sum era5 fluxes to produce the total heating term in the MSE budget equation
# Manipulated to represent atmospheric gains
# 
# - NCSU Large Scale and Tropical Dynamics
# 
# Versions
# - A. Thornton, Aug 26, 2023

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


# regridded era5 4x daily data for variables
path_full = '/glade/scratch/athornton/era5_processed_data/2d/2d_updated' 
path_out = '/glade/scratch/athornton/era5_processed_data/budget_terms/unfiltered_terms/'

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


# for each date 
for a_date in dates_list:
    print(a_date)
    # ECMWF defines positive as downward
    # Surface fluxes are atmospheric losses (multiply by negative to fix)
    path_data = path_full + 'ssr_' + a_date.strftime('%Y%m%d') + '.nc'
    ds = xr.open_dataset(path_data)
    ssr = ds.SSR.metpy.quantify()*(-1)
    ds.close()
    
    path_data = path_full + 'str_' + a_date.strftime('%Y%m%d') + '.nc'
    ds = xr.open_dataset(path_data)
    Str = ds.STR.metpy.quantify()*(-1)
    ds.close()
    
    path_data = path_full + 'tsr_' + a_date.strftime('%Y%m%d') + '.nc'
    ds = xr.open_dataset(path_data)
    tsr = ds.TSR.metpy.quantify()
    ds.close()
    
    path_data = path_full + 'ttr_' + a_date.strftime('%Y%m%d') + '.nc'
    ds = xr.open_dataset(path_data)
    ttr = ds.TTR.metpy.quantify()
    ds.close()
    
    Q = ssr + Str + tsr + ttr
    Q.name = 'Q'
    
    file_out = path_out + 'q_unfiltered_' + a_date.strftime("%Y%m%d") + '.nc'
    Q.metpy.dequantify().to_netcdf(path=file_out,  format='NETCDF4', mode='w')






