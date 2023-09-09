#!/usr/bin/env python
# coding: utf-8

# Reads the regridded era5 files and calculates the MSE & DSE 

# NCSU Tropical and Large Scale Dynamics
# A. Aiyyer Jul 2023
# A. Thornton Aug 2023


import numpy as np
import xarray as xr
import pandas as pd
from datetime import date
from numpy import absolute, exp, log

# Any import of metpy will activate the accessors
from metpy.units import units
from metpy.calc import moist_static_energy, dry_static_energy
import metpy.constants

# path for the regridded data
path_data = '/glade/scratch/athornton/era5_processed_data/3d/'

# output path for mse and dse - same as rest of the 3d variables? 
path_out = '/glade/scratch/athornton/era5_processed_data/3d/'

# SELECT SUBSET OF DATES
year_start  = 1998
month_start = 3
day_start   = 1

year_end  = 1998
month_end = 3
day_end   = 1

date_series = [pd.date_range(date(i,month_start,day_start),date(i,month_end,day_end), freq ='D') for i in range(year_start,year_end+1)]
# date_series is a list of lists. Lets unpack it now
dates_list = [element for sublist in date_series for element in sublist]

print(dates_list[0].strftime("%Y%m%d"))
print(dates_list[-1].strftime("%Y%m%d"))


for a_date in dates_list:
    print (a_date)
    # read the variables needed for MSE
    
    infile = path_data + 't_' + a_date.strftime("%Y%m%d") + '.nc'
    ds = xr.open_dataset(infile)
    tempK = ds.T.metpy.quantify()
    ds.close()
    
    infile = path_data + 'z_' + a_date.strftime("%Y%m%d") + '.nc'
    ds = xr.open_dataset(infile)
    
    # read and convert to geopotential height
    geopot = ds.Z.metpy.quantify()/metpy.constants.earth_gravity
    ds.close()
    
    infile = path_data + 'q_' + a_date.strftime("%Y%m%d") + '.nc'
    ds = xr.open_dataset(infile)
    sphum = ds.Q.metpy.quantify()
    ds.close()
    
    mse =  moist_static_energy(geopot, tempK, sphum)
    mse.name = 'MSE'
    
    dse =  dry_static_energy(geopot, tempK)
    dse.name = 'DSE'
    
    #write to file (4x daily)
    file_out = path_out + 'mse_' + a_date.strftime("%Y%m%d") + '.nc'
    
    #dequantify is needed to correctly write out units into the netcdf file
    #see: https://unidata.github.io/MetPy/latest/tutorials/xarray_tutorial.html
    mse.metpy.dequantify().to_netcdf(path=file_out,  format='NETCDF4', mode='w')
 
    #write to file (4x daily)
    file_out = path_out + 'dse_' + a_date.strftime("%Y%m%d") + '.nc'
    dse.metpy.dequantify().to_netcdf(path=file_out,  format='NETCDF4', mode='w')
 




