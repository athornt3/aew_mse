#!/usr/bin/env python
# coding: utf-8

# Calculate the terms of the mse budget
# A. Aiyyer & A. Thornton
# 
# Aug 2023

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


# path to regridded era5 4x daily data for variables
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

# We now calculate the various terms term: e.g., advection of the anomalous MSE by the mean zonal wind

# read the basic state MSE
varNam = 'mse'
variab = 'MSE'

# path to averages
path_bar = '/glade/scratch/athornton/era5_processed_data/jas_means/'+varNam+'_jas_1998_2022.nc'
ds = xr.open_dataset(path_bar)

# the basic state MSE [calculated using mse.ipynb]
hbar = ds[variab].metpy.convert_units('joule/kilogram')
print_minmax('hbar', hbar)

# read the basic state Zonal Wind
varNam = 'u'
variab = 'U'
# path to averages
path_bar = '/glade/scratch/athornton/era5_processed_data/jas_means/'+varNam+'_jas_1998_2022.nc'
ds = xr.open_dataset(path_bar)
# Extract variables
ubar = ds[variab]*units['meter / second']
print_minmax('ubar', ubar)


# read the basic state meridional Wind
varNam = 'v'
variab = 'V'
# path to averages
path_bar = '/glade/scratch/athornton/era5_processed_data/jas_means/'+varNam+'_jas_1998_2022.nc'
ds = xr.open_dataset(path_bar)
# Extract variables
vbar = ds[variab]*units['meter / second']
print_minmax('vbar', vbar)


# read the basic state meridional Wind
varNam = 'w'
variab = 'W'
# path to averages
path_bar = '/glade/scratch/athornton/era5_processed_data/jas_means/'+varNam+'_jas_1998_2022.nc'
ds = xr.open_dataset(path_bar)
# Extract variables
wbar = ds[variab]*units['pascal / second']
print_minmax('wbar', wbar)


for a_date in dates_list:
    print (a_date)
    # get the 4x daily mse
    path_data = path_full + 'mse_' + a_date.strftime("%Y%m%d") + '.nc'
    ds = xr.open_dataset(path_data)
    h = ds.MSE.metpy.convert_units('joule/kilogram')
    
    # get the 4x daily zonal wind
    path_data = path_full + 'u_' + a_date.strftime("%Y%m%d") + '.nc'
    ds = xr.open_dataset(path_data)
    u = ds.U.metpy.quantify()
    
    
    # get the 4x daily zonal wind
    path_data = path_full + 'v_' + a_date.strftime("%Y%m%d") + '.nc'
    ds = xr.open_dataset(path_data)
    v = ds.V.metpy.quantify()
    
    
    # calculate the anomaly fields
    hp = h - hbar
    print_minmax('mse anomaly', hp)
   

    # calculate the anomaly fields
    up = u - ubar
    print_minmax('u anomaly', up)
   
    vp = v - vbar
    print_minmax('v anomaly', vp)
 
    # calculate the gradient of hp
    hp_dx, hp_dy = geospatial_gradient(hp)
    # Create  new xarray DataArrays with the calculated gradients
    hp_dx = xr.DataArray(hp_dx, coords=hp.coords, dims=hp.dims)
    hp_dy = xr.DataArray(hp_dy, coords=hp.coords, dims=hp.dims)

    # calculate the gradient of the mean MSE (hbar)
    hbar_dx, hbar_dy = geospatial_gradient(hbar)
    # Create new xarray DataArrays with the calculated gradients
    hbar_dx = xr.DataArray(hbar_dx, coords=hbar.coords, dims=hbar.dims)
    hbar_dy = xr.DataArray(hbar_dy, coords=hbar.coords, dims=hbar.dims)
  
    
    #----------------------------------------------------------------------------------------
    # terms = ubar dh'/dx  & ar dh'/dy    
    # vertically integrated 
    path = '/glade/scratch/athornton/era5_processed_data/budget_terms/unfiltered_terms/'
    
    ubar_hp_dx_vint = mass_weighted_vert_integral(ubar*hp_dx).metpy.convert_units('watt/meter**2')
    ubar_hp_dx_vint.name = 'ubar_hp_dx_vint'
    file_out = path + 'ubar_hp_dx_vint' + '_unfiltered_' + a_date.strftime("%Y%m%d") + '.nc'
    ubar_hp_dx_vint.metpy.dequantify().to_netcdf(path=file_out, format='NETCDF4', mode='w')
    
    vbar_hp_dy_vint = mass_weighted_vert_integral(vbar*hp_dy).metpy.convert_units('watt/meter**2')
    vbar_hp_dy_vint.name = 'vbar_hp_dy_vint'    
    file_out = path + 'vbar_hp_dy_vint' + '_unfiltered_' + a_date.strftime("%Y%m%d") + '.nc'
    vbar_hp_dy_vint.metpy.dequantify().to_netcdf(path=file_out, format='NETCDF4', mode='w')
    
    print_minmax('ubar_hp_dx_vint', ubar_hp_dx_vint)

    
    #----------------------------------------------------------------------------------------
    # terms = u' d/dx (hbar) &    v' d/dy (hbar)
    

    up_hbar_dx_vint = mass_weighted_vert_integral(up*hbar_dx).metpy.convert_units('watt/meter**2')
    up_hbar_dx_vint.name = 'up_hbar_dx_vint'    
    file_out = path + 'up_hbar_dx_vint' + '_unfiltered_' + a_date.strftime("%Y%m%d") + '.nc'    
    up_hbar_dx_vint.metpy.dequantify().to_netcdf(path=file_out, format='NETCDF4', mode='w')

    vp_hbar_dy_vint = mass_weighted_vert_integral(vp*hbar_dy).metpy.convert_units('watt/meter**2')
    vp_hbar_dy_vint.name = 'vp_hbar_dy_vint'    
    file_out = path + 'vp_hbar_dy_vint' + '_unfiltered_' + a_date.strftime("%Y%m%d") + '.nc'    
    vp_hbar_dy_vint.metpy.dequantify().to_netcdf(path=file_out, format='NETCDF4', mode='w')
   
    print_minmax('up_hbar_dx_vint',up_hbar_dx_vint)

    
    #----------------------------------------------------------------------------------------
    # terms = u' dh'/dx    & v' dh'/dy    
    up_hp_dx_vint = mass_weighted_vert_integral(up*hp_dx).metpy.convert_units('watt/meter**2')
    up_hp_dx_vint.name = 'up_hp_dx_vint'    
    file_out = path + 'up_hp_dx_vint' + '_unfiltered_' + a_date.strftime("%Y%m%d") + '.nc'    
    up_hp_dx_vint.metpy.dequantify().to_netcdf(path=file_out, format='NETCDF4', mode='w')
    
    vp_hp_dy_vint = mass_weighted_vert_integral(vp*hp_dy).metpy.convert_units('watt/meter**2')
    vp_hp_dy_vint.name = 'vp_hp_dy_vint'    
    file_out = path + 'vp_hp_dy_vint' + '_unfiltered_' + a_date.strftime("%Y%m%d") + '.nc'    
    vp_hp_dy_vint.metpy.dequantify().to_netcdf(path=file_out, format='NETCDF4', mode='w')
    
    print_minmax('up_hp_dx_vint',up_hp_dx_vint)
  
    #----------------------------------------------------------------------------------------
    # vertical advection terms
    
    # get the 4x daily vertical wind
    path_data = path_full + 'w_' + a_date.strftime("%Y%m%d") + '.nc'
    ds = xr.open_dataset(path_data)
    w = ds.W.metpy.quantify()
    print_minmax('w',w)

    wp = w - wbar
    print_minmax('wp',wp)
 
    wp_hp_dp_vint = mass_weighted_vert_integral( advection(hp, w=wp)).metpy.convert_units('watt/meter**2')
    wp_hp_dp_vint.name = 'wp_hp_dp_vint'    
    file_out = path + 'wp_hp_dp_vint' + '_unfiltered_' + a_date.strftime("%Y%m%d") + '.nc'    
    wp_hp_dp_vint.metpy.dequantify().to_netcdf(path=file_out, format='NETCDF4', mode='w')
    print_minmax('wp_hp_dp_vint',wp_hp_dp_vint)
  
    wbar_hp_dp_vint = mass_weighted_vert_integral( advection(hp, w=wbar)).metpy.convert_units('watt/meter**2')
    wbar_hp_dp_vint.name = 'wbar_hp_dp_vint'    
    file_out = path + 'wbar_hp_dp_vint' + '_unfiltered_' + a_date.strftime("%Y%m%d") + '.nc'    
    wbar_hp_dp_vint.metpy.dequantify().to_netcdf(path=file_out, format='NETCDF4', mode='w')    
    print_minmax('wbar_hp_dp_vint',wbar_hp_dp_vint)
  
    wp_hbar_dp_vint = mass_weighted_vert_integral( advection(hbar, w=wp)).metpy.convert_units('watt/meter**2')
    wp_hbar_dp_vint.name = 'wp_hbar_dp_vint'    
    file_out = path + 'wp_hbar_dp_vint' + '_unfiltered_' + a_date.strftime("%Y%m%d") + '.nc'    
    wp_hbar_dp_vint.metpy.dequantify().to_netcdf(path=file_out, format='NETCDF4', mode='w')    
    print_minmax('wp_hbar_dp_vint',wp_hbar_dp_vint)
  






