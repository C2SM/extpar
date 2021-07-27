#! /usr/bin/env python

import numpy as np
import netCDF4 as nc

#_________________________________________________________________________________________________
#
# read the cdo processed data and retrieve all attributes needed

ei = nc.Dataset("ea.nc", "r")
cells = len(ei.dimensions['cell'])

#_________________________________________________________________________________________________
#
# create extpar BUFFER

extpar_ei = nc.Dataset("ea_BUFFER.nc", "w", format='NETCDF4')

extpar_ei.createDimension('ie', cells)
extpar_ei.createDimension('je', 1)
extpar_ei.createDimension('ke', 1)
extpar_ei.createDimension('time', None)

extpar_time = extpar_ei.createVariable('time', np.float64, ('time'))
extpar_time.standard_name = 'time'
extpar_time.units = 'hours since 1111-1-11 00:00:00'
extpar_time.calendar = 'proleptic_gregorian'
extpar_time[:] = ei.variables['time'][:]

extpar_lon  = extpar_ei.createVariable('lon', np.float64, ('ke', 'je', 'ie',))
extpar_lon.standard_name = 'longitude'
extpar_lon.long_name = 'geographical longitude'
extpar_lon.units = 'degrees_north'
extpar_lon[:,:,:] =  np.rad2deg(np.reshape(ei.variables['clon'][:], (1,1,cells)))

extpar_lat  = extpar_ei.createVariable('lat', np.float64, ('ke', 'je', 'ie',))
extpar_lat.standard_name = 'latitude'
extpar_lat.long_name = 'geographical latitude'
extpar_lat.units = 'degrees_east'
extpar_lat[:,:,:] = np.rad2deg(np.reshape(ei.variables['clat'][:], (1,1,cells)))

extpar_ei_t_s = extpar_ei.createVariable('T_S', np.float32, ('time', 'ke', 'je', 'ie',))
extpar_ei_t_s.standard_name = 'surface_skin_temperature'
extpar_ei_t_s.long_name = 'Skin temperature'
extpar_ei_t_s.units = 'K'
extpar_ei_t_s[:,:,:,:] = np.reshape(ei.variables['T_S'][:,:], (12, 1, 1, cells))

extpar_ei_w_snow = extpar_ei.createVariable('W_SNOW', np.float32, ('time', 'ke', 'je', 'ie',))
extpar_ei_w_snow.standard_name = 'lwe_thickness_of_surface_snow_amount'
extpar_ei_w_snow.long_name ='Snow depth'
extpar_ei_w_snow.units = 'm of water equivalent'
extpar_ei_w_snow[:,:,:,:] = np.reshape(ei.variables['W_SNOW'][:,:], (12, 1, 1, cells))

extpar_ei_t_sea = extpar_ei.createVariable('T_SEA', np.float32, ('time', 'ke', 'je', 'ie',))
extpar_ei_t_sea.standard_name = 'sea_surface_temperature'
extpar_ei_t_sea.long_name = 'Sea surface temperature'
extpar_ei_t_sea.units = 'K'
extpar_ei_t_sea[:,:,:,:] = np.reshape(ei.variables['T_SEA'][:,:], (12, 1, 1, cells))

extpar_ei.close()
