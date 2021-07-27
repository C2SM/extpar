#! /usr/bin/env python

import numpy as np
import netCDF4 as nc

#_________________________________________________________________________________________________
#
# read the cdo processed data and retrieve all attributes needed

ei_oro = nc.Dataset("ea-oro.nc", "r")
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

extpar_ei_t2m = extpar_ei.createVariable('T_2M_CLIM', np.float32, ('time', 'ke', 'je', 'ie',))
extpar_ei_t2m.standard_name = '2m_air_temperature'
extpar_ei_t2m.long_name = '2m air temperature'
extpar_ei_t2m.units = 'K'
extpar_ei_t2m[:,:,:,:] = np.reshape(ei.variables['T_2M_CLIM'][:,:], (12, 1, 1, cells))

extpar_ei_topo = extpar_ei.createVariable('TOPO_CLIM', np.float32, ('ke', 'je', 'ie',))
extpar_ei_topo.standard_name = 'height_above_reference_ellipsoid'
extpar_ei_topo.long_name = 'Geometrical height'
extpar_ei_topo.units = 'm'
extpar_ei_topo[:,:,:] = np.reshape(ei_oro.variables['TOPO_CLIM'][:,:], (1, 1, cells))

extpar_ei.close()
