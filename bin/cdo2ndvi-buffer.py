#! /usr/bin/env python

import numpy as np
import netCDF4 as nc

#_________________________________________________________________________________________________
#
# read the cdo processed data and retrieve all attributes needed

ndvi = nc.Dataset("ndvi-ycon.nc", "r")
cells = len(ndvi.dimensions['cell'])

#_________________________________________________________________________________________________
#
# calculate additional values

# TODO: check if things are fine

# read in ndvi field as interpolated by cdo remapycon
ndvi_field = np.reshape(ndvi.variables['ndvi'][:,:], (12, 1, 1, cells))

# calculate maxval over 12 month
ndvi_field_max = np.amax(np.reshape(ndvi.variables['ndvi'][:,:], (12, 1, 1, cells)), axis=0)

# calculate ratio of ndvi/ndvi_max per month and set 'missing value' to -1
ndvi_field_mrat = np.empty((12, 1, 1, cells), dtype=np.float32)
for t in np.arange(12):
    ndvi_field_mrat[t,:,:,:] = np.divide(ndvi_field[t,:,:,:], ndvi_field_max[:,:,:], where=ndvi_field_max[:,:,:] != 0.0)
    ndvi_field_mrat[t,:,:,:] = np.where(ndvi_field_max[:,:,:] <= 0.0, -1.0, ndvi_field_mrat[t,:,:,:])

# print some statistics    
print("Diagnostics:")
print("   NDVI max")
print("   min: {0:8.5f} mean: {1:8.5f} max: {2:8.5f}".format(np.min(ndvi_field_max),
                                              np.mean(ndvi_field_max),
                                              np.max(ndvi_field_max)))
print("   NDVI")
for t in np.arange(12):
    print("   min: {0:8.5f} mean: {1:8.5f} max: {2:8.5f}".format(np.min(ndvi_field[t,:,:,:]),
                                                  np.mean(ndvi_field[t,:,:,:]),
                                                  np.max(ndvi_field[t,:,:,:])))
print("   NDVI mrat")
for t in np.arange(12):
    print("   min: {0:8.5f} mean: {1:8.5f} max: {2:8.5f}".format(np.min(ndvi_field_mrat[t,:,:,:]),
                                                  np.mean(ndvi_field_mrat[t,:,:,:]),
                                                  np.max(ndvi_field_mrat[t,:,:,:])))
    
#_________________________________________________________________________________________________
#
# create extpar BUFFER

extpar_ndvi = nc.Dataset("ndvi-ycon_BUFFER.nc", "w", format='NETCDF4')

extpar_ndvi.createDimension('ie', cells)
extpar_ndvi.createDimension('je', 1)
extpar_ndvi.createDimension('ke', 1)
extpar_ndvi.createDimension('time', None)
extpar_ndvi.createDimension('mlev', 1)

extpar_time = extpar_ndvi.createVariable('time', np.float32, ('time'))
extpar_time[:] = np.array([11110101, 11110201, 11110301, 11110401, 11110501, 11110601,
                           11110701, 11110801, 11110901, 11111001,  11111101, 11111201],
                          dtype=np.float32)

extpar_mlev = extpar_ndvi.createVariable('mlev', np.float32, ('mlev'))
extpar_mlev[:] = np.array([1], dtype=np.float32)

extpar_lon = extpar_ndvi.createVariable('lon', np.float32, ('ke', 'je', 'ie'))
extpar_lon.standard_name = 'longitude'
extpar_lon.long_name = 'geographical longitude'
extpar_lon.units = 'degrees_north'
extpar_lon[:,:,:] =  np.rad2deg(np.reshape(ndvi.variables['clon'][:], (1,1,cells)))

extpar_lat = extpar_ndvi.createVariable('lat', np.float32, ('ke', 'je', 'ie'))
extpar_lat.standard_name = 'latitude'
extpar_lat.long_name = 'geographical latitude'
extpar_lat.units = 'degrees_east'
extpar_lat[:,:,:] = np.rad2deg(np.reshape(ndvi.variables['clat'][:], (1,1,cells)))

extpar_ndvi_field_max = extpar_ndvi.createVariable('NDVI_MAX', np.float32, ('ke', 'je', 'ie'))
extpar_ndvi_field_max.standard_name = ''
extpar_ndvi_field_max.long_name = 'NDVI yearly maximum for climatology 1998-2003'
extpar_ndvi_field_max.units = '-'
extpar_ndvi_field_max[:,:,:] = ndvi_field_max 

extpar_ndvi_field = extpar_ndvi.createVariable('NDVI', np.float32, ('time', 'ke', 'je', 'ie'))
extpar_ndvi_field.standard_name = ''
extpar_ndvi_field.long_name = 'monthly mean NDVI climatology 1998-2003'
extpar_ndvi_field.units = '-'
extpar_ndvi_field[:,:,:,:] = ndvi_field

extpar_ndvi_field_mrat = extpar_ndvi.createVariable('NDVI_MRAT', np.float32, ('time', 'ke', 'je', 'ie'))
extpar_ndvi_field_mrat.standard_name = ''
extpar_ndvi_field_mrat.long_name = 'monthly proportion of actual value/maximum NDVI'
extpar_ndvi_field_mrat.units = '-'
extpar_ndvi_field_mrat[:,:,:,:] = ndvi_field_mrat

extpar_ndvi.close()
