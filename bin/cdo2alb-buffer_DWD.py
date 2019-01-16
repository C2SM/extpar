#! /usr/bin/env python

import numpy as np
import netCDF4 as nc

#_________________________________________________________________________________________________
#
# read the cdo processed data and retrieve all attributes needed

alb = nc.Dataset("alb-dis.nc", "r")
cells = len(alb.dimensions['cell'])

albni = nc.Dataset("alnid-dis.nc", "r")
albuv = nc.Dataset("aluvd-dis.nc", "r")

#_________________________________________________________________________________________________
#
# create extpar BUFFER

extpar_alb = nc.Dataset("alb-dis_BUFFER.nc", "w", format='NETCDF4')

extpar_alb.createDimension('ie', cells)
extpar_alb.createDimension('je', 1)
extpar_alb.createDimension('ke', 1)
extpar_alb.createDimension('time', None)
extpar_alb.createDimension('mlev', 1)

extpar_time = extpar_alb.createVariable('time', np.float32, ('time',))
extpar_time[:] = np.array([11110101, 11110201, 11110301, 11110401, 11110501, 11110601,
                           11110701, 11110801, 11110901, 11111001,  11111101, 11111201],
                          dtype=np.float32)

extpar_mlev = extpar_alb.createVariable('mlev', np.float32, ('mlev',))
extpar_mlev[:] = np.array([1], dtype=np.float32)

extpar_lon  = extpar_alb.createVariable('lon', np.float32, ('ke', 'je', 'ie',))
extpar_lon.standard_name = 'longitude'
extpar_lon.long_name = 'geographical longitude'
extpar_lon.units = 'degrees_north'
extpar_lon[:,:,:] =  np.rad2deg(np.reshape(alb.variables['clon'][:], (1,1,cells)))

extpar_lat  = extpar_alb.createVariable('lat', np.float32, ('ke', 'je', 'ie',))
extpar_lat.standard_name = 'latitude'
extpar_lat.long_name = 'geographical latitude'
extpar_lat.units = 'degrees_east'
extpar_lat[:,:,:] = np.rad2deg(np.reshape(alb.variables['clat'][:], (1,1,cells)))

extpar_alb_dif12 = extpar_alb.createVariable('ALB_DIF12', np.float32, ('time', 'ke', 'je', 'ie',))
extpar_alb_dif12.standard_name = 'Albedo'
extpar_alb_dif12.long_name = 'Albedo'
extpar_alb_dif12.units = '%'
extpar_alb_dif12[:,:,:,:] = np.reshape(alb.variables['al'][:,:], (12, 1, 1, cells))

extpar_alnid12 = extpar_alb.createVariable('ALNID12', np.float32, ('time', 'ke', 'je', 'ie',))
extpar_alnid12.standard_name = 'NI_Albedo'
extpar_alnid12.long_name ='NI_Albedo'
extpar_alnid12.units = '%'
extpar_alnid12[:,:,:,:] = np.reshape(albni.variables['alnid'][:,:], (12, 1, 1, cells))

extpar_aluvd12 = extpar_alb.createVariable('ALUVD12', np.float32, ('time', 'ke', 'je', 'ie',))
extpar_aluvd12.standard_name = 'UV_Albedo'
extpar_aluvd12.long_name = 'UV_Albedo'
extpar_aluvd12.units = '%'
extpar_aluvd12[:,:,:,:] = np.reshape(albuv.variables['aluvd'][:,:], (12, 1, 1, cells))

extpar_alb.close()
