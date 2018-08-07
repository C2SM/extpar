#! /usr/bin/env python

import numpy as np
import netCDF4 as nc

#_________________________________________________________________________________________________
#
# read the cdo processed data and retrieve all attributes needed

tcl = nc.Dataset("t_cl-dis.nc", "r")
cells = len(tcl.dimensions['cell'])

#_________________________________________________________________________________________________
#
# create extpar BUFFER

extpar_tcl = nc.Dataset("t_cl-dis_BUFFER.nc", "w", format='NETCDF4')

extpar_tcl.createDimension('ie', cells)
extpar_tcl.createDimension('je', 1)
extpar_tcl.createDimension('ke', 1)

extpar_lon  = extpar_tcl.createVariable('lon', np.float32, ('ke', 'je', 'ie',))
extpar_lon.standard_name = 'longitude'
extpar_lon.long_name = 'geographical longitude'
extpar_lon.units = 'degrees_north'
extpar_lon[:,:,:] =  np.rad2deg(np.reshape(tcl.variables['clon'][:], (1,1,cells)))

extpar_lat  = extpar_tcl.createVariable('lat', np.float32, ('ke', 'je', 'ie',))
extpar_lat.standard_name = 'latitude'
extpar_lat.long_name = 'geographical latitude'
extpar_lat.units = 'degrees_east'
extpar_lat[:,:,:] = np.rad2deg(np.reshape(tcl.variables['clat'][:], (1,1,cells)))

extpar_hsurf = extpar_tcl.createVariable('HSURF', np.float32, ('ke', 'je', 'ie',))
extpar_hsurf.standard_name = 'surface_altitude'
extpar_hsurf.long_name = 'CRU grid elevation'
extpar_hsurf.units = 'm'
extpar_hsurf[:,:,:] = np.reshape(tcl.variables['HSURF'][:], (1, 1, cells))

extpar_t_cl = extpar_tcl.createVariable('T_CL', np.float32, ('ke', 'je', 'ie',))
extpar_t_cl.standard_name = 'soil_temperature'
extpar_t_cl.long_name = 'CRU near surface temperature climatology'
extpar_t_cl.units = 'K'
extpar_t_cl[:,:,:] = np.reshape(tcl.variables['T_CL'][:], (1, 1, cells))

extpar_tcl.close()
