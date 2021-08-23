import logging
import math
import numpy as np
import netCDF4 as nc

import utilities as utils
from fortran_namelist import read_variable

'''
Module providing classes and functions for target grids,
it contains:

    -class CosmoGrid: grid definition from namelist

    -class IconGrid: grid information from NetCDF Icon-grid
'''


class CosmoGrid:
    '''
    store all Cosmo grid information and provide functions related to

    member functions:

        -__init__

        -create_grid_description

        -lon_rot

        -lat_rot

        -latlon_cosmo_to_latlon_regular

        -rlarot2rla

        -phirot2phi

        -cdo_sellonlat
    '''

    def __init__(self, namelist):
        '''
        init grid from existing Fortran namelist 'namelist'

        the return value of function "read_variable" 
        is converted to the right type (int, float) 
        No check if retrieved values are meaningful is done
        '''

        self.pollon = read_variable(namelist, 'pollon', float)
        self.pollat = read_variable(namelist, 'pollat', float)

        self.dlon = read_variable(namelist, 'dlon', float)
        self.dlat = read_variable(namelist, 'dlat', float)

        self.ie_tot = read_variable(namelist, 'ie_tot', int)
        self.je_tot = read_variable(namelist, 'je_tot', int)

        self.startlon_tot = \
            read_variable(namelist, 'startlon_tot', float)

        self.startlat_tot = \
            read_variable(namelist, 'startlat_tot', float)

        # infer from existing values
        self.ke_tot = 1
        self.gridsize = self.ie_tot * self.je_tot

        self.lats,self.lons = self.latlon_cosmo_to_latlon_regular()

    def create_grid_description(self,name):
        '''
        write grid description required for cdo

        grid decription has name "name"
        it contains all required information
        about the grid needed for the interpolation using CDO
        '''

        with open(name,'w') as f:

            f.write(f'gridtype  = projection\n')
            f.write(f'gridsize  = {self.gridsize}\n')
            f.write(f'xsize     = {self.ie_tot}\n')
            f.write(f'ysize     = {self.je_tot}\n')
            f.write(f'xinc      = {self.dlon}\n')
            f.write(f'yinc      = {self.dlat}\n')
            f.write(f'xfirst    = {self.startlon_tot}\n')
            f.write(f'yfirst    = {self.startlat_tot}\n')

            f.write(f'grid_mapping_name = rotated_latitude_longitude\n')
            f.write(f'grid_north_pole_longitude = {self.pollon}\n') 
            f.write(f'grid_north_pole_latitude  = {self.pollat}\n') 

    def lon_rot(self):

        '''return array with rotated longitude values'''

        lon = np.empty(self.ie_tot)

        for i in range(self.ie_tot):
            lon[i] = self.startlon_tot + i * self.dlon

        return lon

    def lat_rot(self):

        '''return array with rotated latitude values'''

        lat = np.empty([self.je_tot])

        for j in range(self.je_tot):
            lat[j] = self.startlat_tot + j * self.dlat

        return lat

    def latlon_cosmo_to_latlon_regular(self):
        '''
        return array(je_tot,ie_tot) with the regular lat/lon values

        CAUTION: lon and lat geo are stored in j,i-order, because the netCDF
        requires this index-order
        '''
        lat_cosmo = self.lat_rot()
        lon_cosmo = self.lon_rot()

        lat_reg = np.empty([self.je_tot,self.ie_tot])
        lon_reg = np.empty([self.je_tot,self.ie_tot])

        for j in range(self.je_tot):
            for i in range(self.ie_tot):

                lon_reg[j,i] = self.rlarot2rla(lat_cosmo[j], 
                                               lon_cosmo[i], 
                                               self.pollat, 
                                               self.pollon)

                lat_reg[j,i] = self.phirot2phi(lat_cosmo[j], 
                                               lon_cosmo[i], 
                                               self.pollat, 
                                               self.pollon)

        return lat_reg,lon_reg

    def rlarot2rla(self,phirot, rlarot, polphi, pollam): 
        '''
        convert rotated longitude to regular longitude

        functions taken from Fortran module "mo_utilities_extpar.f90"
        results of function are cross-validated 
        with cartopy coordinate transformation
        '''

        zrpi18 = 57.29577951308232
        zpir18 = 0.017453292519943295

        zsinpol = math.sin(zpir18 * polphi)
        zcospol = math.cos(zpir18 * polphi)

        zlampol = zpir18 * pollam
        zphis   = zpir18 * phirot

        if (rlarot > 180.0): 
            zrlas = rlarot - 360.0
        else:
            zrlas = rlarot

        zrlas = zpir18 * zrlas

        zarg1   = (math.sin(zlampol) * (-zsinpol * math.cos(zrlas) * 
                                        math.cos(zphis) +
                                        zcospol * math.sin(zphis)) -
                   math.cos(zlampol) * math.sin(zrlas) * math.cos(zphis))

        zarg2   = (math.cos(zlampol) * (-zsinpol * math.cos(zrlas) * 
                                        math.cos(zphis) +
                                        zcospol * math.sin(zphis)) +
                   math.sin(zlampol) * math.sin(zrlas) * math.cos(zphis))

        if (zarg2 == 0.0): 
            zarg2 = 1.0E-20

        rla = zrpi18 * math.atan2(zarg1,zarg2)

        return rla

    def phirot2phi(self, phirot, rlarot, polphi, pollam): 
        '''
        convert rotated latitude to regular latitude

        functions taken from Fortran module "mo_utilities_extpar.f90"
        results of function are cross-validated 
        with cartopy coordinate transformation
        '''

        zrpi18 = 57.29577951308232
        zpir18 = 0.017453292519943295

        zsinpol     = math.sin(zpir18 * polphi)
        zcospol     = math.cos(zpir18 * polphi)

        zphis       = zpir18 * phirot

        if (rlarot > 180.0):
            zrlas = rlarot - 360.0
        else:
            zrlas = rlarot

        zrlas       = zpir18 * zrlas

        zarg  = (zcospol * math.cos(zphis) * math.cos(zrlas) +
                 zsinpol * math.sin(zphis))

        phi = zrpi18 * math.asin(zarg)

        return phi

    def cdo_sellonlat(self):
        '''
        create the string for the cdo option "-sellonlatbox"

        determine lat/lon extent from grid information and add
        a +/- 1 degree safety margin. Tested for regional and global
        grids
        '''

        extent = {}
        extent['maxlon'] = min( math.ceil(np.amax(self.lons) + 1), 180.0)
        extent['minlon'] = max( math.floor(np.amin(self.lons) - 1), -180.0)

        extent['maxlat'] = min( math.ceil(np.amax(self.lats) + 1), 90.0)
        extent['minlat'] = max(math.floor(np.amin(self.lats) - 1), -90.0)

        cdo_selbox = ('-sellonlatbox,' 
                      f"{extent['minlon']},"
                      f"{extent['maxlon']},"
                      f"{extent['minlat']},"
                      f"{extent['maxlat']}")

        return cdo_selbox


class IconGrid:
    '''
    store relevant Icon grid information and 
    provide functions related to

    member functions:

        -__init__

        -cdo_sellonlat

        -reduce_grid
    '''

    def __init__(self,gridfile):

        self.gridfile = gridfile
        self.grid = nc.Dataset(gridfile, "r")

        self.lons = np.rad2deg(self.grid.variables["clon"][:]) 
        self.lats = np.rad2deg(self.grid.variables["clat"][:])

    def cdo_sellonlat(self):
        '''
        create the string for the cdo option "-sellonlatbox"

        determine lat/lon extent from grid information and add
        a +/- 1 degree safety margin. Tested for regional and global
        grids
        '''

        extent = {}
        extent['maxlon'] = min( math.ceil(np.amax(self.lons) + 1), 180.0)
        extent['minlon'] = max( math.floor(np.amin(self.lons) - 1), -180.0)

        extent['maxlat'] = min( math.ceil(np.amax(self.lats) + 1), 90.0)
        extent['minlat'] = max(math.floor(np.amin(self.lats) - 1), -90.0)

        cdo_selbox = ('-sellonlatbox,' 
                      f"{extent['minlon']},"
                      f"{extent['maxlon']},"
                      f"{extent['minlat']},"
                      f"{extent['maxlat']}")

        return cdo_selbox

    def reduce_grid(self, reduced_grid):
        '''
        reduce icon grid and return name of reduce_grid

        due to inconsistency in the order of dimensions in icon grids,
        only use a reduced subset of the respective icon grid to avoid
        problems with dimension names later in the code
        '''

        logging.info(f'Reduce icon grid {self.gridfile} '
                     f'-> {reduced_grid} to avoid wrong grid usage '
                     'during CDO interpolation\n')

        utils.launch_shell('cdo', '-f', 'nc4',
                           '-selvar,cell_area,clat,clat_vertices,'
                           'clon,clon_vertices',
                           self.gridfile,
                           reduced_grid)

        return reduced_grid
