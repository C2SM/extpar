import numpy as np
'''
Module providing the Meta-Data classes for the buffer file,
it contains:

    -Parent: Coordinates -> Child: Lon, Lat

    -Parent: AlbMeta     -> Child: AL, NI, UV

    -Parent: NdviMeta    -> Child: NDVI, NdviMax, NdviMrat

    -Parent: EmissMeta   -> Child: EmissMean, EmissMax, EmissMrat

    -Parent: ClimMeta    -> Child: TempClim, HsurfClim

Meta-Data that is shared amongs all fields of an Extpar class is defined in
the parent class, for example CoordsMeta 
Meta-Data that is only valid for one specific field is defined 
in the respective child class

The following Meta-Data is required during the write_*d_fields functions:
    -name ( string)
    -short ( string)
    -standard ( string)
    -long ( string)
    -type ( basic datatype (float, integer,...)
    -dim ( dictionary for all dimensions of fields, ie is the last index)
'''

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# Coordinates
# ->Lon
# ->Lat


class CoordsMeta:
    def __init__(self):         
        self.type = np.float32
        self.dim = {0: 'ke',
                    1: 'je',
                    2: 'ie'}


class Lon(CoordsMeta):
    def __init__(self):         
        super().__init__()
        self.name = 'lon'
        self.short = '_'
        self.standard = 'longitude'
        self.long = 'geographical longitude'
        self.units = 'degrees_east'


class Lat(CoordsMeta):
    def __init__(self):         
        super().__init__()
        self.name = 'lat'
        self.short = '-'
        self.standard = 'latitude'
        self.long = 'geographical latitude'
        self.units = 'degrees_north'
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# Albedo
# ->AL
# ->NI
# ->UV
# ->AlbSat
# ->AlbDry


class AlbMeta:
    def __init__(self):
        self.type = np.float32
        self.units = '%'


class AL(AlbMeta):
    def __init__(self):         
        super().__init__()
        self.name = 'ALB_DIF12'
        self.long = 'Albedo'
        self.orig = 'al'
        self.dim = { 0: 'time',
                     1: 'ke', 
                     2: 'je',
                     3: 'ie'}


class NI(AlbMeta):
    def __init__(self):         
        super().__init__()
        self.name = 'ALNID12'
        self.long = 'NI_Albedo'
        self.orig = 'alnid'
        self.dim = { 0: 'time',
                     1: 'ke', 
                     2: 'je',
                     3: 'ie'}


class UV(AlbMeta):
    def __init__(self):         
        super().__init__()
        self.name = 'ALUVD12'
        self.long = 'UV_Albedo'
        self.orig = 'aluvd'
        self.dim = { 0: 'time',
                     1: 'ke', 
                     2: 'je',
                     3: 'ie'}


class AlbDry(AlbMeta):
    def __init__(self):
        super().__init__()
        self.name = 'ALB_DRY'
        self.long = 'soil albedo for dry soil'
        self.standard = 'surface_albedo'
        self.orig = '_'
        self.dim = { 0: 'ke',
                     1: 'je',
                     2: 'ie'}


class AlbSat(AlbMeta):
    def __init__(self):
        super().__init__()
        self.name = 'ALB_SAT'
        self.long = 'soil albedo for saturated soil'
        self.standard = 'surface_albedo'
        self.orig = '_'
        self.dim = { 0: 'ke',
                     1: 'je',
                     2: 'ie'}
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# NDVI
# ->NDVI
# ->NdviMax
# ->NdviMrat


class NdviMeta:
    def __init__(self):
        self.type = np.float32
        self.units = '_'
        self.standard = '_'
        self.short = '_'


class NDVI(NdviMeta):
    def __init__(self):
        super().__init__()
        self.dim = { 0: 'time',
                     1: 'ke', 
                     2: 'je',
                     3: 'ie'}
        self.name = 'NDVI'
        self.long = 'monthly mean NDVI climatology 1998-2003'


class NdviMax(NdviMeta):
    def __init__(self):
        super().__init__()
        self.dim = { 0: 'ke', 
                     1: 'je',
                     2: 'ie'}
        self.name = 'NDVI_MAX'
        self.long = 'NDVI yearly maximum for climatology 1998-2003'


class NdviMrat(NdviMeta):
    def __init__(self):
        super().__init__()
        self.dim = { 0: 'time',
                     1: 'ke', 
                     2: 'je',
                     3: 'ie'}
        self.name = 'NDVI_MRAT'
        self.long = 'monthly proportion of actual value/maximum NDVI'

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# Emiss
# ->EmissMean
# ->EmissMax
# ->EmissMrat


class EmissMeta:
    def __init__(self):
        self.type = np.float32
        self.units = '_'
        self.standard = '_'
        self.short = '_'


class EmissMean(EmissMeta):
    def __init__(self):
        super().__init__()
        self.dim = { 0: 'time',
                     1: 'ke', 
                     2: 'je',
                     3: 'ie'}
        self.name = 'EMISS'
        self.long = 'monthly mean EMISS climatology 1998-2003'


class EmissMax(EmissMeta):
    def __init__(self):
        super().__init__()
        self.dim = { 0: 'ke', 
                     1: 'je',
                     2: 'ie'}
        self.name = 'EMISS_MAX'
        self.long = 'EMISS yearly maximum for climatology 1998-2003'


class EmissMrat(EmissMeta):
    def __init__(self):
        super().__init__()
        self.dim = { 0: 'time',
                     1: 'ke', 
                     2: 'je',
                     3: 'ie'}
        self.name = 'EMISS_MRAT'
        self.long = 'monthly proportion of actual value/maximum EMISS'

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# TCLIM
# ->TempClim
# ->HsurfClim


class ClimMeta:
    def __init__(self):
        self.type = np.float32
        self.dim = { 0: 'ke', 
                     1: 'je',
                     2: 'ie'}


class TempClim(ClimMeta):
    def __init__(self):
        super().__init__()
        self.name = 'T_CL'
        self.standard = 'soil_temperature' 
        self.long = 'CRU near surface temperature climatology'
        self.units = 'K'


class HsurfClim(ClimMeta):
    def __init__(self):
        super().__init__()
        self.name = 'HSURF'
        self.standard = 'surface_altitude' 
        self.long = 'CRU grid elevation'
        self.units = 'm'

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# ERA5
# ->SstEra
# ->T2mEra
# ->OroEra
# ->SdEra

class EraMeta:
    def __init__(self):
        self.type = np.float32
        self.standard = ''

class SstEra(EraMeta):
    def __init__(self):
        super().__init__()
        self.dim = { 0: 'time', 
                     1: 'ke',
                     2: 'je',
                     3: 'ie'}
        self.name = 'T_SEA'
        self.long = 'monthly mean ERA5SST climatology 1990-2019'
        self.units = 'K'

class T2mEra(EraMeta):
    def __init__(self):
        super().__init__()
        self.dim = { 0: 'time', 
                     1: 'ke',
                     2: 'je',
                     3: 'ie'}
        self.name = 'T_2M_CLIM'
        self.long = 'monthly mean ERA52T climatology 1990-2019'
        self.units = 'K'
        
class OroEra(EraMeta):
    def __init__(self):
        super().__init__()
        self.dim = { 0: 'ke', 
                     1: 'je',
                     2: 'ie'}
        self.name = 'TOPO_CLIM'
        self.long = 'monthly mean ERA5ORO climatology 1990-2019'
        self.units = 'm'

class SdEra(EraMeta):
    def __init__(self):
        super().__init__()
        self.dim = { 0: 'time', 
                     1: 'ke',
                     2: 'je',
                     3: 'ie'}
        self.name = 'W_SNOW'
        self.long = 'monthly mean ERA5SD climatology 1990-2019'
        self.units = 'm'
