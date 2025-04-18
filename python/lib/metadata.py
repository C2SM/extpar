import numpy as np
'''
Module providing the Meta-Data classes for the buffer file,
it contains:

    -Parent: Coordinates -> Child: Lon, Lat

    -Parent: AlbMeta     -> Child: AL, NI, UV

    -Parent: NdviMeta    -> Child: NDVI, NdviMax, NdviMrat

    -Parent: EmissMeta   -> Child: EmissMean, EmissMax, EmissMrat

    -Parent: ClimMeta    -> Child: TempClim, HsurfClim

    -Parent: EraMeta     -> Child: SstEra5, T2mEra5, OroEra5, SdEra5
                                   SstEra-I, T2mEra-I, OroEra-I
                                   SdEra-I

    -Parent: AhfMeta     -> Child: Ahf_2min, Ahf_30sec

    -Parent: IsaMeta     -> Child: Isa_30sec, Isa_10sec

    -Parent: EdgarMeta   -> Child: EdgarOC, EdgarSO2, EdgarNOx, EdgarNH3

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
        self.dim = {0: 'ke', 1: 'je', 2: 'ie'}


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
        self.dim = {0: 'time', 1: 'ke', 2: 'je', 3: 'ie'}


class NI(AlbMeta):

    def __init__(self):
        super().__init__()
        self.name = 'ALNID12'
        self.long = 'NI_Albedo'
        self.orig = 'alnid'
        self.dim = {0: 'time', 1: 'ke', 2: 'je', 3: 'ie'}


class UV(AlbMeta):

    def __init__(self):
        super().__init__()
        self.name = 'ALUVD12'
        self.long = 'UV_Albedo'
        self.orig = 'aluvd'
        self.dim = {0: 'time', 1: 'ke', 2: 'je', 3: 'ie'}


class AlbDry(AlbMeta):

    def __init__(self):
        super().__init__()
        self.name = 'ALB_DRY'
        self.long = 'soil albedo for dry soil'
        self.standard = 'surface_albedo'
        self.orig = '_'
        self.dim = {0: 'ke', 1: 'je', 2: 'ie'}


class AlbSat(AlbMeta):

    def __init__(self):
        super().__init__()
        self.name = 'ALB_SAT'
        self.long = 'soil albedo for saturated soil'
        self.standard = 'surface_albedo'
        self.orig = '_'
        self.dim = {0: 'ke', 1: 'je', 2: 'ie'}


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
        self.dim = {0: 'time', 1: 'ke', 2: 'je', 3: 'ie'}
        self.name = 'NDVI'
        self.long = 'monthly mean NDVI climatology 1998-2003'


class NdviMax(NdviMeta):

    def __init__(self):
        super().__init__()
        self.dim = {0: 'ke', 1: 'je', 2: 'ie'}
        self.name = 'NDVI_MAX'
        self.long = 'NDVI yearly maximum for climatology 1998-2003'


class NdviMrat(NdviMeta):

    def __init__(self):
        super().__init__()
        self.dim = {0: 'time', 1: 'ke', 2: 'je', 3: 'ie'}
        self.name = 'NDVI_MRAT'
        self.long = 'monthly proportion of actual value/maximum NDVI'


#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# EDGAR
# ->EdgarBC
# ->EdgarOC
# ->EdgarSO2
# ->EdgarNOx
# ->EdgarNH3


class EdgarMeta:

    def __init__(self):
        self.type = np.float32
        self.units = 'kg m-2 s-1'
        self.standard = '_'
        self.short = '_'


class EdgarBC(EdgarMeta):

    def __init__(self):
        super().__init__()
        self.dim = {0: 'ke', 1: 'je', 2: 'ie'}
        self.name = 'emi_bc'
        self.long = 'Black Carbon for year 2022. Source: European Commission, Joint Research Centre (JRC)/Netherlands Environmental Assessment Agency (PBL). Emission Databasefor Global Atmospheric Research (EDGAR), http://edgar.jrc.ec.europe.eu'


class EdgarOC(EdgarMeta):

    def __init__(self):
        super().__init__()
        self.dim = {0: 'ke', 1: 'je', 2: 'ie'}
        self.name = 'emi_oc'
        self.long = 'Organic Carbon for year 2022. Source: European Commission, Joint Research Centre (JRC)/Netherlands Environmental Assessment Agency (PBL). Emission Databasefor Global Atmospheric Research (EDGAR), http://edgar.jrc.ec.europe.eu'


class EdgarSO2(EdgarMeta):

    def __init__(self):
        super().__init__()
        self.dim = {0: 'ke', 1: 'je', 2: 'ie'}
        self.name = 'emi_so2'
        self.long = 'Sulfur Dioxide for year 2022. Source: European Commission, Joint Research Centre (JRC)/Netherlands Environmental Assessment Agency (PBL). Emission Databasefor Global Atmospheric Research (EDGAR), http://edgar.jrc.ec.europe.eu'


class EdgarNOx(EdgarMeta):

    def __init__(self):
        super().__init__()
        self.dim = {0: 'ke', 1: 'je', 2: 'ie'}
        self.name = 'emi_nox'
        self.long = 'Nitrogen Oxides for year 2022. Source: European Commission, Joint Research Centre (JRC)/Netherlands Environmental Assessment Agency (PBL). Emission Databasefor Global Atmospheric Research (EDGAR), http://edgar.jrc.ec.europe.eu'


class EdgarNH3(EdgarMeta):

    def __init__(self):
        super().__init__()
        self.dim = {0: 'ke', 1: 'je', 2: 'ie'}
        self.name = 'emi_nh3'
        self.long = 'Ammonia for year 2022. Source: European Commission, Joint Research Centre (JRC)/Netherlands Environmental Assessment Agency (PBL). Emission Databasefor Global Atmospheric Research (EDGAR), http://edgar.jrc.ec.europe.eu'


#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# cdnc


class CdncMeta:

    def __init__(self):
        self.type = np.float32
        self.units = 'cm-3'
        self.standard = '_'
        self.short = '_'


class Cdnc(CdncMeta):

    def __init__(self):
        super().__init__()
        self.dim = {0: 'time', 1: 'ke', 2: 'je', 3: 'ie'}
        self.name = 'cdnc'
        self.long = 'cloud droplet number density (characteristic value for atmospheric column). Source: National Aeronautics and Space Administration (NASA). MODerate resolution Imaging Spectroradiometer (MODIS), https://modis.gsfc.nasa.gov/data/'


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
        self.dim = {0: 'time', 1: 'ke', 2: 'je', 3: 'ie'}
        self.name = 'EMISS'
        self.long = 'monthly mean EMISS climatology 1998-2003'


class EmissMax(EmissMeta):

    def __init__(self):
        super().__init__()
        self.dim = {0: 'ke', 1: 'je', 2: 'ie'}
        self.name = 'EMISS_MAX'
        self.long = 'EMISS yearly maximum for climatology 1998-2003'


class EmissMrat(EmissMeta):

    def __init__(self):
        super().__init__()
        self.dim = {0: 'time', 1: 'ke', 2: 'je', 3: 'ie'}
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
        self.dim = {0: 'ke', 1: 'je', 2: 'ie'}


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
# ERA
# ->SstEra5
# ->T2mEra5
# ->OroEra5
# ->SdEra5

# ->SstEraI
# ->T2mEraI
# ->OroEraI
# ->SdEraI


class EraMeta:

    def __init__(self):
        self.type = np.float32
        self.standard = ''


class SstEra5(EraMeta):

    def __init__(self):
        super().__init__()
        self.dim = {0: 'time', 1: 'ke', 2: 'je', 3: 'ie'}
        self.name = 'T_SEA'
        self.long = 'Temperature of sea water near the surface (sst) \
            from monthly mean ERA5 climatology 1990-2019'

        self.units = 'K'


class T2mEra5(EraMeta):

    def __init__(self):
        super().__init__()
        self.dim = {0: 'time', 1: 'ke', 2: 'je', 3: 'ie'}
        self.name = 'T_2M_CLIM'
        self.long = 'Temperature of air at 2m above the surface of land, \
            sea or in-land waters (2t) from monthly mean ERA5 \
            climatology 1990-2019'

        self.units = 'K'


class OroEra5(EraMeta):

    def __init__(self):
        super().__init__()
        self.dim = {0: 'ke', 1: 'je', 2: 'ie'}
        self.name = 'TOPO_CLIM'
        self.long = 'Geometric Height of the earths surface above \
            sea level (hsurf) from monthly mean \
            ERA5 climatology 1990-2019'

        self.units = 'm'


class SdEra5(EraMeta):

    def __init__(self):
        super().__init__()
        self.dim = {0: 'time', 1: 'ke', 2: 'je', 3: 'ie'}
        self.name = 'W_SNOW'
        self.long = 'Snow water equivalent of the snow-covered area \
            of a grid box (sd) from monthly mean \
            ERA5 climatology 1990-2019'

        self.units = 'kg m-2'


class SstEraI(EraMeta):

    def __init__(self):
        super().__init__()
        self.dim = {0: 'time', 1: 'ke', 2: 'je', 3: 'ie'}
        self.name = 'T_SEA'
        self.long = 'Temperature of sea water near the surface (sst) \
            from monthly mean ERA-I climatology 1986-2015'

        self.units = 'K'


class T2mEraI(EraMeta):

    def __init__(self):
        super().__init__()
        self.dim = {0: 'time', 1: 'ke', 2: 'je', 3: 'ie'}
        self.name = 'T_2M_CLIM'
        self.long = 'Temperature of air at 2m above the surface of land, \
             sea or in-land waters (2t) from monthly mean \
             ERA-I climatology 1986-2015'

        self.units = 'K'


class OroEraI(EraMeta):

    def __init__(self):
        super().__init__()
        self.dim = {0: 'ke', 1: 'je', 2: 'ie'}
        self.name = 'TOPO_CLIM'
        self.long = 'Geometric Height of the earths surface above \
            sea level (hsurf) from monthly mean \
            ERA-I climatology 1986-2015'

        self.units = 'm'


class SdEraI(EraMeta):

    def __init__(self):
        super().__init__()
        self.dim = {0: 'time', 1: 'ke', 2: 'je', 3: 'ie'}
        self.name = 'W_SNOW'
        self.long = 'Snow water equivalent of the snow-covered area \
            of a grid box (sd) from monthly mean \
            ERA-I climatology 1986-2015'

        self.units = 'kg m-2'


#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# AHF
# ->Ahf_2min
# ->Ahf_30sec


class AhfMeta:

    def __init__(self):
        self.type = np.float32
        self.standard = ''
        self.dim = {0: 'ke', 1: 'je', 2: 'ie'}
        self.name = 'AHF'
        self.units = 'W m-2'


class Ahf_2min(AhfMeta):

    def __init__(self):
        super().__init__()
        self.long = "Annual-mean anthropogenic heat flux \
            from non-renewable energy sources \
            (coal, petroleum, natural gas, and nuclear), \
            2006 after Flanner(2009) 2,5'"


class Ahf_30sec(AhfMeta):

    def __init__(self):
        super().__init__()
        self.long = 'Annual-mean anthropogenic heat flux \
            from non-renewable energy sources \
            (coal, petroleum, natural gas, and nuclear), \
            2006 after Flanner(2009) 30"'


#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# ISA
# ->Isa_30sec
# ->Isa_10sec


class IsaMeta:

    def __init__(self):
        self.type = np.float32
        self.standard = ''
        self.dim = {0: 'ke', 1: 'je', 2: 'ie'}
        self.name = 'ISA'
        self.long = 'Impervious surface area '
        self.units = '%'


class Isa_30sec(IsaMeta):

    def __init__(self):
        super().__init__()
        self.long = 'NOAA 30sec'


class Isa_10sec(IsaMeta):

    def __init__(self):
        super().__init__()
        self.long = 'European Environmental Agency 10sec'


#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# AOT


class AotMeta:

    def __init__(self):
        self.type = np.float32
        self.standard = ''
        self.dim = {0: 'time', 1: 'ntype', 2: 'ke', 3: 'je', 4: 'ie'}
        self.name = 'AOT_TG'
        self.units = ''


class AotTegen(AotMeta):

    def __init__(self):
        super().__init__()
        self.long = 'aerosol optical thickness; Tegen JGR 1997 (NASA/GISS)'


class AotAeroCom(AotMeta):

    def __init__(self):
        super().__init__()
        self.long = 'aerosol optical thickness; AeroCom1 (MPI_MET)'


# art


class ArtMeta:

    def __init__(self):
        self.type = np.float32
        self.units = '1'
        self.dim = {0: 'ke', 1: 'je', 2: 'ie'}


class ART_hcla(ArtMeta):

    def __init__(self):
        super().__init__()
        self.name = 'fr_hcla'
        self.standard = self.name + '.st'
        self.long = 'Fraction of Heavy Clay'
        self.short = self.name + '.sh'


class ART_silc(ArtMeta):

    def __init__(self):
        super().__init__()
        self.name = 'fr_silc'
        self.standard = self.name + '.st'
        self.long = 'Fraction of Silty Clay'
        self.short = self.name + '.sh'


class ART_lcla(ArtMeta):

    def __init__(self):
        super().__init__()
        self.name = 'fr_lcla'
        self.standard = self.name + '.st'
        self.long = 'Fraction of Light Clay'
        self.short = self.name + '.sh'


class ART_sicl(ArtMeta):

    def __init__(self):
        super().__init__()
        self.name = 'fr_sicl'
        self.standard = self.name + '.st'
        self.long = 'Fraction of Silty Clay Loam'
        self.short = self.name + '.sh'


class ART_cloa(ArtMeta):

    def __init__(self):
        super().__init__()
        self.name = 'fr_cloa'
        self.standard = self.name + '.st'
        self.long = 'Fraction of Clay Loam'
        self.short = self.name + '.sh'


class ART_silt(ArtMeta):

    def __init__(self):
        super().__init__()
        self.name = 'fr_silt'
        self.standard = self.name + '.st'
        self.long = 'Fraction of Silt'
        self.short = self.name + '.sh'


class ART_silo(ArtMeta):

    def __init__(self):
        super().__init__()
        self.name = 'fr_silo'
        self.standard = self.name + '.st'
        self.long = 'Fraction of Silty Loam'
        self.short = self.name + '.sh'


class ART_scla(ArtMeta):

    def __init__(self):
        super().__init__()
        self.name = 'fr_scla'
        self.standard = self.name + '.st'
        self.long = 'Fraction of Sandy Clay'
        self.short = self.name + '.sh'


class ART_loam(ArtMeta):

    def __init__(self):
        super().__init__()
        self.name = 'fr_loam'
        self.standard = self.name + '.st'
        self.long = 'Fraction of Loam'
        self.short = self.name + '.sh'


class ART_sclo(ArtMeta):

    def __init__(self):
        super().__init__()
        self.name = 'fr_sclo'
        self.standard = self.name + '.st'
        self.long = 'Fraction of Sandy Clay Loam'
        self.short = self.name + '.sh'


class ART_sloa(ArtMeta):

    def __init__(self):
        super().__init__()
        self.name = 'fr_sloa'
        self.standard = self.name + '.st'
        self.long = 'Fraction of Sandy Loam'
        self.short = self.name + '.sh'


class ART_lsan(ArtMeta):

    def __init__(self):
        super().__init__()
        self.name = 'fr_lsan'
        self.standard = self.name + '.st'
        self.long = 'Fraction of Loamy Sand'
        self.short = self.name + '.sh'


class ART_sand(ArtMeta):

    def __init__(self):
        super().__init__()
        self.name = 'fr_sand'
        self.standard = self.name + '.st'
        self.long = 'Fraction of Sand'
        self.short = self.name + '.sh'


class ART_udef(ArtMeta):

    def __init__(self):
        super().__init__()
        self.name = 'fr_udef'
        self.standard = self.name + '.st'
        self.long = 'Fraction of Undefined or Water'
        self.short = self.name + '.sh'
