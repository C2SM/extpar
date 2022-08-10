#!/usr/bin/env python3 
import argparse
from WebPep import call_webpep


aerosols = {"NASA/GISS":1, "AeroCom1, MPI_MET":2}
luse = {"GLOBCOVER":1, "GLC2000":2, "ECOCLIMAP":4}
soil_type = {"FAO-DSMW":1, "HWSD TERRA":2}
alb = {"MODIS dry & sat":2, "MODIS12 vis":3} # Missing type 1
orography = {"GLOBE":1, "ASTER":2, "MERIT":3}

def extpar_client( raw_data_path, run_dir, account, host, origin_lon, origin_lat, dlon, dlat,
                   nlon, nlat, startlon, startlat, oro, landuse, soil, aot, albedo, urban,
                   orofilter=False, sgsl=False ):

    # Create namelist INPUT_COSMO_GRID
    with open("INPUT_COSMO_GRID","w") as f:
        f.write('&lmgrid \n')
        f.write('pollon = ' + str(origin_lon) + ', \n')
        f.write('pollat = ' + str(origin_lat) + ', \n')
        f.write('dlon = ' + str(dlon) + ', \n')
        f.write('dlat = ' + str(dlat) + ', \n')
        f.write('ie_tot = ' + str(nlon) + ', \n')
        f.write('je_tot = ' + str(nlat) + ', \n')
        f.write('startlon_tot = ' + str(startlon) + ', \n')
        f.write('startlat_tot = ' + str(startlat) + ', \n')
        f.write('/ \n')

    call_webpep(input_cosmo_grid="INPUT_COSMO_GRID",
            iaot_type=aerosols[aot], ilu_type=luse[landuse], ialb_type=alb[albedo], isoil_type=soil_type[soil],
            itopo_type=orography[oro],  raw_data_path=raw_data_path,
            run_dir=run_dir, account=account, host=host, lurban=urban, lsgsl=sgsl, lfilter_oro=orofilter )
