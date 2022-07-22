#!/usr/bin/env python3 
import argparse
import os
import numpy as np

# extpar modules from lib
from grid_def import CosmoGrid


def main():
    parser = argparse.ArgumentParser(description='Process some integers.')

    parser.add_argument(
        '--input_cosmo_grid',
        type=str,
        required=True,
        help=
        'Fortran Namelist "INPUT_COSMO_GRID"')
    parser.add_argument(
        '--itopo_type',
        type=int,
        required=True,
        help=
        '1: GLOBE, 2: ASTER')
    parser.add_argument(
        '--lsgsl',
        action='store_true',
        help=
        'Compute subgrid-scale slope parameter (S_ORO)')
    parser.add_argument(
        '--raw_data_path',
        type=str,
        required=True,
        help=
        'Path to folder "linked_data" of exptar-input-data repository')
    args = parser.parse_args()

    args.raw_data_path = os.path.abspath(args.raw_data_path)

    tg = CosmoGrid(args.input_cosmo_grid)

    namelist = setup_namelist(tg,args)
    print(namelist)
    write_namelist(namelist)

def setup_oro_namelists(tg,args):
    namelist = {}

    # we need an extra '' -> string convenction in Fortran namelist
    raw_data_path = f"'{args.raw_data_path}'"

    # &orography_io_extpar
    namelist['orography_buffer_file'] = "'oro_buffer.nc'"
    namelist['orography_output_file'] = "'oro_cosmo.nc'"

    # &oro_runcontrol
    if args.lsgsl:
        namelist['lcompute_sgsl'] = ".TRUE."
    else:
        namelist['lcompute_sgsl'] = ".FALSE."

    # &orography_raw_data  
    namelist['itopo_type'] = args.itopo_type
    namelist['raw_data_orography_path'] = raw_data_path

    if args.itopo_type == 1:
        namelist['topo_files'] = [f"'GLOBE_{letter.upper()}10.nc' " for letter in list(map(chr,range(ord('a'),ord('p')+1)))]
        namelist['ntiles_column'] = 4
        namelist['ntiles_row'] = 4

        if args.lsgsl:
            namelist['sgsl_files'] = [f"'S_ORO_{letter.upper()}10.nc' " for letter in list(map(chr,range(ord('a'),ord('p')+1)))]

        if tg.dlon < 0.02 and tg.dlat < 0.02:
            namelist['lscale_separation'] = ".FALSE."
            namelist['lsso_param'] = ".FALSE."
        else:
            namelist['lscale_separation'] = ".TRUE."
            namelist['lsso_param'] = ".TRUE."
            
    elif args.itopo_type == 2:
        namelist.update(compute_aster_tiles(tg,args.lsgsl))
        namelist['lscale_separation'] = ".FALSE."
        namelist['lsso_param'] = ".TRUE."


    # &radtopo
    namelist['nhori'] = 24
    if tg.dlon < 0.05 and tg.dlat < 0.05:
        namelist['lradtopo'] = ".TRUE."
    else:
        namelist['lradtopo'] = ".FALSE."

    # &sgsl_raw_data
    namelist['raw_data_sgsl_path'] = raw_data_path
    namelist['idem_type'] = args.itopo_type

    return namelist

def setup_namelist(tg: CosmoGrid,args) -> dict:

    namelist = {}

    namelist.update(setup_oro_namelists(tg,args))

    return namelist

def write_namelist(namelist):
    templates_dir = '../templates'
    filled_dir = '../filled'
    names = ['INPUT_ORO', 'INPUT_RADTOPO', 'INPUT_OROSMOOTH','INPUT_SGSL']
    namelist_templates = {}

    # read namelist templates
    for name in names:
        with open(os.path.join(templates_dir,name),'r') as f:
            namelist_templates[name] = f.read() 
    
    # replace all @PLACEHOLDERS@ with real values
    for key,value in namelist.items():
        key = f'@{key.upper()}@'
        for name in names:
            if isinstance(value,list):
                namelist_templates[name] = namelist_templates[name].replace(key,str("".join(value)))
            else:
                namelist_templates[name] = namelist_templates[name].replace(key,str(value))


    # write complete namelists to file
    for name in names:
        with open(os.path.join(filled_dir,name),'w') as f:
                f.write(namelist_templates[name])
        print(name,'written')



def compute_aster_tiles(tg: CosmoGrid,lsgsl: bool) -> dict:

    zlonmax = np.amax(tg.lons)
    zlonmin = np.amin(tg.lons)
    zlatmin = np.amin(tg.lats)
    zlatmax = np.amax(tg.lats)

    # safety check
    if zlatmax > 60.0 or zlatmax < -60.0:
        raise ValueError('Domains using Aster cannot exceed 60 N or 60 S')

    aster_tiles_lon = np.empty([12,20])
    aster_tiles_lat = np.empty([12,20])

    aster_lon=-180.0
    aster_lat=60.0
    for j in range(0,20):
        for i in range(0,12):
            aster_tiles_lon[i,j] = aster_lon + float(i * 30)
            aster_tiles_lat[i,j] = aster_lat - float(j * 6)

    ilon_min = 0
    ilon_max = 0
    ilat_min = 0
    ilat_max = 0

    for j in range(0,20):
        for i in range(0,12):
            if aster_tiles_lon[i,j] < zlonmin:
                ilon_min = i
            if aster_tiles_lon[i,j] < zlonmax: 
                ilon_max = i
            if aster_tiles_lat[i,j] > zlatmin: 
                ilat_max = j
            if aster_tiles_lat[i,j] > zlatmax:
                ilat_min = j

    ntiles_column = ilon_max - ilon_min + 1
    ntiles_row = ilat_max - ilat_min + 1

    aster_files = np.empty(240,int)

    icount = 0
    for j in range(ilat_min,ilat_max+1):
        for i in range(ilon_min,ilon_max+1):
            aster_files[icount] = int(1 + i + 12*j)
            icount += 1

    namelist = {}
    namelist['ntiles_column'] = ntiles_column
    namelist['ntiles_row'] = ntiles_row
    namelist['topo_files'] = [f"'ASTER_orig_T{aster_files[idx-1]:03}.nc' " for idx in range(1,icount+1)]

    if lsgsl:
        namelist['sgsl_files'] = [f"'S_ORO_T{aster_files[idx-1]:03}.nc' " for idx in range(1,icount+1)]

    return namelist


if __name__ == '__main__':

    main()
