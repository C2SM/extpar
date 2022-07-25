#!/usr/bin/env python3 
import argparse
import os
import numpy as np

# extpar modules from lib
from grid_def import CosmoGrid
from utilities import launch_shell


def main():
    parser = argparse.ArgumentParser(description='Process some integers.')

    parser.add_argument(
        '--input_cosmo_grid',
        type=str,
        required=True,
        help='Fortran Namelist "INPUT_COSMO_GRID"')
    parser.add_argument(
        '--iaot_type',
        type=int,
        required=True,
        help='1: Global Aerosol Climatology Project, 2: AeroCom1')
    parser.add_argument(
        '--ilu_type',
        type=int,
        required=True,
        help='1: GLOBCOVER 2: GLC2000 4: ECOCLIMAP')
    parser.add_argument(
        '--ialb_type',
        type=int,
        required=True,
        help='1: VIS,UV,NIR 2: SOIL 3: VIS')
    parser.add_argument(
        '--isoil_type',
        type=int,
        required=True,
        help='1: FAO Digital Soil Map of the World, 2: Harmonized World Soil Database (with additional data for deepsoil): Harmonized World Soil Database (HWSD)')
    parser.add_argument(
        '--itopo_type',
        type=int,
        required=True,
        help='1: GLOBE, 2: ASTER')
    parser.add_argument(
        '--lsgsl',
        action='store_true',
        help='Compute subgrid-scale slope parameter (S_ORO)')
    parser.add_argument(
        '--lfilter_oro',
        action='store_true',
        help='Smooth orography')
    parser.add_argument(
        '--lurban',
        action='store_true',
        help='Compute parameters for urban parametrizations (AHF and ISA)')
    parser.add_argument(
        '--raw_data_path',
        type=str,
        required=True,
        help='Path to folder "linked_data" of exptar-input-data repository')
    parser.add_argument(
        '--run_dir',
        type=str,
        required=True,
        help='Folder for running Extpar')
    args = parser.parse_args()

    args.raw_data_path = os.path.abspath(args.raw_data_path)
    args.run_dir = os.path.abspath(args.run_dir)
    os.makedirs(args.run_dir,exist_ok=True)

    tg = CosmoGrid(args.input_cosmo_grid)

    namelist = setup_namelist(tg,args)
    write_namelist(args,namelist)
    create_runscript(args)

def create_runscript(args):
    files = ['submit.daint.sh']
    dir = '../templates'

    # copy modules.env
    launch_shell('cp','../modules.env', os.path.join(args.run_dir,'modules.env'))

    runscript = {}

    extpar_exes = ['"extpar_topo_to_buffer.exe" ',
                   '"extpar_landuse_to_buffer.exe" ',
                   '"extpar_aot_to_buffer.exe" ',
                   '"extpar_flake_to_buffer.exe" ',
                   '"extpar_soil_to_buffer.exe" ',
                   '"extpar_alb_to_buffer.py" ',
                   '"extpar_cru_to_buffer.py" ',
                   '"extpar_ndvi_to_buffer.py" ']

    if args.lurban:
        extpar_exes.append('"extpar_ahf_to_buffer.py" ')
        extpar_exes.append('"extpar_isa_to_buffer.py" ')
                 
    extpar_exes.append('"extpar_consistency_check.exe" ')

    runscript['extpar_executables'] = extpar_exes

    replace_placeholders(args,files,dir,runscript)


    
def setup_oro_namelist(tg,args):
    namelist = {}

    # &orography_io_extpar
    namelist['orography_buffer_file'] = 'oro_buffer.nc'
    namelist['orography_output_file'] = 'oro_cosmo.nc'

    # &oro_runcontrol
    if args.lsgsl:
        namelist['lcompute_sgsl'] = ".TRUE."
    else:
        namelist['lcompute_sgsl'] = ".FALSE."

    # &orography_raw_data  
    namelist['itopo_type'] = args.itopo_type
    namelist['raw_data_orography_path'] = args.raw_data_path

    if args.itopo_type == 1:
        namelist['topo_files'] = [f"'GLOBE_{letter.upper()}10.nc' " 
                                  for letter in 
                                  list(map(chr,range(ord('a'),ord('p') + 1)))]
        namelist['ntiles_column'] = 4
        namelist['ntiles_row'] = 4

        if args.lsgsl:
            namelist['sgsl_files'] = [f"'S_ORO_{letter.upper()}10.nc' " 
                                      for letter in 
                                      list(map(chr,range(ord('a'),ord('p') + 1)))]

        if tg.dlon < 0.02 and tg.dlat < 0.02:
            namelist['lscale_separation'] = ".FALSE."
            namelist['lsso_param'] = ".FALSE."
            namelist['scale_sep_files'] = "'placeholder_file'"
        else:
            namelist['lscale_separation'] = ".TRUE."
            namelist['scale_sep_files'] = [f"'GLOBE_{letter.upper()}_filt_lanczos_window.nc' " 
                                           for letter in 
                                           list(map(chr,range(ord('a'),ord('p') + 1)))]

            namelist['lsso_param'] = ".TRUE."

    elif args.itopo_type == 2:
        namelist.update(compute_aster_tiles(tg,args.lsgsl))
        namelist['lscale_separation'] = ".FALSE."
        namelist['scale_sep_files'] = "'placeholder_file'"
        namelist['lsso_param'] = ".TRUE."

    # &scale_separated_raw_data
    # other paramters of the namelist already set earlier
    namelist['raw_data_scale_sep_path'] = args.raw_data_path

    # &orography_smoothing
    if args.lfilter_oro:
        namelist['lfilter_oro'] = ".TRUE."
    else:
        namelist['lfilter_oro'] = ".FALSE."

    namelist['ilow_pass_oro'] = 4
    namelist['numfilt_oro'] = 1
    namelist['ilow_pass_xso'] = 5
    namelist['lxso_first'] = ".FALSE."
    namelist['numfilt_xso'] = 1
    namelist['rxso_mask'] = 750.0
    namelist['eps_filter'] = 0.1
    namelist['rfill_valley'] = 0.0
    namelist['ifill_valley'] = 1

    # &radtopo
    namelist['nhori'] = 24
    if tg.dlon < 0.05 and tg.dlat < 0.05:
        namelist['lradtopo'] = ".TRUE."
    else:
        namelist['lradtopo'] = ".FALSE."

    # &sgsl_raw_data
    namelist['raw_data_sgsl_path'] = args.raw_data_path
    namelist['idem_type'] = args.itopo_type

    return namelist


def setup_lu_namelist(args):
    namelist = {}
    namelist['i_landuse_data'] = args.ilu_type
    namelist['ilookup_table_lu'] = args.ilu_type
    namelist['raw_data_lu_path'] = args.raw_data_path
    namelist['raw_data_glcc_path'] = args.raw_data_path
    namelist['lu_buffer_file'] = 'lu_buffer.nc'
    namelist['raw_data_glcc_filename'] = 'GLCC_usgs_class_byte'
    namelist['glcc_buffer_file'] = 'glcc_buffer.nc'
    namelist['l_use_corine'] = ".FALSE."
    if args.ilu_type == 1:
        namelist['raw_data_lu_filename'] = [f"'GLOBCOVER_{i}_16bit.nc' " 
                                            for i in range(0,6)]
    elif args.ilu_type == 2:
        # we need "" padding for correct replacement in Fortran namelist
        namelist['raw_data_lu_filename'] = "'GLC2000_byte.nc'"
    elif args.ilu_type == 4:
        # we need "" padding for correct replacement in Fortran namelist
        namelist['raw_data_lu_filename'] = "'GLCC_usgs_class_byte.nc'"

    return namelist


def setup_aot_namelist(args):
    namelist = {}
    namelist['iaot_type'] = args.iaot_type
    namelist['raw_data_aot_path'] = args.raw_data_path
    namelist['aot_buffer_file'] = 'aot_buffer.nc'
    if args.iaot_type == 1:
        namelist['raw_data_aot_filename'] = 'aot_GACP.nc'
    elif args.iaot_type == 2:
        namelist['raw_data_aot_filename'] = 'aod_AeroCom1.nc'

    return namelist


def setup_soil_namelist(args):
    namelist = {}

    # &soil_raw_data 
    namelist['isoil_data'] = args.isoil_type
    namelist['raw_data_soil_path'] = args.raw_data_path

    if args.isoil_type == 1:
        namelist['raw_data_soil_filename'] = 'FAO_DSMW_double.nc'
        namelist['ldeep_soil'] = ".FALSE"
    elif args.isoil_type == 2:
        namelist['raw_data_soil_filename'] = 'HWSD0_30_topsoil.nc'
        namelist['raw_data_deep_soil_filename'] = 'HWSD30_100_subsoil.nc'
        namelist['ldeep_soil'] = ".TRUE"
    elif args.isoil_type == 3:
        namelist['raw_data_soil_filename'] = 'HWSD0_30_topsoil.nc'
        namelist['raw_data_deep_soil_filename'] = 'HWSD30_100_subsoil.nc'
        namelist['ldeep_soil'] = ".FALSE"

    # &soil_io_extpar
    namelist['soil_buffer_file'] = 'soil_buffer.nc'

    # &HWSD_index_files

    # Quickfix, should be in linked_data as well
    namelist['path_HWSD_index_files'] = os.path.join(args.raw_data_path,'../soil')
    namelist['lookup_table_HWSD'] = 'LU_TAB_HWSD_UF.data'
    namelist['HWSD_data'] = 'HWSD_DATA_COSMO.data'
    namelist['HWSD_data_deep'] = 'HWSD_DATA_COSMO_S.data'

    return namelist


def setup_tclim_namelist(args):
    namelist = {}

    namelist['it_cl_type'] = 1
    namelist['raw_data_t_clim_path'] = args.raw_data_path
    namelist['raw_data_tclim_coarse'] = 'absolute_hadcrut3.nc'
    namelist['raw_data_tclim_fine'] = 'CRU_T_SOIL_clim.nc'
    namelist['t_clim_buffer_file'] = 'tclim_buffer.nc'

    return namelist


def setup_flake_namelist(args):
    namelist = {}

    namelist['raw_data_flake_path'] = args.raw_data_path
    namelist['raw_data_flake_filename'] = 'GLDB_lakedepth.nc'
    namelist['flake_buffer_file'] = 'flake_buffer.nc'

    return namelist


def setup_albedo_namelist(args):
    namelist = {}

    namelist['raw_data_alb_path'] = args.raw_data_path
    namelist['ialb_type'] = args.ialb_type
    namelist['alb_buffer_file'] = 'alb_buffer.nc'

    if args.ialb_type == 1:
        namelist['raw_data_alb_filename'] = 'alb_new.nc'
        namelist['raw_data_alnid_filename'] = 'alnid_new.nc'
        namelist['raw_data_aluvd_filename'] = 'aluvd_new.nc'
    elif args.ialb_type == 2:
        namelist['raw_data_alb_filename'] = 'global_soil_albedo.nc'
    elif args.ialb_type == 3:
        namelist['raw_data_alb_filename'] = 'alb_new.nc'

    return namelist


def setup_ndvi_namelist(args):
    namelist = {}

    namelist['raw_data_ndvi_path'] = args.raw_data_path
    namelist['raw_data_ndvi_filename'] = 'NDVI_1998_2003.nc'
    namelist['ndvi_buffer_file'] = 'ndvi_buffer.nc'

    return namelist


def setup_urban_namelist(args):
    namelist = {}

    # input_ahf
    namelist['iahf_type'] = 1
    namelist['raw_data_ahf_path'] = args.raw_data_path
    namelist['raw_data_ahf_filename'] = 'AHF_2006_2.5min_lonlat.nc'
    namelist['ahf_buffer_file'] = 'ahf_buffer.nc'

    # input_isa
    namelist['isa_type'] = 1
    namelist['raw_data_isa_path'] = args.raw_data_path
    namelist['raw_data_isa_filename'] = 'NOAA_ISA_16bit_lonlat.nc'
    namelist['isa_buffer_file'] = 'isa_buffer.nc'

    return namelist


def setup_check_namelist(args):
    namelist = {}

    namelist['netcdf_output_filename'] = 'external_parameter.nc'
    namelist['i_lsm_data'] = 1
    namelist['land_sea_mask_file'] = ""
    namelist['number_special_points'] = 0
    namelist['lflake_correction'] = ".TRUE."

    return namelist


def setup_namelist(tg: CosmoGrid,args) -> dict:

    namelist = {}

    namelist.update(setup_oro_namelist(tg,args))
    namelist.update(setup_albedo_namelist(args))
    namelist.update(setup_aot_namelist(args))
    namelist.update(setup_tclim_namelist(args))
    namelist.update(setup_lu_namelist(args))
    namelist.update(setup_flake_namelist(args))
    namelist.update(setup_ndvi_namelist(args))
    namelist.update(setup_urban_namelist(args))
    namelist.update(setup_soil_namelist(args))
    namelist.update(setup_check_namelist(args))

    return namelist


def write_namelist(args,namelist):
    templates_dir = '../templates'
    files = ['INPUT_ORO',
                 'INPUT_RADTOPO',
                 'INPUT_OROSMOOTH',
                 'INPUT_SGSL',
                 'INPUT_AOT',
                 'INPUT_LU',
                 'INPUT_FLAKE',
                 'INPUT_SCALE_SEP',
                 'INPUT_SOIL',
                 'INPUT_CHECK',
                 'namelist.py']

    replace_placeholders(args,files,templates_dir,namelist)


def replace_placeholders(args,templates,dir,actual_values):
    all_templates = {}

    # read templates
    for template in templates:
        with open(os.path.join(dir,template),'r') as f:
            all_templates[template] = f.read() 

    # replace all @PLACEHOLDERS@ with real values
    for key,value in actual_values.items():
        key = f'@{key.upper()}@'
        for template in templates:
            if isinstance(value,list):
                all_templates[template] = all_templates[template].replace(key,str("".join(value)))
            else:
                all_templates[template] = all_templates[template].replace(key,str(value))

    # write complete template to file
    for template in templates:
        with open(os.path.join(args.run_dir,template),'w') as f:
            f.write(all_templates[template])
        print(template,'written')

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

    aster_lon = -180.0
    aster_lat = 60.0
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
    for j in range(ilat_min,ilat_max + 1):
        for i in range(ilon_min,ilon_max + 1):
            aster_files[icount] = int(1 + i + 12 * j)
            icount += 1

    namelist = {}
    namelist['ntiles_column'] = ntiles_column
    namelist['ntiles_row'] = ntiles_row
    namelist['topo_files'] = [f"'ASTER_orig_T{aster_files[idx-1]:03}.nc' " 
                              for idx in range(1,icount + 1)]

    if lsgsl:
        namelist['sgsl_files'] = [f"'S_ORO_T{aster_files[idx-1]:03}.nc' " 
                                  for idx in range(1,icount + 1)]

    return namelist


if __name__ == '__main__':

    main()
