#!/usr/bin/env python3
import argparse
import os
import stat
import shutil
import logging
import json

import numpy as np

# extpar modules from lib
try:
    from extpar.lib.grid_def import CosmoGrid, IconGrid
    from extpar.lib.utilities import launch_shell
    DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
except ImportError:  # package not installed -> use PYTHONPATH
    from grid_def import CosmoGrid, IconGrid
    from utilities import launch_shell
    DATA_DIR = os.path.join(os.path.dirname(__file__), '..')


def main():

    # initialize logger
    logging.basicConfig(level=logging.INFO, format='%(message)s')

    parser = argparse.ArgumentParser(description='Process some integers.')

    parser.add_argument('--extpar-config',
                        type=str,
                        required=True,
                        help='Path to extpar config file')

    parser.add_argument(
        '--input-grid',
        type=str,
        help='COSMO: Fortran Namelist "INPUT_COSMO_GRID", ICON: Icon grid file'
    )

    parser.add_argument(
        '--raw-data-path',
        type=str,
        required=True,
        help='Path to folder "linked_data" of exptar-input-data repository')
    parser.add_argument('--run-dir',
                        type=str,
                        required=True,
                        help='Folder for running Extpar')
    parser.add_argument('--account',
                        type=str,
                        required=True,
                        help='Account for slurm job')
    parser.add_argument('--host', type=str, required=True, help='Host')
    parser.add_argument('--no-batch-job',
                        action='store_true',
                        help="Run jobscript not as batch job")

    args = parser.parse_args()

    # Read and parse the JSON configuration file
    with open(args.extpar_config, 'r') as f:
        config = json.load(f)

    config = config['extpar']

    # main use case is icon -> igrid_type=1
    igrid_type = config.get('igrid_type', 1)
    iaot_type = config.get('iaot_type')
    ilu_type = config.get('ilu_type')
    ialb_type = config.get('ialb_type')
    isoil_type = config.get('isoil_type')
    itopo_type = config.get('itopo_type')
    it_cl_type = config.get('it_cl_type')
    iera_type = config.get('iera_type')
    iemiss_type = config.get('iemiss_type')
    ilookup_table_lu = config.get('ilookup_table_lu')
    enable_cdnc = config.get('enable_cdnc', False)
    enable_edgar = config.get('enable_edgar', False)
    enable_art = config.get('enable_art', False)
    use_array_cache = config.get('use_array_cache', False)
    lsgsl = config.get('lsgsl', False)
    lfilter_oro = config.get('lfilter_oro', False)
    lurban = config.get('lurban', False)
    l_use_corine = config.get('l_use_corine', False)
    lradtopo = config.get('lradtopo', False)
    nhori = config.get('nhori', 24)
    radtopo_radius = config.get('radtopo_radius', 40000.0)
    tcorr_lapse_rate = config.get('tcorr_lapse_rate', 0.0065)
    tcorr_offset = config.get('tcorr_offset', 0.0)

    generate_external_parameters(
        igrid_type, args.input_grid, iaot_type, ilu_type, ialb_type,
        isoil_type, itopo_type, it_cl_type, iera_type, iemiss_type,
        ilookup_table_lu, enable_cdnc, enable_edgar, enable_art,
        use_array_cache, nhori, radtopo_radius, tcorr_lapse_rate, tcorr_offset,
        args.raw_data_path, args.run_dir, args.account, args.host,
        args.no_batch_job, lurban, lsgsl, lfilter_oro, l_use_corine, lradtopo)


def generate_external_parameters(igrid_type,
                                 input_grid,
                                 iaot_type,
                                 ilu_type,
                                 ialb_type,
                                 isoil_type,
                                 itopo_type,
                                 it_cl_type,
                                 iera_type,
                                 iemiss_type,
                                 ilookup_table_lu,
                                 enable_cdnc,
                                 enable_edgar,
                                 enable_art,
                                 use_array_cache,
                                 nhori,
                                 radtopo_radius,
                                 tcorr_lapse_rate,
                                 tcorr_offset,
                                 raw_data_path,
                                 run_dir,
                                 account,
                                 host,
                                 no_batch_job=False,
                                 lurban=False,
                                 lsgsl=False,
                                 lfilter_oro=False,
                                 l_use_corine=False,
                                 lradtopo=False):

    # initialize logger
    logging.basicConfig(level=logging.INFO, format='%(message)s')

    # make all paths absolut
    raw_data_path = os.path.abspath(raw_data_path)
    run_dir = os.path.abspath(run_dir)
    input_grid = os.path.abspath(input_grid)

    args_dict = {
        'input_grid': input_grid,
        'igrid_type': igrid_type,
        'iaot_type': iaot_type,
        'ilu_type': ilu_type,
        'ialb_type': ialb_type,
        'isoil_type': isoil_type,
        'itopo_type': itopo_type,
        'it_cl_type': it_cl_type,
        'iera_type': iera_type,
        'iemiss_type': iemiss_type,
        'ilookup_table_lu': ilookup_table_lu,
        'enable_cdnc': enable_cdnc,
        'enable_edgar': enable_edgar,
        'enable_art': enable_art,
        'use_array_cache': use_array_cache,
        'l_use_corine': l_use_corine,
        'lradtopo': lradtopo,
        'nhori': nhori,
        'radtopo_radius': radtopo_radius,
        'tcorr_lapse_rate': tcorr_lapse_rate,
        'tcorr_offset': tcorr_offset,
        'lsgsl': lsgsl,
        'lfilter_oro': lfilter_oro,
        'lurban': lurban,
        'raw_data_path': raw_data_path,
        'run_dir': run_dir,
        'account': account,
        'host': host,
        'no_batch_job': no_batch_job
    }

    namelist = setup_namelist(args_dict)

    runscript = setup_runscript(args_dict)

    prepare_sandbox(args_dict, namelist, runscript)

    run_extpar(args_dict)


def run_extpar(args):

    exec_file = os.path.join(args['run_dir'], f'submit.{args["host"]}.sh')
    os.chmod(exec_file, 0o755)

    if args['no_batch_job']:
        logging.info("run jobscript")
        cwd = os.getcwd()
        os.chdir(args['run_dir'])
        logging.info("run job script")
        launch_shell(exec_file)
        os.chdir(cwd)

    else:
        logging.info("submit job and wait for it's completion")
        launch_shell('sbatch', '--chdir', args['run_dir'], '--wait', exec_file)
    logging.info("job finished")


def prepare_sandbox(args, namelist, runscript, test_run=False):

    os.makedirs(args['run_dir'], exist_ok=True)
    write_namelist(args, namelist)
    write_runscript(args, runscript)
    if not test_run:
        copy_required_files(args, runscript['extpar_executables'])


def write_runscript(args, runscript):
    dir = os.path.join(DATA_DIR, "templates")
    files = [f'submit.{args["host"]}.sh']

    replace_placeholders(args, files, dir, runscript)


def write_namelist(args, namelist):
    templates_dir = os.path.join(DATA_DIR, "templates")
    files = [
        'INPUT_ORO', 'INPUT_RADTOPO', 'INPUT_OROSMOOTH', 'INPUT_SGSL',
        'INPUT_LU', 'INPUT_FLAKE', 'INPUT_SCALE_SEP', 'INPUT_SOIL',
        'INPUT_CHECK', 'namelist'
    ]

    replace_placeholders(args, files, templates_dir, namelist)

    igrid_type = args['igrid_type']
    # INPUT_grid_org
    with open(os.path.join(args['run_dir'], 'INPUT_grid_org'), 'w') as f:
        f.write('&GRID_DEF \n')
        f.write(f'igrid_type = {igrid_type} \n')
        f.write("domain_def_namelist = 'INPUT_GRID' \n")
        f.write('/ \n')

    # COSMO_GRID
    if igrid_type == 2:
        shutil.copy(args['input_grid'],
                    os.path.join(args['run_dir'], 'INPUT_GRID'))
    # ICON_GRID
    else:
        grid_path = os.path.dirname(args['input_grid'])
        grid_file = os.path.basename(args['input_grid'])
        with open(os.path.join(args['run_dir'], 'INPUT_GRID'), 'w') as f:
            f.write('&icon_grid_info \n')
            f.write(f'icon_grid_dir = "{grid_path}" \n')
            f.write(f"icon_grid_nc_file = '{grid_file}' \n")
            f.write('/ \n')


def copy_required_files(args, executables):

    shutil.copy(os.path.join(DATA_DIR, 'modules.env'),
                os.path.join(args['run_dir'], 'modules.env'))

    shutil.copy(
        os.path.join(DATA_DIR, 'test/testsuite/bin/runcontrol_functions.sh'),
        os.path.join(args['run_dir'], 'runcontrol_functions.sh'))

    for exe in executables:
        exe = exe.replace('"', '')
        exe = exe.strip()
        target = os.path.join(args['run_dir'], exe)
        shutil.copy(os.path.join(DATA_DIR, 'bin', exe), target)
        st = os.stat(target)
        os.chmod(target, st.st_mode | stat.S_IEXEC)


def setup_oro_namelist(args):

    igrid_type = args['igrid_type']

    if igrid_type == 2:
        return setup_oro_namelist_cosmo(args)
    else:
        tg = IconGrid(args['input_grid'])
        lonmax = np.amax(tg.lons)
        lonmin = np.amin(tg.lons)
        latmin = np.amin(tg.lats)
        latmax = np.amax(tg.lats)
        return setup_oro_namelist_icon(args, lonmax, lonmin, latmax, latmin)


def setup_oro_namelist_cosmo(args):

    tg = CosmoGrid(args['input_grid'])

    namelist = {}

    # &orography_io_extpar
    namelist['orography_buffer_file'] = 'oro_buffer.nc'
    namelist['orography_output_file'] = 'oro_grid.nc'

    # &oro_runcontrol
    if args['lsgsl']:
        namelist['lcompute_sgsl'] = ".TRUE."
        namelist['sgsl_buffer_file'] = 'sgsl_buffer.nc'
    else:
        namelist['lcompute_sgsl'] = ".FALSE."
        namelist['sgsl_buffer_file'] = 'placeholder_file'

    # &orography_raw_data
    namelist['itopo_type'] = args['itopo_type']
    namelist['raw_data_orography_path'] = args['raw_data_path']

    if args['lradtopo']:
        tg_ext = extend_cosmo_grid_for_radtopo(args["run_dir"], tg)
        tg_for_extent = tg_ext
    else:
        tg_for_extent = tg

    lonmax = np.amax(tg_for_extent.lons)
    lonmin = np.amin(tg_for_extent.lons)
    latmin = np.amin(tg_for_extent.lats)
    latmax = np.amax(tg_for_extent.lats)

    if args['itopo_type'] == 1:
        namelist['topo_files'] = generate_globe_filenames()
        namelist['ntiles_column'] = 4
        namelist['ntiles_row'] = 4

        if args['lsgsl']:
            namelist['sgsl_files'] = [
                f"'S_ORO_{letter.upper()}10.nc' "
                for letter in list(map(chr, range(ord('a'),
                                                  ord('p') + 1)))
            ]
            namelist['lpreproc_oro'] = ".FALSE."
        else:
            namelist['lpreproc_oro'] = ".FALSE."
            namelist['sgsl_files'] = 'placeholder_file'

        if tg.dlon < 0.02 and tg.dlat < 0.02:
            namelist['lscale_separation'] = ".FALSE."
            namelist['lsso_param'] = ".FALSE."
            namelist['scale_sep_files'] = "'placeholder_file'"
        else:
            namelist['lscale_separation'] = ".TRUE."
            namelist['scale_sep_files'] = [
                f"'GLOBE_{letter.upper()}_filt_lanczos_window.nc' "
                for letter in list(map(chr, range(ord('a'),
                                                  ord('p') + 1)))
            ]

            namelist['lsso_param'] = ".TRUE."

    elif args['itopo_type'] == 2:
        namelist.update(
            compute_aster_tiles(lonmax=lonmax,
                                lonmin=lonmin,
                                latmax=latmax,
                                latmin=latmin,
                                lsgsl=args['lsgsl']))
        namelist['lscale_separation'] = ".FALSE."
        namelist['scale_sep_files'] = "'placeholder_file'"
        namelist['lsso_param'] = ".TRUE."
    elif args['itopo_type'] == 3:
        namelist.update(
            compute_merit_tiles(lonmax=lonmax,
                                lonmin=lonmin,
                                latmax=latmax,
                                latmin=latmin,
                                lsgsl=args['lsgsl']))
        namelist['lscale_separation'] = ".FALSE."
        namelist['scale_sep_files'] = "'placeholder_file'"
        namelist['lsso_param'] = ".TRUE."

    else:
        logging.error(f'Unknown itopo_type {args["itopo_type"]}')
        raise ValueError(f'Unknown itopo_type {args["itopo_type"]}')

    # &scale_separated_raw_data
    # other paramters of the namelist already set earlier
    namelist['raw_data_scale_sep_path'] = args['raw_data_path']

    # &orography_smoothing
    if args['lfilter_oro']:
        namelist['lfilter_oro'] = ".TRUE."
    else:
        namelist['lfilter_oro'] = ".FALSE."

    namelist.update(orography_smoothing_params())

    # &radtopo
    namelist['nhori'] = args['nhori']
    if args['lradtopo']:
        namelist['lradtopo'] = ".TRUE."
    else:
        namelist['lradtopo'] = ".FALSE."

    # not relevant for COSMO grid, but required for namelist
    namelist['max_missing'] = 0.95
    namelist['min_circ_cov'] = 1
    namelist['radius'] = args['radtopo_radius']
    namelist['itype_scaling'] = 0

    # &sgsl_raw_data
    namelist['raw_data_sgsl_path'] = args['raw_data_path']
    namelist['idem_type'] = args['itopo_type']

    return namelist


def orography_smoothing_params():
    namelist = {}
    namelist['ilow_pass_oro'] = 4
    namelist['numfilt_oro'] = 1
    namelist['ilow_pass_xso'] = 5
    namelist['lxso_first'] = ".FALSE."
    namelist['numfilt_xso'] = 1
    namelist['rxso_mask'] = 750.0
    namelist['eps_filter'] = 0.1
    namelist['rfill_valley'] = 0.0
    namelist['ifill_valley'] = 1
    return namelist


def setup_oro_namelist_icon(args, lonmax, lonmin, latmax, latmin):

    namelist = {}

    # &orography_io_extpar
    namelist['orography_buffer_file'] = 'oro_buffer.nc'
    namelist['orography_output_file'] = 'oro_grid.nc'

    namelist['lcompute_sgsl'] = ".FALSE."
    namelist['sgsl_buffer_file'] = 'placeholder_file'

    # &orography_raw_data
    namelist['itopo_type'] = args['itopo_type']
    namelist['raw_data_orography_path'] = args['raw_data_path']

    if args['itopo_type'] == 1:
        namelist['topo_files'] = generate_globe_filenames()
        namelist['ntiles_column'] = 4
        namelist['ntiles_row'] = 4

        namelist['lscale_separation'] = ".FALSE."
        namelist['scale_sep_files'] = "'placeholder_file'"
        namelist['lsso_param'] = ".TRUE."

    elif args['itopo_type'] == 2:
        namelist.update(
            compute_aster_tiles(lonmax=lonmax,
                                lonmin=lonmin,
                                latmax=latmax,
                                latmin=latmin,
                                lsgsl=False))
        namelist['lscale_separation'] = ".FALSE."
        namelist['scale_sep_files'] = "'placeholder_file'"
        namelist['lsso_param'] = ".TRUE."
    elif args['itopo_type'] == 3:
        namelist.update(
            compute_merit_tiles(lonmax=lonmax,
                                lonmin=lonmin,
                                latmax=latmax,
                                latmin=latmin,
                                lsgsl=False))
        namelist['lscale_separation'] = ".FALSE."
        namelist['scale_sep_files'] = "'placeholder_file'"
        namelist['lsso_param'] = ".TRUE."

    else:
        logging.error(f'Unknown itopo_type {args["itopo_type"]}')
        raise ValueError(f'Unknown itopo_type {args["itopo_type"]}')

    # &scale_separated_raw_data
    # other paramters of the namelist already set earlier
    namelist['raw_data_scale_sep_path'] = args['raw_data_path']

    # &orography_smoothing
    namelist['lfilter_oro'] = ".FALSE."

    # not relevant for ICON grid, but required for namelist
    namelist.update(orography_smoothing_params())
    namelist['lpreproc_oro'] = ".FALSE."
    namelist['sgsl_files'] = 'placeholder_file'
    namelist['idem_type'] = args['itopo_type']
    namelist['raw_data_sgsl_path'] = args['raw_data_path']

    # &radtopo
    if args['lradtopo']:
        namelist['lradtopo'] = ".TRUE."
    else:
        namelist['lradtopo'] = ".FALSE."

    # only relevant if lradtopo=.TRUE., but needed for namelist
    namelist['nhori'] = args['nhori']
    namelist['max_missing'] = 0.95
    namelist['min_circ_cov'] = 1
    namelist['radius'] = args['radtopo_radius']
    namelist['itype_scaling'] = 0

    return namelist


def generate_globe_filenames():
    return [
        f"'GLOBE_{letter.upper()}10.nc' "
        for letter in list(map(chr, range(ord('a'),
                                          ord('p') + 1)))
    ]


def setup_lu_namelist(args):
    namelist = {}
    namelist['i_landuse_data'] = args['ilu_type']
    namelist['ilookup_table_lu'] = args['ilookup_table_lu']
    namelist['raw_data_lu_path'] = args['raw_data_path']
    namelist['raw_data_glcc_path'] = args['raw_data_path']
    namelist['lu_buffer_file'] = 'lu_buffer.nc'
    namelist['raw_data_glcc_filename'] = 'GLCC_usgs_class_byte.nc'
    namelist['glcc_buffer_file'] = 'glcc_buffer.nc'
    namelist['ntiles_globcover'] = 6
    namelist['l_terra_urb'] = ".FALSE."

    if args['l_use_corine']:
        namelist['l_use_corine'] = ".TRUE."
    else:
        namelist['l_use_corine'] = ".FALSE."

    if args['ilu_type'] == 1:
        if args['l_use_corine']:
            namelist['raw_data_lu_filename'] = "'CORINE_globcover.nc'"
            namelist['ntiles_globcover'] = 1
        else:
            namelist['raw_data_lu_filename'] = [
                f"'GLOBCOVER_{i}_16bit.nc' " for i in range(0, 6)
            ]
    elif args['ilu_type'] == 2:
        # we need "" padding for correct replacement in Fortran namelist
        namelist['raw_data_lu_filename'] = "'GLC2000_byte.nc'"
    elif args['ilu_type'] == 6:
        # we need "" padding for correct replacement in Fortran namelist
        namelist['raw_data_lu_filename'] = "'ECOCLIMAP_SG.nc'"
        namelist['l_terra_urb'] = ".TRUE."
    else:
        logging.error(f'Unknown ilu_type {args["ilu_type"]}')
        raise ValueError(f'Unknown ilu_type {args["ilu_type"]}')

    return namelist


def setup_aot_namelist(args):
    namelist = {}
    namelist['iaot_type'] = args['iaot_type']
    namelist['raw_data_aot_path'] = args['raw_data_path']
    namelist['aot_buffer_file'] = 'aot_buffer.nc'
    if args['iaot_type'] == 1:
        namelist['raw_data_aot_filename'] = 'aot_GACP.nc'
    elif args['iaot_type'] == 2:
        namelist['raw_data_aot_filename'] = 'aod_AeroCom1.nc'
    else:
        logging.error(f'Unknown iaot_type {args["iaot_type"]}')
        raise ValueError(f'Unknown iaot_type {args["iaot_type"]}')

    return namelist


def setup_soil_namelist(args):
    namelist = {}

    # &soil_raw_data
    namelist['isoil_data'] = args['isoil_type']
    namelist['raw_data_soil_path'] = args['raw_data_path']

    if args['isoil_type'] == 1:
        namelist['raw_data_soil_filename'] = 'FAO_DSMW_double.nc'
    elif args['isoil_type'] == 2:
        namelist['raw_data_soil_filename'] = 'HWSD0_30_topsoil.nc'
    elif args['isoil_type'] == 3:
        namelist['raw_data_soil_filename'] = 'HWSD0_30_topsoil.nc'

    else:
        logging.error(f'Unknown isoil_type {args["isoil_type"]}')
        raise ValueError(f'Unknown isoil_type {args["isoil_type"]}')

    # &soil_io_extpar
    namelist['soil_buffer_file'] = 'soil_buffer.nc'

    # &HWSD_index_files

    # Quickfix, should be in linked_data as well
    namelist['path_HWSD_index_files'] = os.path.join(args['raw_data_path'],
                                                     '../soil')
    namelist['lookup_table_HWSD'] = 'LU_TAB_HWSD_UF.data'
    namelist['HWSD_data'] = 'HWSD_DATA_COSMO.data'

    return namelist


def setup_tclim_namelist(args):
    namelist = {}

    namelist['raw_data_t_clim_path'] = args['raw_data_path']
    namelist['t_clim_buffer_file'] = 'tclim_buffer.nc'

    namelist['it_cl_type'] = args['it_cl_type']
    if args['it_cl_type'] > 2:
        raise ValueError(f'Unknown it_cl_type {args["it_cl_type"]}')

    namelist['raw_data_tclim_coarse'] = 'absolute_hadcrut3.nc'
    namelist['raw_data_tclim_fine'] = 'CRU_T_SOIL_clim.nc'

    namelist['tcorr_lapse_rate'] = args['tcorr_lapse_rate']
    namelist['tcorr_offset'] = args['tcorr_offset']

    return namelist


def setup_flake_namelist(args):
    namelist = {}

    namelist['raw_data_flake_path'] = args['raw_data_path']
    namelist['raw_data_flake_filename'] = 'GLDB_lakedepth.nc'
    namelist['flake_buffer_file'] = 'flake_buffer.nc'

    return namelist


def setup_art_namelist(args):
    namelist = {}

    namelist['raw_data_art_path'] = args['raw_data_path']
    namelist['raw_data_art_filename'] = 'HWSD0_USDA.nc'
    namelist['art_buffer_file'] = 'art_buffer.nc'

    return namelist


def setup_albedo_namelist(args):
    namelist = {}

    namelist['raw_data_alb_path'] = args['raw_data_path']
    namelist['ialb_type'] = args['ialb_type']
    namelist['alb_buffer_file'] = 'alb_buffer.nc'

    if args['ialb_type'] == 1:
        namelist['raw_data_alb_filename'] = 'alb_new.nc'
        namelist['raw_data_alnid_filename'] = 'alnid_new.nc'
        namelist['raw_data_aluvd_filename'] = 'aluvd_new.nc'
    elif args['ialb_type'] == 2:
        namelist['raw_data_alb_filename'] = 'global_soil_albedo.nc'
    elif args['ialb_type'] == 3:
        namelist['raw_data_alb_filename'] = 'alb_new.nc'
    else:
        logging.error(f'Unknown ialb_type {args["ialb_type"]}')
        raise ValueError(f'Unknown ialb_type {args["ialb_type"]}')

    return namelist


def setup_ndvi_namelist(args):
    namelist = {}

    namelist['raw_data_ndvi_path'] = args['raw_data_path']
    namelist['raw_data_ndvi_filename'] = 'NDVI_1998_2003.nc'
    namelist['ndvi_buffer_file'] = 'ndvi_buffer.nc'

    return namelist


def setup_era_namelist(args):
    namelist = {}
    iera_type = args['iera_type']

    namelist['iera_type'] = iera_type
    namelist['era_buffer_file'] = 'era_buffer.nc'
    namelist['raw_data_era_path'] = args['raw_data_path']

    if iera_type == 1:
        namelist['raw_data_era_ORO'] = 'ERA5_ORO_1990.nc'
        namelist['raw_data_era_SD'] = 'ERA5_SD_1990_2019.nc'
        namelist['raw_data_era_T2M'] = 'ERA5_T2M_1990_2019.nc'
        namelist['raw_data_era_SST'] = 'ERA5_SST_1990_2019.nc'
    elif iera_type == 2:
        namelist['raw_data_era_ORO'] = 'ERA-I_ORO_1986.nc'
        namelist['raw_data_era_SD'] = 'ERA-I_SD_1986_2015.nc'
        namelist['raw_data_era_T2M'] = 'ERA-I_T2M_1986_2015.nc'
        namelist['raw_data_era_SST'] = 'ERA-I_SST_1986_2015.nc'
    else:
        raise ValueError(f'Unknown iera_type {iera_type}')

    return namelist


def setup_emiss_namelist(args):
    namelist = {}
    iemiss_type = args['iemiss_type']

    namelist['iemiss_type'] = iemiss_type
    namelist['emiss_buffer_file'] = 'emiss_buffer.nc'
    namelist['raw_data_emiss_path'] = args['raw_data_path']
    if iemiss_type == 1:
        namelist['raw_data_emiss_filename'] = 'CAMEL_bbe_full_2010-2015.nc'
    elif iemiss_type == 2:
        namelist['raw_data_emiss_filename'] = 'CAMEL_bbe_lw_2010-2015.nc'
    else:
        raise ValueError(f'Unknown iemiss_type {iemiss_type}')

    return namelist


def setup_urban_namelist(args):
    namelist = {}

    # input_ahf
    namelist['iahf_type'] = 1
    namelist['raw_data_ahf_path'] = args['raw_data_path']
    namelist['raw_data_ahf_filename'] = 'AHF_2006_2.5min_lonlat.nc'
    namelist['ahf_buffer_file'] = 'ahf_buffer.nc'

    # input_isa
    namelist['isa_type'] = 1
    namelist['raw_data_isa_path'] = args['raw_data_path']
    namelist['raw_data_isa_filename'] = 'NOAA_ISA_16bit_lonlat.nc'
    namelist['isa_buffer_file'] = 'isa_buffer.nc'

    return namelist


def setup_cdnc_namelist(args):
    namelist = {}
    icdnc_type = args['icdnc_type']

    namelist['raw_data_cdnc_path'] = args['raw_data_path']
    namelist['cdnc_buffer_file'] = 'cdnc_buffer.nc'
    if icdnc_type == 1:
        namelist['raw_data_cdnc_filename'] = 'modis_cdnc_climatology_Q06.nc'
    elif icdnc_type == 2:
        namelist[
            'raw_data_cdnc_filename'] = 'modis_cdnc_climatology_BR17_37mu_adjusted.nc'

    return namelist


def setup_edgar_namelist(args):
    namelist = {}

    namelist['raw_data_edgar_path'] = args['raw_data_path']
    namelist[
        'raw_data_edgar_filename_bc'] = 'v8.1_FT2022_AP_BC_2022_TOTALS_flx.nc'
    namelist[
        'raw_data_edgar_filename_oc'] = 'v8.1_FT2022_AP_OC_2022_TOTALS_flx.nc'
    namelist[
        'raw_data_edgar_filename_so2'] = 'v8.1_FT2022_AP_SO2_2022_TOTALS_flx.nc'
    namelist[
        'raw_data_edgar_filename_nox'] = 'v8.1_FT2022_AP_NOx_2022_TOTALS_flx.nc'
    namelist[
        'raw_data_edgar_filename_nh3'] = 'v8.1_FT2022_AP_NH3_2022_TOTALS_flx.nc'
    namelist['edgar_buffer_file'] = 'edgar_buffer.nc'

    return namelist


def setup_check_namelist(args):
    namelist = {}

    namelist['netcdf_output_filename'] = 'external_parameter.nc'
    namelist['i_lsm_data'] = 1
    namelist['land_sea_mask_file'] = ""
    namelist['number_special_points'] = 0
    namelist['lflake_correction'] = ".TRUE."

    if args['igrid_type'] == 1:  # ICON
        namelist['tile_mode'] = 1
    elif args['igrid_type'] == 2:  # COSMO
        namelist['tile_mode'] = 0

    if args['use_array_cache']:
        namelist['l_use_array_cache'] = ".TRUE."
    else:
        namelist['l_use_array_cache'] = ".FALSE."

    return namelist


def setup_namelist(args) -> dict:

    namelist = {}

    namelist.update(setup_oro_namelist(args))
    namelist.update(setup_albedo_namelist(args))
    namelist.update(setup_aot_namelist(args))
    namelist.update(setup_tclim_namelist(args))
    namelist.update(setup_lu_namelist(args))
    namelist.update(setup_flake_namelist(args))
    namelist.update(setup_ndvi_namelist(args))
    namelist.update(setup_urban_namelist(args))
    namelist.update(setup_soil_namelist(args))
    namelist.update(setup_era_namelist(args))
    namelist.update(setup_emiss_namelist(args))
    namelist.update(setup_cdnc_namelist(args))
    namelist.update(setup_edgar_namelist(args))
    namelist.update(setup_art_namelist(args))
    namelist.update(setup_check_namelist(args))

    return namelist


def setup_runscript(args):
    runscript = {}

    runscript['account'] = args['account']
    runscript['pythonpath'] = os.path.join(os.getcwd(), 'lib')

    executables = [
        '"extpar_landuse_to_buffer.exe" ', '"extpar_topo_to_buffer.exe" ',
        '"extpar_cru_to_buffer.py" ', '"extpar_aot_to_buffer.py" ',
        '"extpar_flake_to_buffer.exe" ', '"extpar_soil_to_buffer.exe" ',
        '"extpar_alb_to_buffer.py" ', '"extpar_ndvi_to_buffer.py" '
    ]

    if args['lurban']:
        executables.append('"extpar_ahf_to_buffer.py" ')
        executables.append('"extpar_isa_to_buffer.py" ')

    if args['igrid_type'] == 1:
        executables.append('"extpar_era_to_buffer.py" ')
        executables.append('"extpar_emiss_to_buffer.py" ')

        # ICON only executables
        if args['enable_cdnc']:
            executables.append('"extpar_cdnc_to_buffer.py" ')
        if args['enable_edgar']:
            executables.append('"extpar_edgar_to_buffer.py" ')
        if args['enable_art']:
            executables.append('"extpar_art_to_buffer.py" ')

    executables.append('"extpar_consistency_check.exe" ')

    runscript['extpar_executables'] = executables

    return runscript


def replace_placeholders(args, templates, dir, actual_values):
    all_templates = {}

    # read templates
    for template in templates:
        with open(os.path.join(dir, template), 'r') as f:
            all_templates[template] = f.read()

    # replace all @PLACEHOLDERS@ with real values
    for key, value in actual_values.items():
        key = f'@{key.upper()}@'
        for template in templates:
            if isinstance(value, list):
                all_templates[template] = all_templates[template].replace(
                    key, str("".join(value)))
            else:
                if value is not None:
                    all_templates[template] = all_templates[template].replace(
                        key, str(value))
                else:
                    raise ValueError(
                        f'The placeholder {key} in {all_templates[template]} was replaced with None.'
                        'This likely means that it was not specified in the config file and has no default value.'
                    )

    # check that no @PLACEHOLDERS@ are left
    for template in templates:
        if '@' in all_templates[template]:
            raise ValueError(
                f'Not all placeholders in {all_templates[template]} were replaced'
            )

    # write complete template to file
    for template in templates:
        # special case for namelist.py, which is a python file but not valid with placeholders
        file = template
        if template == 'namelist':
            file = 'namelist.py'
        with open(os.path.join(args['run_dir'], file), 'w') as f:
            f.write(all_templates[template])
        logging.info(f'{template} written to {args["run_dir"]}')


def extend_cosmo_grid_for_radtopo(run_dir: str, tg: CosmoGrid):

    circum_earth = 40075160.0
    horizon_radius = 40000.0
    res_in = tg.dlon * (circum_earth / 360.0)
    nborder = max(int(horizon_radius / res_in), 4)

    startlon_tot = tg.startlon_tot - nborder * tg.dlon
    startlat_tot = tg.startlat_tot - nborder * tg.dlat
    ie_tot = tg.ie_tot + 2 * nborder
    je_tot = tg.je_tot + 2 * nborder

    extended_grid = os.path.join(run_dir, '.extended_grid_radtopo')

    with open(extended_grid, 'w') as f:

        f.write(f'&lmgrid\n')
        f.write(f'pollon = {tg.pollon}\n')
        f.write(f'pollat = {tg.pollat}\n')
        f.write(f'dlon = {tg.dlon}\n')
        f.write(f'dlat = {tg.dlat}\n')
        f.write(f'startlon_tot = {startlon_tot}\n')
        f.write(f'startlat_tot = {startlat_tot}\n')
        f.write(f'ie_tot = {ie_tot}\n')
        f.write(f'je_tot = {je_tot}\n')
        f.write(f'/\n')

    return CosmoGrid(extended_grid)


def compute_merit_tiles(lonmax: float, lonmin: float, latmax: float,
                        latmin: float, lsgsl: bool) -> dict:

    name_lon = [
        'W180-W150', 'W150-W120', 'W120-W090', 'W090-W060', 'W060-W030',
        'W030-E000', 'E000-E030', 'E030-E060', 'E060-E090', 'E090-E120',
        'E120-E150', 'E150-E180'
    ]

    name_lat = [
        'N90-N60', 'N60-N30', 'N30-N00', 'N00-S30', 'S30-S60', 'S60-S90'
    ]

    prefix_lat = ['MERIT', 'MERIT', 'MERIT', 'MERIT', 'MERIT', 'REMA_BKG']

    merit_tiles_lon = np.empty([12, 6])
    merit_tiles_lat = np.empty([12, 6])

    merit_lon = -180.0
    merit_lat = 90.0
    for j in range(0, 6):
        for i in range(0, 12):
            merit_tiles_lon[i, j] = merit_lon + float(i * 30)
            merit_tiles_lat[i, j] = merit_lat - float(j * 30)

    ilon_min = 0
    ilon_max = 0
    ilat_min = 0
    ilat_max = 0

    for j in range(0, 6):
        for i in range(0, 12):
            if merit_tiles_lon[i, j] < lonmin:
                ilon_min = i
            if merit_tiles_lon[i, j] < lonmax:
                ilon_max = i
            if merit_tiles_lat[i, j] > latmin:
                ilat_max = j
            if merit_tiles_lat[i, j] > latmax:
                ilat_min = j

    ntiles_column = ilon_max - ilon_min + 1
    ntiles_row = ilat_max - ilat_min + 1

    name_tiles = []
    for j in range(ilat_min, ilat_max + 1):
        for i in range(ilon_min, ilon_max + 1):
            name_tiles.append(
                f"'{prefix_lat[j]}_{name_lat[j]}_{name_lon[i]}.nc' ")

    namelist = {}
    namelist['ntiles_column'] = ntiles_column
    namelist['ntiles_row'] = ntiles_row
    namelist['topo_files'] = name_tiles

    if lsgsl:
        namelist['sgsl_files'] = [
            f"'S_ORO_{tile_nr}' " for tile_nr in range(0, len(name_tiles))
        ]
        namelist['lpreproc_oro'] = ".TRUE."

    return namelist


def compute_aster_tiles(lonmax: float, lonmin: float, latmax: float,
                        latmin: float, lsgsl: bool) -> dict:

    # safety check
    if latmax > 60.0 or latmax < -60.0:
        logging.error('Domains using Aster cannot exceed 60 N or 60 S')
        raise ValueError('Domains using Aster cannot exceed 60 N or 60 S')

    aster_tiles_lon = np.empty([12, 20])
    aster_tiles_lat = np.empty([12, 20])

    aster_lon = -180.0
    aster_lat = 60.0
    for j in range(0, 20):
        for i in range(0, 12):
            aster_tiles_lon[i, j] = aster_lon + float(i * 30)
            aster_tiles_lat[i, j] = aster_lat - float(j * 6)

    ilon_min = 0
    ilon_max = 0
    ilat_min = 0
    ilat_max = 0

    for j in range(0, 20):
        for i in range(0, 12):
            if aster_tiles_lon[i, j] < lonmin:
                ilon_min = i
            if aster_tiles_lon[i, j] < lonmax:
                ilon_max = i
            if aster_tiles_lat[i, j] > latmin:
                ilat_max = j
            if aster_tiles_lat[i, j] > latmax:
                ilat_min = j

    ntiles_column = ilon_max - ilon_min + 1
    ntiles_row = ilat_max - ilat_min + 1

    aster_files = np.empty(240, int)

    icount = 0
    for j in range(ilat_min, ilat_max + 1):
        for i in range(ilon_min, ilon_max + 1):
            aster_files[icount] = int(1 + i + 12 * j)
            icount += 1

    namelist = {}
    namelist['ntiles_column'] = ntiles_column
    namelist['ntiles_row'] = ntiles_row
    namelist['topo_files'] = [
        f"'ASTER_orig_T{aster_files[idx-1]:03}.nc' "
        for idx in range(1, icount + 1)
    ]

    if lsgsl:
        namelist['sgsl_files'] = [
            f"'S_ORO_T{aster_files[idx-1]:03}.nc' "
            for idx in range(1, icount + 1)
        ]
        namelist['lpreproc_oro'] = ".FALSE."

    return namelist


if __name__ == '__main__':

    main()
