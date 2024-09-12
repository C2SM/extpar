#!/usr/bin/env python3
import argparse
import os
import stat
import sys
import shutil
import logging

import numpy as np

# extpar modules from lib
try:
    from extpar.lib.grid_def import CosmoGrid
    from extpar.lib.utilities import launch_shell
    DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
except ImportError:  # package not installed -> use PYTHONPATH
    from grid_def import CosmoGrid
    from utilities import launch_shell
    DATA_DIR = ".."


def main():

    # initialize logger
    logging.basicConfig(level=logging.INFO, format='%(message)s')

    parser = argparse.ArgumentParser(description='Process some integers.')

    parser.add_argument('--input_cosmo_grid',
                        type=str,
                        required=True,
                        help='Fortran Namelist "INPUT_COSMO_GRID"')
    parser.add_argument(
        '--iaot_type',
        type=int,
        required=True,
        help='1: Global Aerosol Climatology Project, 2: AeroCom1')
    parser.add_argument('--ilu_type',
                        type=int,
                        required=True,
                        help='1: GLOBCOVER 2: GLC2000')
    parser.add_argument('--ialb_type',
                        type=int,
                        required=True,
                        help='1: VIS,UV,NIR 2: SOIL 3: VIS')
    parser.add_argument(
        '--isoil_type',
        type=int,
        required=True,
        help=
        '1: FAO Digital Soil Map of the World, 2: Harmonized World Soil Database (with additional data for deepsoil): Harmonized World Soil Database (HWSD)'
    )
    parser.add_argument('--itopo_type',
                        type=int,
                        required=True,
                        help='1: GLOBE, 2: ASTER, 3: MERIT')
    parser.add_argument('--lsgsl',
                        action='store_true',
                        help='Compute subgrid-scale slope parameter (S_ORO)')
    parser.add_argument('--lfilter_oro',
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
    parser.add_argument('--run_dir',
                        type=str,
                        required=True,
                        help='Folder for running Extpar')
    parser.add_argument('--account',
                        type=str,
                        required=True,
                        help='Account for slurm job')
    parser.add_argument('--host',
                        type=str,
                        required=True,
                        choices=('levante'),
                        help='Host')
    parser.add_argument('--no_batch_job',
                        action='store_true',
                        help="Run jobscript not as batch job")

    args = parser.parse_args()

    generate_external_parameters(args.input_cosmo_grid, args.iaot_type,
                                 args.ilu_type, args.ialb_type,
                                 args.isoil_type, args.itopo_type,
                                 args.raw_data_path, args.run_dir,
                                 args.account, args.host, args.no_batch_job,
                                 args.lurban, args.lsgsl, args.lfilter_oro)


def generate_external_parameters(input_cosmo_grid,
                                 iaot_type,
                                 ilu_type,
                                 ialb_type,
                                 isoil_type,
                                 itopo_type,
                                 raw_data_path,
                                 run_dir,
                                 account,
                                 host,
                                 no_batch_job=False,
                                 lurban=False,
                                 lsgsl=False,
                                 lfilter_oro=False):

    # initialize logger
    logging.basicConfig(level=logging.INFO, format='%(message)s')

    # make all paths absolut
    raw_data_path = os.path.abspath(raw_data_path)
    run_dir = os.path.abspath(run_dir)
    input_cosmo_grid = os.path.abspath(input_cosmo_grid)

    args_dict = {
        'input_cosmo_grid': input_cosmo_grid,
        'iaot_type': iaot_type,
        'ilu_type': ilu_type,
        'ialb_type': ialb_type,
        'isoil_type': isoil_type,
        'itopo_type': itopo_type,
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


def prepare_sandbox(args, namelist, runscript):

    os.makedirs(args['run_dir'], exist_ok=True)
    write_namelist(args, namelist)
    write_runscript(args, runscript)
    copy_required_files(args, runscript['extpar_executables'])


def write_runscript(args, runscript):
    dir = os.path.join(DATA_DIR, "templates")
    files = [f'submit.{args["host"]}.sh']

    replace_placeholders(args, files, dir, runscript)


def write_namelist(args, namelist):
    templates_dir = os.path.join(DATA_DIR, "templates")
    files = [
        'INPUT_ORO', 'INPUT_RADTOPO', 'INPUT_OROSMOOTH', 'INPUT_SGSL',
        'INPUT_AOT', 'INPUT_LU', 'INPUT_FLAKE', 'INPUT_SCALE_SEP',
        'INPUT_SOIL', 'INPUT_CHECK', 'namelist.py'
    ]

    replace_placeholders(args, files, templates_dir, namelist)

    # INPUT_grid_org
    with open(os.path.join(args['run_dir'], 'INPUT_grid_org'), 'w') as f:
        f.write('&GRID_DEF \n')
        f.write('igrid_type = 2 \n')
        f.write("domain_def_namelist = 'INPUT_COSMO_GRID' \n")
        f.write('/ \n')

    # INPUT_COSMO_GRID
    shutil.copy(args['input_cosmo_grid'],
                os.path.join(args['run_dir'], 'INPUT_COSMO_GRID'))


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

    tg = CosmoGrid(args['input_cosmo_grid'])

    namelist = {}

    if tg.dlon < 0.05 and tg.dlat < 0.05:
        lradtopo = True
    else:
        lradtopo = False

    # &orography_io_extpar
    namelist['orography_buffer_file'] = 'oro_buffer.nc'
    namelist['orography_output_file'] = 'oro_cosmo.nc'

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

    if args['itopo_type'] == 1:
        namelist['topo_files'] = [
            f"'GLOBE_{letter.upper()}10.nc' "
            for letter in list(map(chr, range(ord('a'),
                                              ord('p') + 1)))
        ]
        namelist['ntiles_column'] = 4
        namelist['ntiles_row'] = 4

        if args['lsgsl']:
            namelist['sgsl_files'] = [
                f"'S_ORO_{letter.upper()}10.nc' "
                for letter in list(map(chr, range(ord('a'),
                                                  ord('p') + 1)))
            ]
            namelist['lpreproc_oro'] = ".FALSE."

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
            compute_aster_tiles(args["run_dir"], tg, args['lsgsl'], lradtopo))
        namelist['lscale_separation'] = ".FALSE."
        namelist['scale_sep_files'] = "'placeholder_file'"
        namelist['lsso_param'] = ".TRUE."
    elif args['itopo_type'] == 3:
        namelist.update(
            compute_merit_tiles(args["run_dir"], tg, args['lsgsl'], lradtopo))
        namelist['lscale_separation'] = ".FALSE."
        namelist['scale_sep_files'] = "'placeholder_file'"
        namelist['lsso_param'] = ".TRUE."

    # &scale_separated_raw_data
    # other paramters of the namelist already set earlier
    namelist['raw_data_scale_sep_path'] = args['raw_data_path']

    # &orography_smoothing
    if args['lfilter_oro']:
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
    if lradtopo:
        namelist['lradtopo'] = ".TRUE."
    else:
        namelist['lradtopo'] = ".FALSE."

    # &sgsl_raw_data
    namelist['raw_data_sgsl_path'] = args['raw_data_path']
    namelist['idem_type'] = args['itopo_type']

    return namelist


def setup_lu_namelist(args):
    namelist = {}
    namelist['i_landuse_data'] = args['ilu_type']
    namelist['ilookup_table_lu'] = args['ilu_type']
    namelist['raw_data_lu_path'] = args['raw_data_path']
    namelist['raw_data_glcc_path'] = args['raw_data_path']
    namelist['lu_buffer_file'] = 'lu_buffer.nc'
    namelist['raw_data_glcc_filename'] = 'GLCC_usgs_class_byte.nc'
    namelist['glcc_buffer_file'] = 'glcc_buffer.nc'
    namelist['l_use_corine'] = ".FALSE."
    if args['ilu_type'] == 1:
        namelist['raw_data_lu_filename'] = [
            f"'GLOBCOVER_{i}_16bit.nc' " for i in range(0, 6)
        ]
    elif args['ilu_type'] == 2:
        # we need "" padding for correct replacement in Fortran namelist
        namelist['raw_data_lu_filename'] = "'GLC2000_byte.nc'"

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

    return namelist


def setup_soil_namelist(args):
    namelist = {}

    # &soil_raw_data
    namelist['isoil_data'] = args['isoil_type']
    namelist['raw_data_soil_path'] = args['raw_data_path']

    if args['isoil_type'] == 1:
        namelist['raw_data_soil_filename'] = 'FAO_DSMW_double.nc'
        namelist['ldeep_soil'] = ".FALSE"
    elif args['isoil_type'] == 2:
        namelist['raw_data_soil_filename'] = 'HWSD0_30_topsoil.nc'
        namelist['raw_data_deep_soil_filename'] = 'HWSD30_100_subsoil.nc'
        namelist['ldeep_soil'] = ".TRUE"
    elif args['isoil_type'] == 3:
        namelist['raw_data_soil_filename'] = 'HWSD0_30_topsoil.nc'
        namelist['raw_data_deep_soil_filename'] = 'HWSD30_100_subsoil.nc'
        namelist['ldeep_soil'] = ".FALSE"

    # &soil_io_extpar
    namelist['soil_buffer_file'] = 'soil_buffer.nc'

    # &HWSD_index_files

    # Quickfix, should be in linked_data as well
    namelist['path_HWSD_index_files'] = os.path.join(args['raw_data_path'],
                                                     '../soil')
    namelist['lookup_table_HWSD'] = 'LU_TAB_HWSD_UF.data'
    namelist['HWSD_data'] = 'HWSD_DATA_COSMO.data'
    namelist['HWSD_data_deep'] = 'HWSD_DATA_COSMO_S.data'

    return namelist


def setup_tclim_namelist(args):
    namelist = {}

    namelist['it_cl_type'] = 1
    namelist['raw_data_t_clim_path'] = args['raw_data_path']
    namelist['raw_data_tclim_coarse'] = 'absolute_hadcrut3.nc'
    namelist['raw_data_tclim_fine'] = 'CRU_T_SOIL_clim.nc'
    namelist['t_clim_buffer_file'] = 'tclim_buffer.nc'

    return namelist


def setup_flake_namelist(args):
    namelist = {}

    namelist['raw_data_flake_path'] = args['raw_data_path']
    namelist['raw_data_flake_filename'] = 'GLDB_lakedepth.nc'
    namelist['flake_buffer_file'] = 'flake_buffer.nc'

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

    return namelist


def setup_ndvi_namelist(args):
    namelist = {}

    namelist['raw_data_ndvi_path'] = args['raw_data_path']
    namelist['raw_data_ndvi_filename'] = 'NDVI_1998_2003.nc'
    namelist['ndvi_buffer_file'] = 'ndvi_buffer.nc'

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


def setup_check_namelist(args):
    namelist = {}

    namelist['netcdf_output_filename'] = 'external_parameter.nc'
    namelist['i_lsm_data'] = 1
    namelist['land_sea_mask_file'] = ""
    namelist['number_special_points'] = 0
    namelist['lflake_correction'] = ".TRUE."

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
    namelist.update(setup_check_namelist(args))

    return namelist


def setup_runscript(args):
    runscript = {}

    runscript['account'] = args['account']
    runscript['pythonpath'] = os.path.join(os.getcwd(), 'lib')

    executables = [
        '"extpar_landuse_to_buffer.exe" ', '"extpar_topo_to_buffer.exe" ',
        '"extpar_cru_to_buffer.py" ', '"extpar_aot_to_buffer.exe" ',
        '"extpar_flake_to_buffer.exe" ', '"extpar_soil_to_buffer.exe" ',
        '"extpar_alb_to_buffer.py" ', '"extpar_ndvi_to_buffer.py" '
    ]

    if args['lurban']:
        executables.append('"extpar_ahf_to_buffer.py" ')
        executables.append('"extpar_isa_to_buffer.py" ')

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
                all_templates[template] = all_templates[template].replace(
                    key, str(value))

    # write complete template to file
    for template in templates:
        with open(os.path.join(args['run_dir'], template), 'w') as f:
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


def compute_merit_tiles(run_dir: str, tg: CosmoGrid, lsgsl: bool,
                        lradtopo: bool) -> dict:

    name_lon = [
        'W180-W150', 'W150-W120', 'W120-W090', 'W090-W060', 'W060-W030',
        'W030-E000', 'E000-E030', 'E030-E060', 'E060-E090', 'E090-E120',
        'E120-E150', 'E150-E180'
    ]

    name_lat = [
        'N90-N60', 'N60-N30', 'N30-N00', 'N00-S30', 'S30-S60', 'S60-S90'
    ]

    prefix_lat = ['MERIT', 'MERIT', 'MERIT', 'MERIT', 'MERIT', 'REMA']
    if lradtopo:
        tg = extend_cosmo_grid_for_radtopo(run_dir, tg)

    zlonmax = np.amax(tg.lons)
    zlonmin = np.amin(tg.lons)
    zlatmin = np.amin(tg.lats)
    zlatmax = np.amax(tg.lats)

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
            if merit_tiles_lon[i, j] < zlonmin:
                ilon_min = i
            if merit_tiles_lon[i, j] < zlonmax:
                ilon_max = i
            if merit_tiles_lat[i, j] > zlatmin:
                ilat_max = j
            if merit_tiles_lat[i, j] > zlatmax:
                ilat_min = j

    ntiles_column = ilon_max - ilon_min + 1
    ntiles_row = ilat_max - ilat_min + 1

    merit_files = np.empty(72, int)

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


def compute_aster_tiles(run_dir: str, tg: CosmoGrid, lsgsl: bool,
                        lradtopo: bool) -> dict:

    if lradtopo:
        tg = extend_cosmo_grid_for_radtopo(run_dir, tg)

    zlonmax = np.amax(tg.lons)
    zlonmin = np.amin(tg.lons)
    zlatmin = np.amin(tg.lats)
    zlatmax = np.amax(tg.lats)

    # safety check
    if zlatmax > 60.0 or zlatmax < -60.0:
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
            if aster_tiles_lon[i, j] < zlonmin:
                ilon_min = i
            if aster_tiles_lon[i, j] < zlonmax:
                ilon_max = i
            if aster_tiles_lat[i, j] > zlatmin:
                ilat_max = j
            if aster_tiles_lat[i, j] > zlatmax:
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
