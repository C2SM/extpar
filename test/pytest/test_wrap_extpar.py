import os
import tempfile
import shutil
import requests
import pytest

try:
    from extpar.WrapExtpar import *
except ImportError:
    import sys
    sys.path.append(
        os.path.join(os.path.dirname(__file__), '..', '..', 'python'))
    sys.path.append(
        os.path.join(os.path.dirname(__file__), '..', '..', 'python/lib'))
    from WrapExtpar import *


@pytest.fixture
def tmp_dir():
    tmp_dir = tempfile.mkdtemp()
    yield tmp_dir
    shutil.rmtree(tmp_dir)


@pytest.fixture
def icon_grid(tmp_dir):

    # URL of the file to download
    url = 'http://icon-downloads.mpimet.mpg.de/grids/public/mpim/0013/icon_grid_0013_R02B04_G.nc'

    # Local file path to save the downloaded file
    output_file = os.path.join(tmp_dir, 'icon_grid_0013_R02B04_G.nc')

    # Download the file
    response = requests.get(url, stream=True)
    response.raise_for_status()  # Raise an error for HTTP codes 4xx/5xx

    # Write the content to a local file
    with open(output_file, 'wb') as file:
        for chunk in response.iter_content(chunk_size=8192):
            file.write(chunk)
    yield output_file


def test_setup_flake_namelist():
    args = {'raw_data_path': '/path/to/raw/data'}
    expected_namelist = {
        'raw_data_flake_path': '/path/to/raw/data',
        'raw_data_flake_filename': 'GLDB_lakedepth.nc',
        'flake_buffer_file': 'flake_buffer.nc'
    }

    result = setup_flake_namelist(args)

    assert result == expected_namelist


def test_setup_art_namelist():
    args = {'raw_data_path': '/path/to/raw/data'}
    expected_namelist = {
        'raw_data_art_path': '/path/to/raw/data',
        'raw_data_art_filename': 'HWSD0_USDA.nc',
        'art_buffer_file': 'art_buffer.nc'
    }

    result = setup_art_namelist(args)

    assert result == expected_namelist


def test_setup_tclim_1_namelist():
    args = {
        'raw_data_path': '/path/to/raw/data',
        'it_cl_type': 1,
        'tcorr_lapse_rate': 0.0065,
        'tcorr_offset': 0.0
    }
    expected_output = {
        'it_cl_type': 1,
        'raw_data_t_clim_path': '/path/to/raw/data',
        'raw_data_tclim_coarse': 'absolute_hadcrut3.nc',
        'raw_data_tclim_fine': 'CRU_T_SOIL_clim.nc',
        't_clim_buffer_file': 'tclim_buffer.nc',
        'tcorr_lapse_rate': 0.0065,
        'tcorr_offset': 0.0
    }

    result = setup_tclim_namelist(args)
    assert result == expected_output


def test_setup_tclim_2_namelist():
    args = {
        'raw_data_path': '/path/to/raw/data',
        'it_cl_type': 2,
        'tcorr_lapse_rate': 0.0065,
        'tcorr_offset': 0.0
    }
    expected_output = {
        'it_cl_type': 2,
        'raw_data_t_clim_path': '/path/to/raw/data',
        'raw_data_tclim_coarse': 'absolute_hadcrut3.nc',
        'raw_data_tclim_fine': 'CRU_T_SOIL_clim.nc',
        't_clim_buffer_file': 'tclim_buffer.nc',
        'tcorr_lapse_rate': 0.0065,
        'tcorr_offset': 0.0
    }

    result = setup_tclim_namelist(args)
    assert result == expected_output


def test_setup_tclim_2_custom_tcorrection_namelist():
    args = {
        'raw_data_path': '/path/to/raw/data',
        'it_cl_type': 2,
        'tcorr_lapse_rate': 0.0051,
        'tcorr_offset': 2.0
    }
    expected_output = {
        'it_cl_type': 2,
        'raw_data_t_clim_path': '/path/to/raw/data',
        'raw_data_tclim_coarse': 'absolute_hadcrut3.nc',
        'raw_data_tclim_fine': 'CRU_T_SOIL_clim.nc',
        't_clim_buffer_file': 'tclim_buffer.nc',
        'tcorr_lapse_rate': 0.0051,
        'tcorr_offset': 2.0
    }

    result = setup_tclim_namelist(args)
    assert result == expected_output


def test_setup_tclim_unknown_type():
    args = {
        'raw_data_path': '/path/to/raw/data',
        'it_cl_type': 9,
        'tcorr_lapse_rate': 0.0065,
        'tcorr_offset': 0.0
    }
    with pytest.raises(ValueError, match='Unknown it_cl_type 9'):
        setup_tclim_namelist(args)


def test_setup_albedo_namelist_type_1():
    args = {'raw_data_path': '/path/to/data', 'ialb_type': 1}
    expected_namelist = {
        'raw_data_alb_path': '/path/to/data',
        'ialb_type': 1,
        'alb_buffer_file': 'alb_buffer.nc',
        'raw_data_alb_filename': 'alb_new.nc',
        'raw_data_alnid_filename': 'alnid_new.nc',
        'raw_data_aluvd_filename': 'aluvd_new.nc'
    }
    assert setup_albedo_namelist(args) == expected_namelist


def test_setup_albedo_namelist_type_2():
    args = {'raw_data_path': '/path/to/data', 'ialb_type': 2}
    expected_namelist = {
        'raw_data_alb_path': '/path/to/data',
        'ialb_type': 2,
        'alb_buffer_file': 'alb_buffer.nc',
        'raw_data_alb_filename': 'global_soil_albedo.nc'
    }
    assert setup_albedo_namelist(args) == expected_namelist


def test_setup_albedo_namelist_type_3():
    args = {'raw_data_path': '/path/to/data', 'ialb_type': 3}
    expected_namelist = {
        'raw_data_alb_path': '/path/to/data',
        'ialb_type': 3,
        'alb_buffer_file': 'alb_buffer.nc',
        'raw_data_alb_filename': 'alb_new.nc'
    }
    assert setup_albedo_namelist(args) == expected_namelist


def test_setup_ndvi_namelist():
    args = {'raw_data_path': '/path/to/data'}
    expected_namelist = {
        'raw_data_ndvi_path': '/path/to/data',
        'raw_data_ndvi_filename': 'NDVI_1998_2003.nc',
        'ndvi_buffer_file': 'ndvi_buffer.nc'
    }
    assert setup_ndvi_namelist(args) == expected_namelist


def test_setup_urban_namelist():
    args = {'raw_data_path': '/path/to/data'}
    expected_namelist = {
        'iahf_type': 1,
        'raw_data_ahf_path': '/path/to/data',
        'raw_data_ahf_filename': 'AHF_2006_2.5min_lonlat.nc',
        'ahf_buffer_file': 'ahf_buffer.nc',
        'isa_type': 1,
        'raw_data_isa_path': '/path/to/data',
        'raw_data_isa_filename': 'NOAA_ISA_16bit_lonlat.nc',
        'isa_buffer_file': 'isa_buffer.nc'
    }
    assert setup_urban_namelist(args) == expected_namelist


def test_setup_check_namelist():
    args = {'igrid_type': 1, 'use_array_cache': False}
    expected_namelist = {
        'netcdf_output_filename': 'external_parameter.nc',
        'i_lsm_data': 1,
        'land_sea_mask_file': "",
        'number_special_points': 0,
        'lflake_correction': ".TRUE.",
        'tile_mode': 1,
        'l_use_array_cache': ".FALSE."
    }
    assert setup_check_namelist(args) == expected_namelist


def test_setup_check_namelist_array_cache():
    args = {'igrid_type': 1, 'use_array_cache': True}
    expected_namelist = {
        'netcdf_output_filename': 'external_parameter.nc',
        'i_lsm_data': 1,
        'land_sea_mask_file': "",
        'number_special_points': 0,
        'lflake_correction': ".TRUE.",
        'tile_mode': 1,
        'l_use_array_cache': ".TRUE."
    }
    assert setup_check_namelist(args) == expected_namelist


def test_setup_check_namelist_tile_mode():
    args = {'igrid_type': 2, 'use_array_cache': False}
    expected_namelist = {
        'netcdf_output_filename': 'external_parameter.nc',
        'i_lsm_data': 1,
        'land_sea_mask_file': "",
        'number_special_points': 0,
        'lflake_correction': ".TRUE.",
        'tile_mode': 0,
        'l_use_array_cache': ".FALSE."
    }
    assert setup_check_namelist(args) == expected_namelist


def test_setup_lu_namelist_type_1():
    args = {
        'ilu_type': 1,
        'raw_data_path': '/path/to/data',
        'l_use_corine': False,
        'ilookup_table_lu': 1
    }
    expected_namelist = {
        'i_landuse_data':
        1,
        'ilookup_table_lu':
        1,
        'raw_data_lu_path':
        '/path/to/data',
        'raw_data_glcc_path':
        '/path/to/data',
        'lu_buffer_file':
        'lu_buffer.nc',
        'raw_data_glcc_filename':
        'GLCC_usgs_class_byte.nc',
        'glcc_buffer_file':
        'glcc_buffer.nc',
        'ntiles_globcover':
        6,
        'l_terra_urb':
        ".FALSE.",
        'l_use_corine':
        ".FALSE.",
        'raw_data_lu_filename': [
            "'GLOBCOVER_0_16bit.nc' ", "'GLOBCOVER_1_16bit.nc' ",
            "'GLOBCOVER_2_16bit.nc' ", "'GLOBCOVER_3_16bit.nc' ",
            "'GLOBCOVER_4_16bit.nc' ", "'GLOBCOVER_5_16bit.nc' "
        ]
    }
    assert setup_lu_namelist(args) == expected_namelist


def test_setup_lu_namelist_corine():
    args = {
        'ilu_type': 1,
        'raw_data_path': '/path/to/data',
        'l_use_corine': True,
        'ilookup_table_lu': 1
    }
    expected_namelist = {
        'i_landuse_data': 1,
        'ilookup_table_lu': 1,
        'raw_data_lu_path': '/path/to/data',
        'raw_data_glcc_path': '/path/to/data',
        'lu_buffer_file': 'lu_buffer.nc',
        'raw_data_glcc_filename': 'GLCC_usgs_class_byte.nc',
        'glcc_buffer_file': 'glcc_buffer.nc',
        'ntiles_globcover': 1,
        'l_terra_urb': ".FALSE.",
        'l_use_corine': ".TRUE.",
        'raw_data_lu_filename': "'CORINE_globcover.nc'"
    }
    assert setup_lu_namelist(args) == expected_namelist


def test_setup_lu_namelist_type_2():
    args = {
        'ilu_type': 2,
        'raw_data_path': '/path/to/data',
        'l_use_corine': False,
        'ilookup_table_lu': 2
    }
    expected_namelist = {
        'i_landuse_data': 2,
        'ilookup_table_lu': 2,
        'raw_data_lu_path': '/path/to/data',
        'raw_data_glcc_path': '/path/to/data',
        'lu_buffer_file': 'lu_buffer.nc',
        'raw_data_glcc_filename': 'GLCC_usgs_class_byte.nc',
        'glcc_buffer_file': 'glcc_buffer.nc',
        'ntiles_globcover': 6,
        'l_terra_urb': ".FALSE.",
        'l_use_corine': ".FALSE.",
        'raw_data_lu_filename': "'GLC2000_byte.nc'"
    }
    assert setup_lu_namelist(args) == expected_namelist


def test_setup_lu_namelist_type_6():
    args = {
        'ilu_type': 6,
        'raw_data_path': '/path/to/data',
        'l_use_corine': False,
        'ilookup_table_lu': 1
    }
    expected_namelist = {
        'i_landuse_data': 6,
        'ilookup_table_lu': 1,
        'raw_data_lu_path': '/path/to/data',
        'raw_data_glcc_path': '/path/to/data',
        'lu_buffer_file': 'lu_buffer.nc',
        'raw_data_glcc_filename': 'GLCC_usgs_class_byte.nc',
        'glcc_buffer_file': 'glcc_buffer.nc',
        'ntiles_globcover': 6,
        'l_terra_urb': ".TRUE.",
        'l_use_corine': ".FALSE.",
        'raw_data_lu_filename': "'ECOCLIMAP_SG.nc'"
    }
    assert setup_lu_namelist(args) == expected_namelist


def test_setup_aot_namelist_type_1():
    args = {'iaot_type': 1, 'raw_data_path': '/path/to/data'}
    expected_namelist = {
        'iaot_type': 1,
        'raw_data_aot_path': '/path/to/data',
        'aot_buffer_file': 'aot_buffer.nc',
        'raw_data_aot_filename': 'aot_GACP.nc'
    }
    assert setup_aot_namelist(args) == expected_namelist


def test_setup_aot_namelist_type_2():
    args = {'iaot_type': 2, 'raw_data_path': '/path/to/data'}
    expected_namelist = {
        'iaot_type': 2,
        'raw_data_aot_path': '/path/to/data',
        'aot_buffer_file': 'aot_buffer.nc',
        'raw_data_aot_filename': 'aod_AeroCom1.nc'
    }
    assert setup_aot_namelist(args) == expected_namelist


def test_generate_globe_filenames():
    expected_filenames = [
        "'GLOBE_A10.nc' ", "'GLOBE_B10.nc' ", "'GLOBE_C10.nc' ",
        "'GLOBE_D10.nc' ", "'GLOBE_E10.nc' ", "'GLOBE_F10.nc' ",
        "'GLOBE_G10.nc' ", "'GLOBE_H10.nc' ", "'GLOBE_I10.nc' ",
        "'GLOBE_J10.nc' ", "'GLOBE_K10.nc' ", "'GLOBE_L10.nc' ",
        "'GLOBE_M10.nc' ", "'GLOBE_N10.nc' ", "'GLOBE_O10.nc' ",
        "'GLOBE_P10.nc' "
    ]
    assert generate_globe_filenames() == expected_filenames


def test_setup_runscript_with_urban_cosmo():
    args = {
        'account': 'test_account',
        'lurban': True,
        'igrid_type': 2,
        'enable_cdnc': False,
        'enable_art': True,
        'enable_edgar': False
    }
    expected_runscript = {
        'account':
        'test_account',
        'pythonpath':
        os.path.join(os.getcwd(), 'lib'),
        'extpar_executables': [
            '"extpar_landuse_to_buffer.exe" ', '"extpar_topo_to_buffer.exe" ',
            '"extpar_cru_to_buffer.py" ', '"extpar_aot_to_buffer.py" ',
            '"extpar_flake_to_buffer.exe" ', '"extpar_soil_to_buffer.exe" ',
            '"extpar_alb_to_buffer.py" ', '"extpar_ndvi_to_buffer.py" ',
            '"extpar_ahf_to_buffer.py" ', '"extpar_isa_to_buffer.py" ',
            '"extpar_consistency_check.exe" '
        ]
    }
    assert setup_runscript(args) == expected_runscript


def test_setup_runscript_without_urban_icon():
    args = {
        'account': 'test_account',
        'lurban': False,
        'igrid_type': 1,
        'enable_cdnc': False,
        'enable_art': False,
        'enable_edgar': False
    }
    expected_runscript = {
        'account':
        'test_account',
        'pythonpath':
        os.path.join(os.getcwd(), 'lib'),
        'extpar_executables': [
            '"extpar_landuse_to_buffer.exe" ', '"extpar_topo_to_buffer.exe" ',
            '"extpar_cru_to_buffer.py" ', '"extpar_aot_to_buffer.py" ',
            '"extpar_flake_to_buffer.exe" ', '"extpar_soil_to_buffer.exe" ',
            '"extpar_alb_to_buffer.py" ', '"extpar_ndvi_to_buffer.py" ',
            '"extpar_era_to_buffer.py" ', '"extpar_emiss_to_buffer.py" ',
            '"extpar_consistency_check.exe" '
        ]
    }
    assert setup_runscript(args) == expected_runscript


def test_setup_runscript_without_urban_with_edgar_and_cdnc_icon():
    args = {
        'account': 'test_account',
        'lurban': False,
        'igrid_type': 1,
        'enable_cdnc': True,
        'enable_art': False,
        'enable_edgar': True
    }
    expected_runscript = {
        'account':
        'test_account',
        'pythonpath':
        os.path.join(os.getcwd(), 'lib'),
        'extpar_executables': [
            '"extpar_landuse_to_buffer.exe" ', '"extpar_topo_to_buffer.exe" ',
            '"extpar_cru_to_buffer.py" ', '"extpar_aot_to_buffer.py" ',
            '"extpar_flake_to_buffer.exe" ', '"extpar_soil_to_buffer.exe" ',
            '"extpar_alb_to_buffer.py" ', '"extpar_ndvi_to_buffer.py" ',
            '"extpar_era_to_buffer.py" ', '"extpar_emiss_to_buffer.py" ',
            '"extpar_cdnc_to_buffer.py" ', '"extpar_edgar_to_buffer.py" ',
            '"extpar_consistency_check.exe" '
        ]
    }
    assert setup_runscript(args) == expected_runscript


def test_setup_runscript_with_art_icon():
    args = {
        'account': 'test_account',
        'lurban': False,
        'igrid_type': 1,
        'enable_cdnc': False,
        'enable_art': True,
        'enable_edgar': False
    }
    expected_runscript = {
        'account':
        'test_account',
        'pythonpath':
        os.path.join(os.getcwd(), 'lib'),
        'extpar_executables': [
            '"extpar_landuse_to_buffer.exe" ', '"extpar_topo_to_buffer.exe" ',
            '"extpar_cru_to_buffer.py" ', '"extpar_aot_to_buffer.py" ',
            '"extpar_flake_to_buffer.exe" ', '"extpar_soil_to_buffer.exe" ',
            '"extpar_alb_to_buffer.py" ', '"extpar_ndvi_to_buffer.py" ',
            '"extpar_era_to_buffer.py" ', '"extpar_emiss_to_buffer.py" ',
            '"extpar_art_to_buffer.py" ', '"extpar_consistency_check.exe" '
        ]
    }
    assert setup_runscript(args) == expected_runscript


def test_compute_aster_tiles_1_1():
    lonmax = 30.0
    lonmin = 29.0
    latmax = 30.0
    latmin = 31.0
    lsgsl = False

    expected_namelist = {
        'ntiles_column': 1,
        'ntiles_row': 1,
        'topo_files': ["'ASTER_orig_T055.nc' "],
    }
    assert compute_aster_tiles(lonmax, lonmin, latmax, latmin,
                               lsgsl) == expected_namelist


def test_compute_aster_tiles_2_1():
    lonmax = 40.0
    lonmin = 29.0
    latmax = 30.0
    latmin = 31.0
    lsgsl = False

    expected_namelist = {
        'ntiles_column': 2,
        'ntiles_row': 1,
        'topo_files': ["'ASTER_orig_T055.nc' ", "'ASTER_orig_T056.nc' "],
    }
    assert compute_aster_tiles(lonmax, lonmin, latmax, latmin,
                               lsgsl) == expected_namelist


def test_compute_aster_tiles_2_1_lsgsl():
    lonmax = 40.0
    lonmin = 29.0
    latmax = 30.0
    latmin = 31.0
    lsgsl = True

    expected_namelist = {
        'ntiles_column': 2,
        'ntiles_row': 1,
        'topo_files': ["'ASTER_orig_T055.nc' ", "'ASTER_orig_T056.nc' "],
        'sgsl_files': ["'S_ORO_T055.nc' ", "'S_ORO_T056.nc' "],
        'lpreproc_oro': '.FALSE.',
    }
    assert compute_aster_tiles(lonmax, lonmin, latmax, latmin,
                               lsgsl) == expected_namelist


def test_compute_aster_tiles_2_2():
    lonmax = 40.0
    lonmin = 29.0
    latmax = 40.0
    latmin = 31.0
    lsgsl = False

    expected_namelist = {
        'ntiles_column':
        2,
        'ntiles_row':
        2,
        'topo_files': [
            "'ASTER_orig_T043.nc' ", "'ASTER_orig_T044.nc' ",
            "'ASTER_orig_T055.nc' ", "'ASTER_orig_T056.nc' "
        ],
    }
    assert compute_aster_tiles(lonmax, lonmin, latmax, latmin,
                               lsgsl) == expected_namelist


def test_compute_aster_tiles_12_20():
    lonmax = 180.0
    lonmin = -180.0
    latmax = 60.0
    latmin = -60.0
    lsgsl = False

    expected_namelist = {
        'ntiles_column':
        12,
        'ntiles_row':
        20,
        'topo_files':
        [f"'ASTER_orig_T{number:03d}.nc' " for number in range(1, 241)]
    }
    assert compute_aster_tiles(lonmax, lonmin, latmax, latmin,
                               lsgsl) == expected_namelist


def test_compute_merit_tiles_1_1():
    lonmax = 30.0
    lonmin = 29.0
    latmax = 30.0
    latmin = 31.0
    lsgsl = False

    expected_namelist = {
        'ntiles_column': 1,
        'ntiles_row': 1,
        'topo_files': ["'MERIT_N60-N30_E000-E030.nc' "],
    }
    assert compute_merit_tiles(lonmax, lonmin, latmax, latmin,
                               lsgsl) == expected_namelist


def test_compute_merit_tiles_2_1():
    lonmax = 30.0
    lonmin = -10.0
    latmax = 30.0
    latmin = 31.0
    lsgsl = False

    expected_namelist = {
        'ntiles_column':
        2,
        'ntiles_row':
        1,
        'topo_files':
        ["'MERIT_N60-N30_W030-E000.nc' ", "'MERIT_N60-N30_E000-E030.nc' "],
    }
    assert compute_merit_tiles(lonmax, lonmin, latmax, latmin,
                               lsgsl) == expected_namelist


def test_compute_merit_tiles_2_1_sgsl():
    lonmax = 30.0
    lonmin = -10.0
    latmax = 30.0
    latmin = 31.0
    lsgsl = True

    expected_namelist = {
        'ntiles_column':
        2,
        'ntiles_row':
        1,
        'topo_files':
        ["'MERIT_N60-N30_W030-E000.nc' ", "'MERIT_N60-N30_E000-E030.nc' "],
        'sgsl_files': ["'S_ORO_0' ", "'S_ORO_1' "],
        'lpreproc_oro':
        '.TRUE.',
    }
    assert compute_merit_tiles(lonmax, lonmin, latmax, latmin,
                               lsgsl) == expected_namelist


def test_compute_merit_tiles_2_3():
    lonmax = 30.0
    lonmin = -10.0
    latmax = 60.0
    latmin = 31.0
    lsgsl = False

    expected_namelist = {
        'ntiles_column':
        2,
        'ntiles_row':
        2,
        'topo_files': [
            "'MERIT_N90-N60_W030-E000.nc' ", "'MERIT_N90-N60_E000-E030.nc' ",
            "'MERIT_N60-N30_W030-E000.nc' ", "'MERIT_N60-N30_E000-E030.nc' "
        ],
    }
    assert compute_merit_tiles(lonmax, lonmin, latmax, latmin,
                               lsgsl) == expected_namelist


def test_setup_soil_namelist_type_1():
    args = {'isoil_type': 1, 'raw_data_path': '/path/to/raw_data'}
    expected_namelist = {
        'isoil_data': 1,
        'raw_data_soil_path': '/path/to/raw_data',
        'raw_data_soil_filename': 'FAO_DSMW_double.nc',
        'soil_buffer_file': 'soil_buffer.nc',
        'path_HWSD_index_files': os.path.join('/path/to/raw_data', '../soil'),
        'lookup_table_HWSD': 'LU_TAB_HWSD_UF.data',
        'HWSD_data': 'HWSD_DATA_COSMO.data',
    }
    assert setup_soil_namelist(args) == expected_namelist


def test_setup_soil_namelist_type_2():
    args = {'isoil_type': 2, 'raw_data_path': '/path/to/raw_data'}
    expected_namelist = {
        'isoil_data': 2,
        'raw_data_soil_path': '/path/to/raw_data',
        'raw_data_soil_filename': 'HWSD0_30_topsoil.nc',
        'soil_buffer_file': 'soil_buffer.nc',
        'path_HWSD_index_files': os.path.join('/path/to/raw_data', '../soil'),
        'lookup_table_HWSD': 'LU_TAB_HWSD_UF.data',
        'HWSD_data': 'HWSD_DATA_COSMO.data',
    }
    assert setup_soil_namelist(args) == expected_namelist


def test_setup_soil_namelist_type_3():
    args = {'isoil_type': 3, 'raw_data_path': '/path/to/raw_data'}
    expected_namelist = {
        'isoil_data': 3,
        'raw_data_soil_path': '/path/to/raw_data',
        'raw_data_soil_filename': 'HWSD0_30_topsoil.nc',
        'soil_buffer_file': 'soil_buffer.nc',
        'path_HWSD_index_files': os.path.join('/path/to/raw_data', '../soil'),
        'lookup_table_HWSD': 'LU_TAB_HWSD_UF.data',
        'HWSD_data': 'HWSD_DATA_COSMO.data',
    }
    assert setup_soil_namelist(args) == expected_namelist


def test_setup_soil_namelist_unknown_type():
    args = {'isoil_type': 99, 'raw_data_path': '/path/to/raw_data'}
    with pytest.raises(ValueError, match='Unknown isoil_type 99'):
        setup_soil_namelist(args)


def test_setup_era_namelist_era5():
    args = {'iera_type': 1, 'raw_data_path': '/path/to/raw/data'}
    expected_namelist = {
        'iera_type': 1,
        'era_buffer_file': 'era_buffer.nc',
        'raw_data_era_path': '/path/to/raw/data',
        'raw_data_era_ORO': 'ERA5_ORO_1990.nc',
        'raw_data_era_SD': 'ERA5_SD_1990_2019.nc',
        'raw_data_era_T2M': 'ERA5_T2M_1990_2019.nc',
        'raw_data_era_SST': 'ERA5_SST_1990_2019.nc'
    }
    assert setup_era_namelist(args) == expected_namelist


def test_setup_era_namelist_erai():
    args = {'iera_type': 2, 'raw_data_path': '/path/to/raw/data'}
    expected_namelist = {
        'iera_type': 2,
        'era_buffer_file': 'era_buffer.nc',
        'raw_data_era_path': '/path/to/raw/data',
        'raw_data_era_ORO': 'ERA-I_ORO_1986.nc',
        'raw_data_era_SD': 'ERA-I_SD_1986_2015.nc',
        'raw_data_era_T2M': 'ERA-I_T2M_1986_2015.nc',
        'raw_data_era_SST': 'ERA-I_SST_1986_2015.nc'
    }
    assert setup_era_namelist(args) == expected_namelist


def test_setup_era_namelist_invalid():
    args = {'iera_type': 99, 'raw_data_path': '/path/to/raw/data'}
    with pytest.raises(ValueError, match='Unknown iera_type 99'):
        setup_era_namelist(args)


def test_setup_emiss_namelist_type_1():
    args = {'iemiss_type': 1, 'raw_data_path': '/path/to/raw/data'}
    expected_namelist = {
        'iemiss_type': 1,
        'emiss_buffer_file': 'emiss_buffer.nc',
        'raw_data_emiss_path': '/path/to/raw/data',
        'raw_data_emiss_filename': 'CAMEL_bbe_full_2010-2015.nc'
    }
    assert setup_emiss_namelist(args) == expected_namelist


def test_setup_emiss_namelist_type_2():
    args = {'iemiss_type': 2, 'raw_data_path': '/path/to/raw/data'}
    expected_namelist = {
        'iemiss_type': 2,
        'emiss_buffer_file': 'emiss_buffer.nc',
        'raw_data_emiss_path': '/path/to/raw/data',
        'raw_data_emiss_filename': 'CAMEL_bbe_lw_2010-2015.nc'
    }
    assert setup_emiss_namelist(args) == expected_namelist


def test_setup_emiss_namelist_invalid():
    args = {'iemiss_type': 99, 'raw_data_path': '/path/to/raw/data'}
    with pytest.raises(ValueError, match='Unknown iemiss_type 99'):
        setup_emiss_namelist(args)


def test_setup_oro_namelist_icon_globe():
    args = {
        'itopo_type': 1,
        'raw_data_path': '/path/to/raw/data',
        'lradtopo': False,
        'nhori': 24,
        'radtopo_radius': 40000.0
    }
    lonmax, lonmin, latmax, latmin = 180.0, -180.0, 90.0, -90.0
    expected_namelist = {
        'orography_buffer_file':
        'oro_buffer.nc',
        'orography_output_file':
        'oro_grid.nc',
        'lcompute_sgsl':
        ".FALSE.",
        'sgsl_buffer_file':
        'placeholder_file',
        'sgsl_files':
        'placeholder_file',
        'itopo_type':
        1,
        'raw_data_orography_path':
        '/path/to/raw/data',
        'topo_files': [
            "'GLOBE_A10.nc' ", "'GLOBE_B10.nc' ", "'GLOBE_C10.nc' ",
            "'GLOBE_D10.nc' ", "'GLOBE_E10.nc' ", "'GLOBE_F10.nc' ",
            "'GLOBE_G10.nc' ", "'GLOBE_H10.nc' ", "'GLOBE_I10.nc' ",
            "'GLOBE_J10.nc' ", "'GLOBE_K10.nc' ", "'GLOBE_L10.nc' ",
            "'GLOBE_M10.nc' ", "'GLOBE_N10.nc' ", "'GLOBE_O10.nc' ",
            "'GLOBE_P10.nc' "
        ],
        'ntiles_column':
        4,
        'ntiles_row':
        4,
        'lscale_separation':
        ".FALSE.",
        'lsso_param':
        ".TRUE.",
        'scale_sep_files':
        "'placeholder_file'",
        'raw_data_scale_sep_path':
        '/path/to/raw/data',
        'raw_data_sgsl_path':
        '/path/to/raw/data',
        'lfilter_oro':
        ".FALSE.",
        'lpreproc_oro':
        '.FALSE.',
        'ilow_pass_oro':
        4,
        'numfilt_oro':
        1,
        'ilow_pass_xso':
        5,
        'lxso_first':
        ".FALSE.",
        'numfilt_xso':
        1,
        'rxso_mask':
        750.0,
        'eps_filter':
        0.1,
        'idem_type':
        1,
        'rfill_valley':
        0.0,
        'ifill_valley':
        1,
        'lradtopo':
        ".FALSE.",
        'itype_scaling':
        0,
        'max_missing':
        0.95,
        'min_circ_cov':
        1,
        'nhori':
        24,
        'radius':
        40000.0,
    }
    assert setup_oro_namelist_icon(args, lonmax, lonmin, latmax,
                                   latmin) == expected_namelist


def test_setup_oro_namelist_icon_merit_lradtopo():
    args = {
        'itopo_type': 3,
        'raw_data_path': '/path/to/raw/data',
        'lradtopo': True,
        'nhori': 24,
        'radtopo_radius': 60000.0
    }
    lonmax = 30.0
    lonmin = -10.0
    latmax = 60.0
    latmin = 31.0
    expected_namelist = {
        'orography_buffer_file':
        'oro_buffer.nc',
        'orography_output_file':
        'oro_grid.nc',
        'lcompute_sgsl':
        ".FALSE.",
        'sgsl_buffer_file':
        'placeholder_file',
        'sgsl_files':
        'placeholder_file',
        'itopo_type':
        3,
        'raw_data_orography_path':
        '/path/to/raw/data',
        'topo_files': [
            "'MERIT_N90-N60_W030-E000.nc' ", "'MERIT_N90-N60_E000-E030.nc' ",
            "'MERIT_N60-N30_W030-E000.nc' ", "'MERIT_N60-N30_E000-E030.nc' "
        ],
        'ntiles_column':
        2,
        'ntiles_row':
        2,
        'lscale_separation':
        ".FALSE.",
        'lsso_param':
        ".TRUE.",
        'scale_sep_files':
        "'placeholder_file'",
        'raw_data_scale_sep_path':
        '/path/to/raw/data',
        'raw_data_sgsl_path':
        '/path/to/raw/data',
        'lfilter_oro':
        ".FALSE.",
        'lpreproc_oro':
        '.FALSE.',
        'ilow_pass_oro':
        4,
        'numfilt_oro':
        1,
        'ilow_pass_xso':
        5,
        'lxso_first':
        ".FALSE.",
        'numfilt_xso':
        1,
        'rxso_mask':
        750.0,
        'eps_filter':
        0.1,
        'idem_type':
        3,
        'rfill_valley':
        0.0,
        'ifill_valley':
        1,
        'lradtopo':
        ".TRUE.",
        'radius':
        60000.0,
        'nhori':
        24,
        'max_missing':
        0.95,
        'min_circ_cov':
        1,
        'itype_scaling':
        0,
    }
    assert setup_oro_namelist_icon(args, lonmax, lonmin, latmax,
                                   latmin) == expected_namelist


def test_setup_oro_namelist_icon_invalid():
    args = {
        'itopo_type': 4,
        'raw_data_path': '/path/to/raw/data',
        'lradtopo': False,
        'radtopo_radius': 40000.0
    }
    lonmax, lonmin, latmax, latmin = 180.0, -180.0, 90.0, -90.0
    with pytest.raises(ValueError, match='Unknown itopo_type 4'):
        setup_oro_namelist_icon(args, lonmax, lonmin, latmax, latmin)


def test_setup_cdnc_namelist():
    args = {'raw_data_path': '/path/to/raw/data'}
    expected_namelist = {
        'raw_data_cdnc_path': '/path/to/raw/data',
        'cdnc_buffer_file': 'cdnc_buffer.nc',
        'raw_data_cdnc_filename': 'modis_cdnc_climatology_Q06.nc'
    }
    assert setup_cdnc_namelist(args) == expected_namelist


def test_setup_edgar_namelist():
    args = {'raw_data_path': '/path/to/raw/data'}
    expected_namelist = {
        'raw_data_edgar_path': '/path/to/raw/data',
        'raw_data_edgar_filename_bc': 'v8.1_FT2022_AP_BC_2022_TOTALS_flx.nc',
        'raw_data_edgar_filename_oc': 'v8.1_FT2022_AP_OC_2022_TOTALS_flx.nc',
        'raw_data_edgar_filename_so2': 'v8.1_FT2022_AP_SO2_2022_TOTALS_flx.nc',
        'raw_data_edgar_filename_nox': 'v8.1_FT2022_AP_NOx_2022_TOTALS_flx.nc',
        'raw_data_edgar_filename_nh3': 'v8.1_FT2022_AP_NH3_2022_TOTALS_flx.nc',
        'edgar_buffer_file': 'edgar_buffer.nc'
    }
    assert setup_edgar_namelist(args) == expected_namelist


def test_all_placeholders_replaced_cosmo(tmp_dir):
    args = {
        "igrid_type": 2,
        "input_grid": 'test/testsuite/data/clm/12km_globe/INPUT_COSMO_GRID',
        "iaot_type": 1,
        "ilu_type": 1,
        "ialb_type": 1,
        "isoil_type": 1,
        "itopo_type": 1,
        "it_cl_type": 1,
        "iera_type": 1,
        "iemiss_type": 1,
        'ilookup_table_lu': 1,
        "enable_cdnc": False,
        "enable_edgar": False,
        "enable_art": False,
        'tile_mode': 0,
        "use_array_cache": False,
        "nhori": 24,
        "radtopo_radius": 40000.0,
        'tcorr_lapse_rate': 0.0065,
        'tcorr_offset': 0.0,
        "raw_data_path": '/dummy/raw_data_path',
        "run_dir": str(tmp_dir),
        "account": 'dummy_account',
        "host": 'docker',
        "no_batch_job": True,
        "lurban": False,
        "lsgsl": False,
        "lfilter_oro": False,
        "l_use_corine": False,
        "lradtopo": False
    }

    namelist = setup_namelist(args)
    runscript = setup_runscript(args)
    # in this function all placeholders from the template are replaced,
    # it raises an error if there are any placeholders left
    prepare_sandbox(args, namelist, runscript, test_run=True)


def test_all_placeholders_replaced_icon(tmp_dir, icon_grid):
    args = {
        "igrid_type": 1,
        "input_grid": icon_grid,
        "iaot_type": 1,
        "ilu_type": 1,
        "ialb_type": 1,
        "isoil_type": 1,
        "itopo_type": 1,
        "it_cl_type": 1,
        "iera_type": 1,
        "iemiss_type": 1,
        'ilookup_table_lu': 1,
        "enable_cdnc": False,
        "enable_edgar": False,
        "enable_art": False,
        'tile_mode': 1,
        "use_array_cache": False,
        "nhori": 24,
        "radtopo_radius": 40000.0,
        'tcorr_lapse_rate': 0.0065,
        'tcorr_offset': 0.0,
        "raw_data_path": '/dummy/raw_data_path',
        "run_dir": str(tmp_dir),
        "account": 'dummy_account',
        "host": 'docker',
        "no_batch_job": True,
        "lurban": False,
        "lsgsl": False,
        "lfilter_oro": False,
        "l_use_corine": False,
        "lradtopo": False
    }

    namelist = setup_namelist(args)
    runscript = setup_runscript(args)
    # in this function all placeholders from the template are replaced,
    # it raises an error if there are any placeholders left
    prepare_sandbox(args, namelist, runscript, test_run=True)
