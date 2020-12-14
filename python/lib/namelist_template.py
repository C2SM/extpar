'''
This is a template for the namelist -> rename to namelist.py for use

Make sure, the namelist.py is contained in the run-directory of Extpar.
All Python scripts use the same namelist.py file, but use different
dictionaries:

    -extpar_alb_to_buffer.py:     input_alb

    -extpar_cru_to_buffer.py:     input_tclim

    -extpar_emiss_to_buffer.py:   input_emiss

    -extpar_ndvi_to_buffer.py:    input_ndvi
'''

input_tclim = {
        'raw_data_t_clim_path': '',
        'raw_data_tclim_coarse': 'absolute_hadcrut3.nc',
        'raw_data_tclim_fine': 'CRU_T2M_SURF_clim.nc',
        't_clim_buffer_file': 'cru_buffer.nc',
        'it_cl_type': 1
        }

input_alb = {
        'ialb_type': 2,
        'raw_data_alb_path': '',
        'raw_data_alb_filename': 'global_soil_albedo.nc',
        'raw_data_alnid_filename': 'month_alnid_new.nc',
        'raw_data_aluvd_filename': 'month_aluvd_new.nc',
        'alb_buffer_file': 'albedo_buffer.nc',
        'alb_output_file': 'albedo_cosmo.nc',
        }

input_emiss = { 
        'iemiss_type': 1, 
        'raw_data_emiss_path': '',
        'raw_data_emiss_filename': 'CAM_bbe_int_2010-2015_full.nc',
        'emiss_buffer_file': 'emiss_buffer.nc',
        'emiss_output_file': 'emiss_icon.nc'
        }

input_ndvi = {
        'raw_data_ndvi_path': '',
        'raw_data_ndvi_filename': 'NDVI_1998_2003.nc',
        'ndvi_buffer_file': 'ndvi_buffer.nc',
        'ndvi_output_file': 'ndvi_extpar_cosmo.nc'
        }

input_era = {
        'raw_data_era_path': '',
        'raw_data_era_ORO': 'ERA5_ORO_1990.nc',
        'raw_data_era_T2M': 'ERA5_T2M_1990_2019.nc',
        'raw_data_era_SST': 'ERA5_SST_1990_2019.nc.',
        'raw_data_era_SD': 'ERA5_SD_1990_2019.nc',
        'era_buffer_file': 'era_buffer.nc',
        }
