'''
This is a template for the namelist -> rename to namelist.py for use

Make sure, the namelist.py is contained in the run-directory of Extpar.
All Python scripts use the same namelist.py file, but use different
dictionaries:

    -extpar_alb_to_buffer.py:     input_alb

    -extpar_cru_to_buffer.py:     input_tclim

    -extpar_emiss_to_buffer.py:   input_emiss

    -extpar_ndvi_to_buffer.py:    input_ndvi

    -extpar_era_to_buffer.py:     input_era

    -extpar_ahf_to_buffer.py:     input_ahf

    -extpar_isa_to_buffer.py:     input_isa  
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

input_edgar = {
    'raw_data_edgar_path': '',
    'raw_data_edgar_filename_bc': 'EDGARv6.1_BC_2018_TOTALS.0.1x0.1.nc',
    'raw_data_edgar_filename_oc': 'EDGARv6.1_OC_2018_TOTALS.0.1x0.1.nc',
    'raw_data_edgar_filename_so2': 'EDGARv6.1_SO2_2018_TOTALS.0.1x0.1.nc',
    'edgar_buffer_file': 'edgar_buffer.nc',
}

input_era = {
    'iera_type': 1,
    'raw_data_era_path': '',
    'raw_data_era_ORO': 'ERA5_ORO_1990.nc',
    'raw_data_era_T2M': 'ERA5_T2M_1990_2019.nc',
    'raw_data_era_SST': 'ERA5_SST_1990_2019.nc.',
    'raw_data_era_SD': 'ERA5_SD_1990_2019.nc',
    'era_buffer_file': 'era_buffer.nc',
}

input_ahf = {
    'iahf_type': 1,
    'raw_data_ahf_path': '',
    'raw_data_ahf_filename': 'AHF_2006_CDO.nc',
    'ahf_buffer_file': 'ahf_buffer.nc',
}

input_isa = {
    'raw_data_isa_path': '',
    'raw_data_isa_filename': 'NOAA_ISA_CDO.nc',
    'isa_buffer_file': 'isa_buffer.nc',
}
