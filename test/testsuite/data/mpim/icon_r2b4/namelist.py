input_tclim = {
    'raw_data_t_clim_path': '',
    'raw_data_tclim_coarse': 'absolute_hadcrut3.nc',
    'raw_data_tclim_fine': 'CRU_T_SOIL_clim.nc',
    't_clim_buffer_file': 'cru_buffer.nc',
    'it_cl_type': 2
}

input_alb = {
    'ialb_type': 1,
    'raw_data_alb_path': '',
    'raw_data_alb_filename': 'alb_new.nc',
    'raw_data_alnid_filename': 'alnid_new.nc',
    'raw_data_aluvd_filename': 'aluvd_new.nc',
    'alb_buffer_file': 'albedo_buffer.nc'
}

input_emiss = {
    'iemiss_type': 2,
    'raw_data_emiss_path': '',
    'raw_data_emiss_filename': 'CAMEL_bbe_lw_2010-2015.nc',
    'emiss_buffer_file': 'emiss_buffer.nc'
}

input_ndvi = {
    'raw_data_ndvi_path': '',
    'raw_data_ndvi_filename': 'NDVI_1998_2003.nc',
    'ndvi_buffer_file': 'ndvi_buffer.nc'
}

input_edgar = {
    'raw_data_edgar_path': '',
    'raw_data_edgar_filename_bc': 'v8.1_FT2022_AP_BC_2022_TOTALS_flx.nc',
    'raw_data_edgar_filename_oc': 'v8.1_FT2022_AP_OC_2022_TOTALS_flx.nc',
    'raw_data_edgar_filename_so2': 'v8.1_FT2022_AP_SO2_2022_TOTALS_flx.nc',
    'raw_data_edgar_filename_nox': 'v8.1_FT2022_AP_NOx_2022_TOTALS_flx.nc',
    'raw_data_edgar_filename_nh3': 'v8.1_FT2022_AP_NH3_2022_TOTALS_flx.nc',
    'edgar_buffer_file': 'edgar_buffer.nc',
}
