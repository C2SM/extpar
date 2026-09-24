input_tclim = {
        'raw_data_t_clim_path': '',
        'raw_data_tclim_coarse': 'absolute_hadcrut3.nc',
        'raw_data_tclim_fine': 'CRU_T_SOIL_clim.nc',
        't_clim_buffer_file': 'tclim_buffer.nc',
        'it_cl_type': 2,
        'tcorr_lapse_rate': 0.006,
        'tcorr_offset': 1
        }

input_alb = {
        'ialb_type': 1,
        'raw_data_alb_path': '',
        'raw_data_alb_filename': 'alb_new.nc',
        'raw_data_alnid_filename': 'alnid_new.nc',
        'raw_data_aluvd_filename': 'aluvd_new.nc',
        'alb_buffer_file': 'alb_buffer.nc'
        }

input_ndvi = {
        'raw_data_ndvi_path': '',
        'raw_data_ndvi_filename': 'NDVI_1998_2003.nc',
        'ndvi_buffer_file': 'ndvi_buffer.nc'
        }

input_ahf = {
        'iahf_type': 1,
        'raw_data_ahf_path': '',
        'raw_data_ahf_filename': 'AHF_2006_2.5min_lonlat.nc',
        'ahf_buffer_file': 'ahf_buffer.nc'
        }

input_isa = {
        'isa_type': 1,
        'raw_data_isa_path': '',
        'raw_data_isa_filename': 'NOAA_ISA_16bit_lonlat.nc',
        'isa_buffer_file': 'isa_buffer.nc'
        }

input_era = {
    'iera_type': 1,
    'raw_data_era_path': '',
    'raw_data_era_ORO': 'ERA5_ORO_1990.nc',
    'raw_data_era_T2M': 'ERA5_T2M_1990_2019.nc',
    'raw_data_era_SST': 'ERA5_SST_1990_2019.nc',
    'raw_data_era_SD': 'ERA5_SD_1990_2019.nc',
    'era_buffer_file': 'era_buffer.nc',
}

input_emiss = {
    'iemiss_type': 2,
    'raw_data_emiss_path': '',
    'raw_data_emiss_filename': 'CAMEL_bbe_lw_2010-2015.nc',
    'emiss_buffer_file': 'emiss_buffer.nc',

}
input_cdnc = {
    'icdnc_type': 1,
    'raw_data_cdnc_path': '',
    'raw_data_cdnc_filename': 'modis_cdnc_climatology_Q06.nc',
    'cdnc_buffer_file': 'cdnc_buffer.nc',
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

input_gfasclim = {
    'raw_data_gfasclim_path': '',
    'raw_data_gfasclim_filename': 'gfasclim2015-2024.nc',
    'gfasclim_buffer_file': 'gfasclim_buffer.nc',
}

input_aot = {
    'iaot_type': 1,
    'raw_data_aot_path': '',
    'raw_data_aot_filename': 'aot_GACP_sea_salt_fixed.nc',
    'aot_buffer_file': 'aot_buffer.nc'
}

input_art = {
    'raw_data_art_path': '',
    'raw_data_art_filename': 'HWSD0_USDA.nc',
    'art_buffer_file': 'art_buffer.nc',
}

input_oro = {
    'itopo_type': 3,
    'orography_buffer_file': 'oro_buffer.nc',
    'raw_data_path': '',
}

input_radtopo = {
    'lradtopo': True,
    'itype_scaling': 0,
    'nhori': 24,
    'radius': 40000,
    'radtopo_type': 3,
    'radtopo_buffer_file': 'radtopo_buffer.nc',
    'ray_origin_elev': 0.2,
    'num_nodes': 7,
}
