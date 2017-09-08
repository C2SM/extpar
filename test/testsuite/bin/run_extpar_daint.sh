#!/usr/bin/ksh

ulimit -s unlimited
ulimit -c unlimited


# Variables which are automatically set
currentdir=`pwd`
rootdir=${currentdir}/..
progdir=${rootdir}/bin
scriptpath=$0
scriptname=${scriptpath##*/}
logfile=${scriptname%.*}_`date +%Y%m%d%H%M%S`.log

#---------------------------------------------------------------------------------------------------------
# NetCDF raw data for external parameter; adjust the path setting!
data_dir=/store/s83/tsm/extpar/raw_data_nc/

# Output file format and names; adjust!
grib_sample='rotated_ll_pl_grib1'
grib_output_filename='external_parameter_mch_cosmo7.g1'
netcdf_output_filename='external_parameter_mch_cosmo7.nc'
#---------------------------------------------------------------------------------------------------------


# Names of executables
binary_alb=extpar_alb_to_buffer.exe
binary_lu=extpar_landuse_to_buffer.exe
binary_topo=extpar_topo_to_buffer.exe
binary_aot=extpar_aot_to_buffer.exe
binary_tclim=extpar_cru_to_buffer.exe
binary_ndvi=extpar_ndvi_to_buffer.exe
binary_soil=extpar_soil_to_buffer.exe
binary_flake=extpar_flake_to_buffer.exe

binary_consistency_check=extpar_consistency_check.exe

if [[ -e ${logfile} ]] ; then
  rm -f ${logfile}
fi
echo "\n>>>> Data will be processed and produced in `pwd` <<<<\n"

#---

raw_data_alb='MODIS_month_alb.nc'
raw_data_alnid='MODIS_month_alnid.nc'
raw_data_aluvd='MODIS_month_aluvd.nc'
buffer_alb='month_alb_buffer.nc'
output_alb='month_alb_extpar_cosmo.nc'

raw_data_aot='aerosol_optical_thickness.nc'
buffer_aot='extpar_buffer_aot.nc'
output_aot='aot_extpar_cosmo.nc'

raw_data_tclim_coarse='CRU_T2M_SURF_clim_coarse.nc'
raw_data_tclim_fine='CRU_T2M_SURF_clim_fine.nc'
buffer_tclim='crutemp_clim_extpar_buffer.nc'
output_tclim='crutemp_clim_extpar_cosmo.nc'

raw_data_glc2000='GLC2000_byte.nc'
buffer_glc2000='extpar_landuse_buffer.nc'
output_glc2000='extpar_landuse_cosmo.nc'
raw_data_glcc='GLCC_usgs_class_byte.nc'
buffer_glcc='glcc_landuse_buffer.nc'
output_glcc='glcc_landuse_cosmo.nc'

# raw_data_globcover='GLOBCOVER_L4_200901_200912_V2.3_int16.nc'
raw_data_globcover_0='GLOBCOVER_0_16bit.nc'
raw_data_globcover_1='GLOBCOVER_1_16bit.nc'
raw_data_globcover_2='GLOBCOVER_2_16bit.nc'
raw_data_globcover_3='GLOBCOVER_3_16bit.nc'
raw_data_globcover_4='GLOBCOVER_4_16bit.nc'
raw_data_globcover_5='GLOBCOVER_5_16bit.nc'

buffer_lu='extpar_landuse_buffer.nc'
output_lu='extpar_landuse_cosmo.nc'

raw_data_globe_A10='GLOBE_A10.nc'
raw_data_globe_B10='GLOBE_B10.nc'
raw_data_globe_C10='GLOBE_C10.nc'
raw_data_globe_D10='GLOBE_D10.nc'
raw_data_globe_E10='GLOBE_E10.nc'
raw_data_globe_F10='GLOBE_F10.nc'
raw_data_globe_G10='GLOBE_G10.nc'
raw_data_globe_H10='GLOBE_H10.nc'
raw_data_globe_I10='GLOBE_I10.nc'
raw_data_globe_J10='GLOBE_J10.nc'
raw_data_globe_K10='GLOBE_K10.nc'
raw_data_globe_L10='GLOBE_L10.nc'
raw_data_globe_M10='GLOBE_M10.nc'
raw_data_globe_N10='GLOBE_N10.nc'
raw_data_globe_O10='GLOBE_O10.nc'
raw_data_globe_P10='GLOBE_P10.nc'

# lanczos filter is recommended when activating scale separation
raw_filt_globe_A10='GLOBE_A_filt_lanczos_window.nc'
raw_filt_globe_B10='GLOBE_B_filt_lanczos_window.nc'
raw_filt_globe_C10='GLOBE_C_filt_lanczos_window.nc'
raw_filt_globe_D10='GLOBE_D_filt_lanczos_window.nc'
raw_filt_globe_E10='GLOBE_E_filt_lanczos_window.nc'
raw_filt_globe_F10='GLOBE_F_filt_lanczos_window.nc'
raw_filt_globe_G10='GLOBE_G_filt_lanczos_window.nc'
raw_filt_globe_H10='GLOBE_H_filt_lanczos_window.nc'
raw_filt_globe_I10='GLOBE_I_filt_lanczos_window.nc'
raw_filt_globe_J10='GLOBE_J_filt_lanczos_window.nc'
raw_filt_globe_K10='GLOBE_K_filt_lanczos_window.nc'
raw_filt_globe_L10='GLOBE_L_filt_lanczos_window.nc'
raw_filt_globe_M10='GLOBE_M_filt_lanczos_window.nc'
raw_filt_globe_N10='GLOBE_N_filt_lanczos_window.nc'
raw_filt_globe_O10='GLOBE_O_filt_lanczos_window.nc'
raw_filt_globe_P10='GLOBE_P_filt_lanczos_window.nc'

raw_data_aster_T01='ASTER_eu_T01.nc'
raw_data_aster_T02='ASTER_eu_T02.nc'
raw_data_aster_T03='ASTER_eu_T03.nc'
raw_data_aster_T04='ASTER_eu_T04.nc'
raw_data_aster_T05='ASTER_eu_T05.nc'
raw_data_aster_T06='ASTER_eu_T06.nc'
raw_data_aster_T07='ASTER_eu_T07.nc'
raw_data_aster_T08='ASTER_eu_T08.nc'
raw_data_aster_T09='ASTER_eu_T09.nc'
raw_data_aster_T10='ASTER_eu_T10.nc'
raw_data_aster_T11='ASTER_eu_T11.nc'
raw_data_aster_T12='ASTER_eu_T12.nc'
raw_data_aster_T13='ASTER_eu_T13.nc'
raw_data_aster_T14='ASTER_eu_T14.nc'
raw_data_aster_T15='ASTER_eu_T15.nc'
raw_data_aster_T16='ASTER_eu_T16.nc'
raw_data_aster_T17='ASTER_eu_T17.nc'
raw_data_aster_T18='ASTER_eu_T18.nc'
raw_data_aster_T19='ASTER_eu_T19.nc'
raw_data_aster_T20='ASTER_eu_T20.nc'
raw_data_aster_T21='ASTER_eu_T21.nc'
raw_data_aster_T22='ASTER_eu_T22.nc'
raw_data_aster_T23='ASTER_eu_T23.nc'
raw_data_aster_T24='ASTER_eu_T24.nc'
raw_data_aster_T25='ASTER_eu_T25.nc'
raw_data_aster_T26='ASTER_eu_T26.nc'
raw_data_aster_T27='ASTER_eu_T27.nc'
raw_data_aster_T28='ASTER_eu_T28.nc'
raw_data_aster_T29='ASTER_eu_T29.nc'
raw_data_aster_T30='ASTER_eu_T30.nc'
raw_data_aster_T31='ASTER_eu_T31.nc'
raw_data_aster_T32='ASTER_eu_T32.nc'
raw_data_aster_T33='ASTER_eu_T33.nc'
raw_data_aster_T34='ASTER_eu_T34.nc'
raw_data_aster_T35='ASTER_eu_T35.nc'
raw_data_aster_T36='ASTER_eu_T36.nc'

buffer_topo='topography_buffer.nc'
output_topo='topography_COSMO.nc'

raw_data_ndvi='NDVI_1998_2003.nc'
buffer_ndvi='ndvi_buffer.nc'
output_ndvi='ndvi_extpar_cosmo.nc'

raw_data_soil_FAO='FAO_DSMW_DP.nc'
raw_data_soil_HWSD='HWSD0_30_texture_2.nc'
raw_data_deep_soil='HWSD30_100_texture_2.nc'
buffer_soil='soil_buffer.nc'
output_soil='soil_COSMO.nc'

raw_lookup_table_HWSD='LU_TAB_HWSD_UF.data'
raw_HWSD_data='HWSD_DATA_COSMO.data'
raw_HWSD_data_deep='HWSD_DATA_COSMO_S.data'
raw_HWSD_data_extpar='HWSD_DATA_COSMO_EXTPAR.asc'

raw_data_flake='GLDB_lakedepth.nc'
buffer_flake='flake_buffer.nc'
output_flake='ext_par_flake_cosmo.nc'

# link raw data files to local workdir
ln -s -f ${data_dir}/*.nc .

# run the programs
# the next seven programs can run independent of each other
echo "\n>> Run ${binary_alb} ..."  ;  time ./${binary_alb} 2>&1 >> ${logfile}
echo "\n>> Run ${binary_aot} ..."  ;  time ./${binary_aot} 2>&1 >> ${logfile}
echo "\n>> Run ${binary_tclim} ..."  ;  time ./${binary_tclim} 2>&1 >> ${logfile}
echo "\n>> Run ${binary_lu} ..."  ;  time ./${binary_lu} 2>&1 >> ${logfile}
echo "\n>> Run ${binary_topo} ..."  ;  time ./${binary_topo} 2>&1 >> ${logfile}
echo "\n>> Run ${binary_ndvi} ..."  ;  time ./${binary_ndvi} 2>&1 >> ${logfile}
echo "\n>> Run ${binary_soil} ..."  ;  time ./${binary_soil} 2>&1 >> ${logfile}
echo "\n>> Run ${binary_flake} ..."  ;  time ./${binary_flake} 2>&1 >> ${logfile}

# the consistency check requires the output of 
# ${binary_aot}, ${binary_tclim}, ${binary_lu}, ${binary_globe}, 
# ${binary_ndvi}, ${binary_soil} and ${binary_flake}
echo "\n>> Run ${binary_consistency_check} ..."  ;  time ./${binary_consistency_check} 2>&1 >> ${logfile}

ls
echo "\n>>>> External parameters for COSMO model generated <<<<\n"




