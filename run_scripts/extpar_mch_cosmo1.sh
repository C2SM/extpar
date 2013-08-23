#!/usr/bin/ksh

ulimit -s unlimited
ulimit -c unlimited


# GRIB API resources; adjust the path setting!
export GRIB_DEFINITION_PATH=/oprusers/osm/lib/libgrib_api_1.9.9.1_pgi12.2.0/share/definitions:/oprusers/osm/lib/libgrib_api_1.9.9.1_pgi12.2.0/share/definitions

export GRIB_SAMPLES_PATH=/oprusers/osm/lib/libgrib_api_1.9.9.1_pgi12.2.0/share/samples

# NetCDF raw data for external parameter; adjust the path setting!
data_dir=/store/s83/tsm/extpar/raw_data_nc/

# Sandbox; adjust the path setting (make sure you have enough disk place at that location)!
sandboxdir=/store/s83/tsm/extpar/sandbox


# Names of executables
binary_alb=extpar_alb_to_buffer.exe
binary_lu=extpar_landuse_to_buffer.exe
binary_topo=extpar_topo_to_buffer.exe
binary_aot=extpar_aot_to_buffer.exe
binary_tclim=extpar_cru_to_buffer.exe
binary_ndvi=extpar_ndvi_to_buffer.exe
binary_soil=extpar_soil_to_buffer.exe
binary_flake=extpar_flake_to_buffer.exe

binary_extpar_consistency_check=extpar_consistency_check.exe


currentdir=`pwd`
progdir=${currentdir}/../bin

if [[ ! -d ${sandboxdir} ]] ; then
  mkdir -p ${sandboxdir}
fi
cd ${sandboxdir}
echo "\n>>>> Data will be processed and produced in `pwd` <<<<\n"



#---

grib_output_filename='external_parameter_mch_cosmo1.g1'
stf_output_filename='external_parameter_mch_cosmo1.stf'
netcdf_output_filename='external_parameter_mch_cosmo1.nc'
grib_sample='DWD_rotated_ll_7km_G_grib1'

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


#raw_data_globcover_gz='GLOBCOVER_L4_200901_200912_V2.3_int16.nc.gz'
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
buffer_ndvi='NDVI_buffer.nc'
output_ndvi='ndvi_extpar_cosmo.nc'

raw_data_soil_FAO='FAO_DSMW_DP.nc'
raw_data_soil_HWSD='HWSD0_30_texture_2.nc'
raw_data_deep_soil='HWSD30_100_texture_2.nc'
buffer_soil='FAO_DSMW_buffer.nc'
output_soil='FAO_DSMW_COSMO.nc'

raw_lookup_table_HWSD='LU_TAB_HWSD_UF.data'
raw_HWSD_data='HWSD_DATA_COSMO.data'
raw_HWSD_data_deep='HWSD_DATA_COSMO_S.data'
raw_HWSD_data_extpar='HWSD_DATA_COSMO_EXTPAR.asc'

raw_data_flake='GLDB_lakedepth.nc'
buffer_flake='flake_buffer.nc'
output_flake='ext_par_flake_cosmo.nc'

# set target grid definition 
cat > INPUT_grid_org << EOF_go
&GRID_DEF 
 igrid_type = 2,
 domain_def_namelist='INPUT_COSMO_GRID'
/ 
EOF_go

#---
cat > INPUT_COSMO_GRID << EOF_grid
&lmgrid
 pollon=-170.0, 
 pollat=43.0, 
 startlon_tot=-10.0, 
 startlat_tot=-7.0,
 dlon=0.01,
 dlat=0.01,
 ie_tot=1701,
 je_tot=1401,
/
EOF_grid
#---
# create input namelists 
#---
cat > INPUT_ALB << EOF_alb
&alb_raw_data
  raw_data_alb_path='',
  raw_data_alb_filename='${raw_data_alb}'
/
&alnid_raw_data
  raw_data_alb_path='',
  raw_data_alnid_filename='${raw_data_alnid}'
/
&aluvd_raw_data
  raw_data_alb_path='',
  raw_data_aluvd_filename='${raw_data_aluvd}'
/
&alb_io_extpar
  alb_buffer_file='${buffer_alb}',
  alb_output_file='${output_alb}'
/
&alb_source_file
  alb_source='al',
  alnid_source='alnid',
  aluvd_source='aluvd'
/
EOF_alb
#---
cat > INPUT_AOT << EOF_aot
&aerosol_raw_data
  raw_data_aot_path='',
  raw_data_aot_filename='${raw_data_aot}'
/  
&aerosol_io_extpar
  aot_buffer_file='${buffer_aot}',
  aot_output_file='${output_aot}'
/
EOF_aot
#---
cat > INPUT_TCLIM << EOF_tclim
&t_clim_raw_data
  raw_data_t_clim_path='',
  raw_data_t_clim_filename='${raw_data_tclim_fine}',
  raw_data_t_id = 1
/  

&t_clim_io_extpar
  t_clim_buffer_file='${buffer_tclim}',
  t_clim_output_file='${output_tclim}'
/  
EOF_tclim
#---
cat > INPUT_LU << EOF_lu
&lu_raw_data
   raw_data_lu_path='',
   raw_data_lu_filename='${raw_data_globcover_0}' '${raw_data_globcover_1}' '${raw_data_globcover_2}' '${raw_data_globcover_3}' '${raw_data_globcover_4}' '${raw_data_globcover_5}',
   i_landuse_data=1,
   ilookup_table_lu=1 
/
&lu_io_extpar
   lu_buffer_file='${buffer_lu}',
   lu_output_file='${output_lu}'
/
&glcc_raw_data
   raw_data_glcc_path='',
   raw_data_glcc_filename='${raw_data_glcc}'
/
&glcc_io_extpar
   glcc_buffer_file='${buffer_glcc}',
   glcc_output_file='${output_glcc}'
/
EOF_lu
#---
cat > INPUT_ORO << EOF_oro
&orography_io_extpar
  orography_buffer_file='${buffer_topo}',
  orography_output_file='${output_topo}'
/
&orography_raw_data
 itopo_type = 2,
 lsso_param = .TRUE.,
 raw_data_orography_path='',
 ntiles_column = 2,
 ntiles_row = 4,
 topo_files = 'ASTER_orig_T006.nc' 'ASTER_orig_T007.nc' 'ASTER_orig_T018.nc' 'ASTER_orig_T019.nc' 'ASTER_orig_T030.nc' 'ASTER_orig_T031.nc' 'ASTER_orig_T042.nc' 'ASTER_orig_T043.nc' 
/
EOF_oro
#--- topo_FILES = '${raw_data_globe_A10}' '${raw_data_globe_B10}'  '${raw_data_globe_C10}'  '${raw_data_globe_D10}'  '${raw_data_globe_E10}'  '${raw_data_globe_F10}'  '${raw_data_globe_G10}'  '${raw_data_globe_H10}'  '${raw_data_globe_I10}'  '${raw_data_globe_J10}'  '${raw_data_globe_K10}'  '${raw_data_globe_L10}'  '${raw_data_globe_M10}'  '${raw_data_globe_N10}'  '${raw_data_globe_O10}'  '${raw_data_globe_P10}'topo_FILES = '${raw_data_aster_T01}' '${raw_data_aster_T02}'  '${raw_data_aster_T03}'  '${raw_data_aster_T04}'  '${raw_data_aster_T05}'  '${raw_data_aster_T06}'  '${raw_data_aster_T07}'  '${raw_data_aster_T08}'  '${raw_data_aster_T09}'  '${raw_data_aster_T10}'  '${raw_data_aster_T11}'  '${raw_data_aster_T12}'  '${raw_data_aster_T13}'  '${raw_data_aster_T14}'  '${raw_data_aster_T15}'  '${raw_data_aster_T16}'  '${raw_data_aster_T17}'  '${raw_data_aster_T18}'  '${raw_data_aster_T19}'  '${raw_data_aster_T20}' '${raw_data_aster_T21}'  '${raw_data_aster_T22}'  '${raw_data_aster_T23}'  '${raw_data_aster_T24}'  '${raw_data_aster_T25}'  '${raw_data_aster_T26}'  '${raw_data_aster_T27}'  '${raw_data_aster_T28}'  '${raw_data_aster_T29}'  '${raw_data_aster_T30}'  '${raw_data_aster_T31}'  '${raw_data_aster_T32}'  '${raw_data_aster_T33}'  '${raw_data_aster_T34}'  '${raw_data_aster_T35}'  '${raw_data_aster_T36}'  
#---
cat > INPUT_OROSMOOTH << EOF_orosm
&orography_smoothing
  lfilter_oro=.TRUE.,
  ilow_pass_oro=4,
  numfilt_oro=1,
  ilow_pass_xso=5,
  lxso_first=.FALSE.,
  numfilt_xso=1,
  rxso_mask=750.0,
  eps_filter=0.1,
  rfill_valley=0.0,
  ifill_valley=7
/
EOF_orosm
#---
cat > INPUT_RADTOPO << EOF_radtopo
&radtopo
  lradtopo=.TRUE.,
  nhori=24
/
EOF_radtopo
#---
cat > INPUT_SCALE_SEP << EOF_scale_sep
&scale_separated_raw_data
  lscale_separation = .FALSE.,
  raw_data_scale_sep_path = '',
  scale_sep_files = '${raw_data_globe_A10}' '${raw_data_globe_B10}'  '${raw_data_globe_C10}'  '${raw_data_globe_D10}'  '${raw_data_globe_E10}'  '${raw_data_globe_F10}'  '${raw_data_globe_G10}'  '${raw_data_globe_H10}'  '${raw_data_globe_I10}'  '${raw_data_globe_J10}'  '${raw_data_globe_K10}'  '${raw_data_globe_L10}'  '${raw_data_globe_M10}'  '${raw_data_globe_N10}'  '${raw_data_globe_O10}'  '${raw_data_globe_P10}'
/
EOF_scale_sep
#---
cat > INPUT_NDVI << EOF_ndvi
&ndvi_raw_data
  raw_data_ndvi_path='',
  raw_data_ndvi_filename='${raw_data_ndvi}'
/  
&ndvi_io_extpar
 ndvi_buffer_file='${buffer_ndvi}',
 ndvi_output_file='${output_ndvi}'
/
EOF_ndvi
#---
cat > INPUT_SOIL << EOF_soil
&soil_raw_data
 isoil_data = 1,
 ldeep_soil = .FALSE.,
 raw_data_soil_path='',
 raw_data_soil_filename='${raw_data_soil_FAO}',
 raw_data_deep_soil_filename='${raw_data_deep_soil}'
/
&soil_io_extpar
  soil_buffer_file='${buffer_soil}',
  soil_output_file='${output_soil}'
/
&HWSD_index_files
 path_HWSD_index_files='',
 lookup_table_HWSD='${raw_lookup_table_HWSD}', 
 HWSD_data='${raw_HWSD_data}',
 HWSD_data_deep='${raw_HWSD_data_deep}',
 HWSD_data_extpar='${raw_HWSD_data_extpar}'
/
EOF_soil
#---
cat > INPUT_FLAKE << EOF_flake
&flake_raw_data
   raw_data_flake_path='',
   raw_data_flake_filename='${raw_data_flake}'
/
&flake_io_extpar
   flake_buffer_file='${buffer_flake}'
   flake_output_file='${output_flake}'
/
EOF_flake
#---
# consistency check
cat > INPUT_CHECK << EOF_check
&extpar_consistency_check_io
  grib_output_filename="${grib_output_filename}",
  grib_sample="${grib_sample}",
  netcdf_output_filename="${netcdf_output_filename}",
  orography_buffer_file="${buffer_topo}",
  soil_buffer_file="${buffer_soil}",
  lu_buffer_file="${buffer_lu}",
  glcc_buffer_file="${buffer_glcc}",
  flake_buffer_file="${buffer_flake}",
  ndvi_buffer_file="${buffer_ndvi}",
  t_clim_buffer_file="${buffer_tclim}",
  aot_buffer_file="${buffer_aot}",
  alb_buffer_file="${buffer_alb}",
  i_lsm_data=1,
  land_sea_mask_file="",
  number_special_points=0,
/  
EOF_check

# link raw data files to local workdir
ln -s -f ${data_dir}/${raw_data_alb}
ln -s -f ${data_dir}/${raw_data_alnid}
ln -s -f ${data_dir}/${raw_data_aluvd}

ln -s -f ${data_dir}/${raw_data_aot}

ln -s -f ${data_dir}/${raw_data_tclim_coarse}
ln -s -f ${data_dir}/${raw_data_tclim_fine}

ln -s -f ${data_dir}/${raw_data_glc2000}
ln -s -f ${data_dir}/${raw_data_glcc}

ln -s -f ${data_dir}/${raw_data_globcover_0}
ln -s -f ${data_dir}/${raw_data_globcover_1}
ln -s -f ${data_dir}/${raw_data_globcover_2}
ln -s -f ${data_dir}/${raw_data_globcover_3}
ln -s -f ${data_dir}/${raw_data_globcover_4}
ln -s -f ${data_dir}/${raw_data_globcover_5}


ln -s -f ${data_dir}/${raw_data_globe_A10} 
ln -s -f ${data_dir}/${raw_data_globe_B10}
ln -s -f ${data_dir}/${raw_data_globe_C10} 
ln -s -f ${data_dir}/${raw_data_globe_D10} 
ln -s -f ${data_dir}/${raw_data_globe_E10}
ln -s -f ${data_dir}/${raw_data_globe_F10} 
ln -s -f ${data_dir}/${raw_data_globe_G10} 
ln -s -f ${data_dir}/${raw_data_globe_H10} 
ln -s -f ${data_dir}/${raw_data_globe_I10} 
ln -s -f ${data_dir}/${raw_data_globe_J10} 
ln -s -f ${data_dir}/${raw_data_globe_K10} 
ln -s -f ${data_dir}/${raw_data_globe_L10} 
ln -s -f ${data_dir}/${raw_data_globe_M10} 
ln -s -f ${data_dir}/${raw_data_globe_N10} 
ln -s -f ${data_dir}/${raw_data_globe_O10} 
ln -s -f ${data_dir}/${raw_data_globe_P10}

ln -s -f ${data_dir}/${raw_data_aster_T01} 
ln -s -f ${data_dir}/${raw_data_aster_T02}
ln -s -f ${data_dir}/${raw_data_aster_T03} 
ln -s -f ${data_dir}/${raw_data_aster_T04} 
ln -s -f ${data_dir}/${raw_data_aster_T05}
ln -s -f ${data_dir}/${raw_data_aster_T06} 
ln -s -f ${data_dir}/${raw_data_aster_T07} 
ln -s -f ${data_dir}/${raw_data_aster_T08} 
ln -s -f ${data_dir}/${raw_data_aster_T09} 
ln -s -f ${data_dir}/${raw_data_aster_T10} 
ln -s -f ${data_dir}/${raw_data_aster_T11} 
ln -s -f ${data_dir}/${raw_data_aster_T12} 
ln -s -f ${data_dir}/${raw_data_aster_T13} 
ln -s -f ${data_dir}/${raw_data_aster_T14} 
ln -s -f ${data_dir}/${raw_data_aster_T15}
ln -s -f ${data_dir}/${raw_data_aster_T16} 
ln -s -f ${data_dir}/${raw_data_aster_T17} 
ln -s -f ${data_dir}/${raw_data_aster_T18} 
ln -s -f ${data_dir}/${raw_data_aster_T19} 
ln -s -f ${data_dir}/${raw_data_aster_T20} 
ln -s -f ${data_dir}/${raw_data_aster_T21} 
ln -s -f ${data_dir}/${raw_data_aster_T22} 
ln -s -f ${data_dir}/${raw_data_aster_T23} 
ln -s -f ${data_dir}/${raw_data_aster_T24} 
ln -s -f ${data_dir}/${raw_data_aster_T25}
ln -s -f ${data_dir}/${raw_data_aster_T26} 
ln -s -f ${data_dir}/${raw_data_aster_T27} 
ln -s -f ${data_dir}/${raw_data_aster_T28} 
ln -s -f ${data_dir}/${raw_data_aster_T29} 
ln -s -f ${data_dir}/${raw_data_aster_T30} 
ln -s -f ${data_dir}/${raw_data_aster_T31} 
ln -s -f ${data_dir}/${raw_data_aster_T32} 
ln -s -f ${data_dir}/${raw_data_aster_T33} 
ln -s -f ${data_dir}/${raw_data_aster_T34} 
ln -s -f ${data_dir}/${raw_data_aster_T35} 
ln -s -f ${data_dir}/${raw_data_aster_T36} 
ln -s -f ${data_dir}/ASTER_orig_T006.nc 
ln -s -f ${data_dir}/ASTER_orig_T007.nc
ln -s -f ${data_dir}/ASTER_orig_T018.nc
ln -s -f ${data_dir}/ASTER_orig_T019.nc
ln -s -f ${data_dir}/ASTER_orig_T030.nc
ln -s -f ${data_dir}/ASTER_orig_T031.nc 
ln -s -f ${data_dir}/ASTER_orig_T042.nc
ln -s -f ${data_dir}/ASTER_orig_T043.nc

ln -s -f ${data_dir}/${raw_data_ndvi}

ln -s -f ${data_dir}/${raw_data_soil_FAO}
ln -s -f ${data_dir}/${raw_data_soil_HWSD}
ln -s -f ${data_dir}/${raw_data_deep_soil}

ln -s -f ${data_dir}/${raw_lookup_table_HWSD}
ln -s -f ${data_dir}/${raw_HWSD_data}
ln -s -f ${data_dir}/${raw_HWSD_data_deep}

ln -s -f ${data_dir}/${raw_data_flake}

# run the programs
# the next seven programs can run independent of each other
time ${progdir}/${binary_alb}
time ${progdir}/${binary_aot}
time ${progdir}/${binary_tclim}
time ${progdir}/${binary_lu}
time ${progdir}/${binary_topo}
time ${progdir}/${binary_ndvi}
time ${progdir}/${binary_soil}
time ${progdir}/${binary_flake}

# the consistency check requires the output of 
# ${binary_aot}, ${binary_tclim}, ${binary_lu}, ${binary_globe}, 
# ${binary_ndvi}, ${binary_soil} and ${binary_flake}
time ${progdir}/${binary_extpar_consistency_check}

ls
echo "\n>>>> External parameters for COSMO model generated <<<<\n"
