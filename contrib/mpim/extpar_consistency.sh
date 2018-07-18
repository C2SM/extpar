#! /bin/bash
#-----------------------------------------------------------------------------

#SBATCH --job-name=consistency
#SBATCH --partition=gpu
##SBATCH --exclusive
#SBATCH --nodelist=mg208
#SBATCH --nodes=1
#SBATCH --mem=864G
#SBATCH --time=08:00:00
#SBATCH --mail-type=FAIL
#SBATCH --account=mh0287
#SBATCH --output=consistency.o%j
#SBATCH --error=consistency.e%j

#-----------------------------------------------------------------------------

# disable core dumps
ulimit -c 0
# limit stacksize
ulimit -s unlimited

set -eux

#-----------------------------------------------------------------------------

cd $HOME/preprocessing

data_dir=/home/mpim/m214089/icon_preprocessing/source/extpar_input.2016/

filename=0017_R02B10_G

workdir=/scratch/m/m214089/preprocessing/extpar_R02B10_new
icon_grid_dir=/home/mpim/m214089/icon_preprocessing/source/preliminary_grids

progdir=$HOME/preprocessing/extpar/bin
export PATH=/home/mpim/m214003/local/bin:$PATH

file_date=$(date  +%Y%m%d)

tile_mode=1
urban_mode=0

if [ "${tile_mode}" = 1 ]
then
    echo "TILE_MODE for EXTPAR is switched ON"
else
    echo "TILE_MODE for EXTPAR is switched OFF"
fi
if [ "${urban_mode}" = 1 ]
then
    echo "URBAN_MODE for EXTPAR is switched ON"
else
    echo "URBAN_MODE for EXTPAR is switched OFF"
fi

export OMP_NUM_THREADS=18
#export OMP_NUM_THREADS=1
ap_command='time '

icon_grid_file=icon_grid_${filename}.nc

#nest_TAG=$(basename ${icon_grid_file} .nc | awk -F '_' '{print $5}' | cut -c 1)
#nest_ID=$(basename ${icon_grid_file} .nc | awk -F '_' '{print $5}' | cut -c 3)

#echo "Check to recognize nested domain with nest_TAG=${nest_TAG} and nest_ID=${nest_ID}"

if [[ ! -d ${workdir} ]] ; then
  mkdir -p ${workdir} 
fi
cd ${workdir}
pwd 

## PreProcessing of ERA-I data now as part of the main script
#cp ${data_dir}/ERAinterim/REMAP/ei_an1986-2015_${filename}_BUFFER.nc ${workdir}
#cp ${data_dir}/ERAinterim/REMAP/ei_2t_an1986-2015_${filename}_BUFFER.nc ${workdir}

#For Binaries in /usr/local/pkg/for0adm/abs/
#--------------------------------------------
#binary_extpar_consistency_check=extpar_consistency_check.new
#binary_aot=extpar_aot_to_buffer.new
#binary_tclim=extpar_cru_to_buffer.new
#binary_lu=extpar_landuse_to_buffer.new
#binary_topo=extpar_topo_to_buffer.new
#binary_ndvi=extpar_ndvi_to_buffer.new
#binary_soil=extpar_soil_to_buffer.new
#binary_flake=extpar_flake_to_buffer.new
#binary_alb=extpar_alb_to_buffer.new
#binary_isa=extpar_isa_to_buffer.new
#binary_ahf=extpar_ahf_to_buffer.new


#For Binaries in /e/uhome/jhelmert/EXTPAR_V2_1
#--------------------------------------------
#binary_extpar_consistency_check=extpar_consistency_check
#binary_aot=tstextpar_aot_to_buffer
#binary_tclim=tstextpar_cru_to_buffer
#binary_lu=tstextpar_landuse_to_buffer
#binary_topo=tstextpar_topo_to_buffer
#binary_ndvi=tstextpar_ndvi_to_buffer
#binary_soil=tstextpar_soil_to_buffer
#binary_flake=tstextpar_flake_to_buffer
#binary_alb=tstextpar_alb_to_buffer
#binary_isa=tstextpar_isa_to_buffer
#binary_ahf=tstextpar_ahf_to_buffer

#For Binaries in /e/uhome/jhelmert/EXTPAR_CSCS_DWD/trunk/bin
#--------------------------------------------
#binary_extpar_consistency_check=extpar_consistency_check.exe
#binary_aot=extpar_aot_to_buffer.exe
#binary_tclim=extpar_cru_to_buffer.exe
#binary_lu=extpar_landuse_to_buffer.exe
#binary_topo=extpar_topo_to_buffer.exe
#binary_ndvi=extpar_ndvi_to_buffer.exe
#binary_soil=extpar_soil_to_buffer.exe
#binary_flake=extpar_flake_to_buffer.exe
#binary_alb=extpar_alb_to_buffer.exe

#For Binaries in /home/mpim/m214089/preprocessing/extpar/bin
#--------------------------------------------
binary_extpar_consistency_check=extpar_consistency_check.exe
binary_aot=extpar_aot_to_buffer.exe
binary_tclim=extpar_cru_to_buffer.exe
binary_lu=extpar_landuse_to_buffer.exe
binary_topo=extpar_topo_to_buffer.exe
binary_ndvi=extpar_ndvi_to_buffer.exe
binary_soil=extpar_soil_to_buffer.exe
binary_flake=extpar_flake_to_buffer.exe
binary_alb=extpar_alb_to_buffer.exe
binary_isa=extpar_isa_to_buffer.exe
binary_ahf=extpar_ahf_to_buffer.exe


cat > INPUT_grid_org << EOF_go
&GRID_DEF 
 igrid_type = 1,
 domain_def_namelist='INPUT_ICON_GRID'
/ 
EOF_go


cat > INPUT_ICON_GRID << EOF5
&icon_grid_info
  icon_grid_dir='${icon_grid_dir}'
  icon_grid_nc_file ='icon_grid_${filename}.nc'
/
EOF5

#---

if [ "${tile_mode}" = 1 ]
then 
    grib_output_filename="icon_extpar_${filename}_${file_date}_tiles.g2"
    netcdf_output_filename="icon_extpar_${filename}_${file_date}_tiles.nc"
    #testfile=$(ls -alrt /e/rhome/routfor/routfox/icon/grids/public/edzw/icon_extpar_${filename}_*.nc| tail -1| awk '{print $9}')
    #testfile=$(ls -alrt /lustre0/uwork1/jhelmert/ICON_EXTPAR_20150723/T1/icon_extpar_${filename}_*.nc| tail -1| awk '{print $9}')
    echo "Output file is $netcdf_output_filename"
else
    grib_output_filename="icon_extpar_${filename}_${file_date}.g2"
    netcdf_output_filename="icon_extpar_${filename}_${file_date}.nc"
    #testfile=$(ls -alrt /e/rhome/routfor/routfox/icon/grids/public/edzw/icon_extpar_${filename}_2015*.nc| tail -1| awk '{print $9}')
    #testfile=$(ls -alrt /lustre0/uwork1/jhelmert/ICON_EXTPAR_20150723/T0/icon_extpar_${filename}_*.nc| tail -1| awk '{print $9}')
    echo "Output file is $netcdf_output_filename"
fi

grib_sample='GRIB2'

echo $netcdf_output_filename
echo $grib_output_filename
#----


# raw data: AVAILABLE
raw_data_alb='month_alb_new.nc'
raw_data_alnid='month_alnid_new.nc'
raw_data_aluvd='month_aluvd_new.nc'
buffer_alb='month_alb_BUFFER.nc'
output_alb='month_alb_extpar_ICON.nc'

# raw data: AVAILABLE
raw_data_aot='aerosol_optical_thickness.nc'
buffer_aot='extpar_aot_BUFFER.nc'
output_aot='aot_extpar_ICON.nc'

# raw data: AVAILABLE
raw_data_tclim_coarse='absolute_hadcrut3.nc'
raw_data_tclim_fine='CRU_T2M_SURF_clim.nc'
buffer_tclim='crutemp_climF_extpar_BUFFER.nc'
output_tclim='crutemp_climF_extpar_ICON.nc'

# raw data: AVAILABLE
raw_data_glc2000='glc2000_byte.nc'
buffer_glc2000='extpar_landuse_BUFFER.nc'
output_glc2000='extpar_landuse_ICON.nc'
raw_data_glcc='glcc_usgs_class_byte.nc'
buffer_glcc='glcc_landuse_BUFFER.nc'
output_glcc='glcc_landuse_ICON.nc'

# raw_data_globcover='GLOBCOVER_L4_200901_200912_V2.3_int16.nc'
# raw data: AVAILABLE
raw_data_globcover_0='GLOBCOVER_0_16bit.nc'
raw_data_globcover_1='GLOBCOVER_1_16bit.nc'
raw_data_globcover_2='GLOBCOVER_2_16bit.nc'
raw_data_globcover_3='GLOBCOVER_3_16bit.nc'
raw_data_globcover_4='GLOBCOVER_4_16bit.nc'
raw_data_globcover_5='GLOBCOVER_5_16bit.nc'
buffer_lu='extpar_landuse_BUFFER.nc'
output_lu='extpar_landuse_ICON.nc'

#raw_data_globe_A10='GLOBE_A_filt_lanczos_window.nc'
#raw_data_globe_B10='GLOBE_B_filt_lanczos_window.nc'
#raw_data_globe_C10='GLOBE_C_filt_lanczos_window.nc'
#raw_data_globe_D10='GLOBE_D_filt_lanczos_window.nc'
#raw_data_globe_E10='GLOBE_E_filt_lanczos_window.nc'
#raw_data_globe_F10='GLOBE_F_filt_lanczos_window.nc'
#raw_data_globe_G10='GLOBE_G_filt_lanczos_window.nc'
#raw_data_globe_H10='GLOBE_H_filt_lanczos_window.nc'
#raw_data_globe_I10='GLOBE_I_filt_lanczos_window.nc'
#raw_data_globe_J10='GLOBE_J_filt_lanczos_window.nc'
#raw_data_globe_K10='GLOBE_K_filt_lanczos_window.nc'
#raw_data_globe_L10='GLOBE_L_filt_lanczos_window.nc'
#raw_data_globe_M10='GLOBE_M_filt_lanczos_window.nc'
#raw_data_globe_N10='GLOBE_N_filt_lanczos_window.nc'
#raw_data_globe_O10='GLOBE_O_filt_lanczos_window.nc'
#raw_data_globe_P10='GLOBE_P_filt_lanczos_window.nc'

# raw data: AVAILABLE
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

buffer_topo='topography_BUFFER.nc'
output_topo='topography_ICON.nc'

# raw data: AVAILABLE
raw_data_ndvi='NDVI_1998_2003.nc'
buffer_ndvi='NDVI_BUFFER.nc'
output_ndvi='ndvi_extpar_ICON.nc'

# raw data: AVAILABLE
raw_data_soil_FAO='FAO_DSMW_DP.nc'
raw_data_soil_HWSD='HWSD0_30_texture_2.nc'
raw_data_deep_soil='HWSD30_100_texture_2.nc'
buffer_soil='SOIL_BUFFER.nc'
output_soil='SOIL_ICON.nc'

# raw data: AVAILABLE
raw_lookup_table_HWSD='LU_TAB_HWSD_UF.data'
raw_HWSD_data='HWSD_DATA_COSMO.data'
raw_HWSD_data_deep='HWSD_DATA_COSMO_S.data'
raw_HWSD_data_extpar='HWSD_DATA_COSMO_EXTPAR.asc'

# raw data: AVAILABLE
raw_data_flake='lakedepth.nc'
buffer_flake='flake_BUFFER.nc'
output_flake='ext_par_flake_ICON.nc'

# NOAA imperviou surface area dataset (default)
# raw data: AVAILABLE
raw_data_isa_0='NOAA_ISA_16bit.nc'
# # EEA impervious surface area dataset
# raw_data_isa_0='EEA_ISA_4_16bit.nc'
buffer_isa='ISA_BUFFER.nc'
output_isa='ISA_extpar_ICON.nc'

# this file is adapted from Flanner (2009)
# raw data: AVAILABLE
raw_data_ahf='AHF_2006_2.5min_latreverse.nc' 
# # AHF is redistributed at 25km scales according to NOAA ISA.
# raw_data_ahf='AHF_2006_NOAAISAredistr.nc'
buffer_ahf='AHF_BUFFER.nc'
output_ahf='AHF_extpar_ICON.nc'

# create input namelists 
cat > INPUT_AOT << EOF_aot
&aerosol_raw_data
  raw_data_aot_path='${data_dir}',
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
  raw_data_t_clim_path='${data_dir}',
  raw_data_t_clim_filename='${raw_data_tclim_coarse}'
  raw_data_t_id = 2
/  

&t_clim_io_extpar
  t_clim_buffer_file='crutemp_climC_extpar_BUFFER.nc',
  t_clim_output_file='crutemp_climC_extpar_ICON.nc'
/  
EOF_tclim
#---
cat > INPUT_LU << EOF_lu
&lu_raw_data
   raw_data_lu_path='${data_dir}',
   raw_data_lu_filename='${raw_data_globcover_0}' '${raw_data_globcover_1}' '${raw_data_globcover_2}' '${raw_data_globcover_3}' '${raw_data_globcover_4}' '${raw_data_globcover_5}',
   i_landuse_data=1,
   ilookup_table_lu=1,
   ntiles_globcover=6
/
&lu_io_extpar
   lu_buffer_file='${buffer_lu}',
   lu_output_file='${output_lu}'
/
&glcc_raw_data
   raw_data_glcc_path='${data_dir}',
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
 itopo_type = 1
 lsso_param = .TRUE., 
 lfilter_topo = .FALSE., 
 lsubtract_mean_slope = .FALSE., 
 raw_data_orography_path='${data_dir}',      
 ntiles_column = 4,
 ntiles_row = 4,
 topo_files =  '${raw_data_globe_A10}' '${raw_data_globe_B10}'  '${raw_data_globe_C10}'  '${raw_data_globe_D10}'  '${raw_data_globe_E10}'  '${raw_data_globe_F10}'  '${raw_data_globe_G10}'  '${raw_data_globe_H10}'  '${raw_data_globe_I10}'  '${raw_data_globe_J10}'  '${raw_data_globe_K10}'  '${raw_data_globe_L10}'  '${raw_data_globe_M10}'  '${raw_data_globe_N10}'  '${raw_data_globe_O10}'  '${raw_data_globe_P10}'   
/
EOF_oro
cat > INPUT_ASTER << EOF_oro
&orography_io_extpar
  orography_buffer_file='${buffer_topo}',
  orography_output_file='${output_topo}'
/
&orography_raw_data
 itopo_type = 2,
 lsso_param = .TRUE.,
 lfilter_topo = .FALSE., 
 lsubtract_mean_slope = .FALSE., 
 raw_data_orography_path='${data_dir}',
 ntiles_column = 3,
 ntiles_row = 4,
 topo_files = 'topo.ASTER_orig_T006.nc' 'topo.ASTER_orig_T007.nc' 'topo.ASTER_orig_T008.nc'  'topo.ASTER_orig_T018.nc' 'topo.ASTER_orig_T019.nc' 'topo.ASTER_orig_T020.nc' 'topo.ASTER_orig_T030.nc' 'topo.ASTER_orig_T031.nc' 'topo.ASTER_orig_T032.nc' 'topo.ASTER_orig_T042.nc' 'topo.ASTER_orig_T043.nc' 'topo.ASTER_orig_T044.nc' 'topo.ASTER_orig_T054.nc' 'topo.ASTER_orig_T055.nc' 'topo.ASTER_orig_T056.nc'
/
EOF_oro

cat > INPUT_OROSMOOTH << EOF_orosmooth
&orography_smoothing
  lfilter_oro=.TRUE.
/
EOF_orosmooth
#---
cat > INPUT_RADTOPO << EOF_rad
&radtopo
  lradtopo=.FALSE.,
  nhori=24,
/
EOF_rad
#---
#---
cat > INPUT_NDVI << EOF_ndvi
&ndvi_raw_data
  raw_data_ndvi_path='${data_dir}',
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
 ldeep_soil = .false.,
 raw_data_soil_path='${data_dir}',
 raw_data_soil_filename='${raw_data_soil_FAO}'
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
if [ "${urban_mode}" = 1 ]
 then 
cat > INPUT_ISA << EOF_isa
&isa_raw_data
   raw_data_isa_path='${data_dir}',
   raw_data_isa_filename='${raw_data_isa_0}'
   ntiles_isa=1
/
&isa_io_extpar
   isa_buffer_file='${buffer_isa}',
   isa_output_file='${output_isa}'
/
EOF_isa


cat > INPUT_AHF << EOF_ahf
&ahf_raw_data
  raw_data_ahf_path='${data_dir}',
  raw_data_ahf_filename='${raw_data_ahf}'
/
&ahf_io_extpar
 ahf_buffer_file='${buffer_ahf}',
 ahf_output_file='${output_ahf}'
/
EOF_ahf
fi
#---

cat > INPUT_FLAKE << EOF_flake
&flake_raw_data
   raw_data_flake_path='${data_dir}',
   raw_data_flake_filename='${raw_data_flake}'
/
&flake_io_extpar
   flake_buffer_file='${buffer_flake}'
   flake_output_file='${output_flake}'
/
EOF_flake

#---

cat > INPUT_ALB << EOF_alb
&alb_raw_data
  raw_data_alb_path='${data_dir}',
  raw_data_alb_filename='${raw_data_alb}'
/
&alnid_raw_data
  raw_data_alb_path='${data_dir}',
  raw_data_alnid_filename='${raw_data_alnid}'
/
&aluvd_raw_data
  raw_data_alb_path='${data_dir}',
  raw_data_aluvd_filename='${raw_data_aluvd}'
/
&alb_io_extpar
  alb_buffer_file='${buffer_alb}',
  alb_output_file='${output_alb}'
/
&alb_source_file     
  alb_source='al'
  alnid_source='alnid'
  aluvd_source='aluvd'
/
EOF_alb


# consistency check
cat > INPUT_CHECK << EOF_check
&extpar_consistency_check_io
  grib_output_filename='${grib_output_filename}',
  netcdf_output_filename='${netcdf_output_filename}',
  orography_buffer_file='${buffer_topo}',
  soil_buffer_file='${buffer_soil}',
  lu_buffer_file='${buffer_lu}',
  glcc_buffer_file='${buffer_glcc}',
  flake_buffer_file='${buffer_flake}',
  ndvi_buffer_file='${buffer_ndvi}',
  sst_icon_file='ei_an1986-2015_${filename}_BUFFER.nc'
  t2m_icon_file='ei_2t_an1986-2015_${filename}_BUFFER.nc'
  t_clim_buffer_file='${buffer_tclim}',
  aot_buffer_file='${buffer_aot}',
  alb_buffer_file='${buffer_alb}',
  i_lsm_data=1,
  land_sea_mask_file="",
  number_special_points=3,
  tile_mode=${tile_mode},
/  
EOF_check

# Modifications for Falkenberg -> k=646652 in 0026_R3B7 ! Caution in Fortran k=1,...N, in IDL k=0,...N-1!
# R3B07 global: 
#    lon_geo_sp1=14.1077
#    lat_geo_sp1=52.1215
# R3B08 Nest + global: 
    lon_geo_sp1=14.1222
    lat_geo_sp1=52.1665

cat > INPUT_SP_1 << EOF_SP
&special_points
    lon_geo_sp=${lon_geo_sp1},
    lat_geo_sp=${lat_geo_sp1},
    soiltype_sp=3.,
    z0_sp=0.03,
    rootdp_sp=0.6,
    plcovmn_sp=0.55,
    plcovmx_sp=0.8,
    laimn_sp=0.5,
    laimx_sp=2.5,
    for_d_sp=-1.,
    for_e_sp=-1.,
    fr_land_sp=-1.,
/
EOF_SP
# Modifications for Waldstation -> k=646649 in 0026_R3B7
# R3B07 global:
#    lon_geo_sp2=13.9894
#    lat_geo_sp2=52.1873
# R3B08 Nest + global : 
    lon_geo_sp2=13.9525
    lat_geo_sp2=52.1817

cat > INPUT_SP_2 << EOF_SP
&special_points
    lon_geo_sp=${lon_geo_sp2},
    lat_geo_sp=${lat_geo_sp2},
    soiltype_sp=3.,
    z0_sp=0.91,
    rootdp_sp=0.6,
    plcovmn_sp=0.79,
    plcovmx_sp=0.81,
    laimn_sp=3.0,
    laimx_sp=4.0,
    for_d_sp=-1.,
    for_e_sp=-1.,
    fr_land_sp=-1.,
/
EOF_SP


# Modifications for MOL -> k=646665 in 0026_R3B7
# R3B07 global :
#    lon_geo_sp3=14.313
#    lat_geo_sp3=52.260
# R3B08 Nest old -> new for IGLO: 
    lon_geo_sp3=14.2320
    lat_geo_sp3=52.2418
# R3B08 Nest new according to meteogram list IEU: 
#    lon_geo_sp3=14.0702
#    lat_geo_sp3=52.2055

cat > INPUT_SP_3 << EOF_SP
&special_points
    lon_geo_sp=${lon_geo_sp3},
    lat_geo_sp=${lat_geo_sp3},
    soiltype_sp=5.,
    z0_sp=0.13,
    rootdp_sp=0.6,
    plcovmn_sp=0.55,
    plcovmx_sp=0.88,
    laimn_sp=0.78,
    laimx_sp=3.18,
    for_d_sp=-1.,
    for_e_sp=-1.,
    fr_land_sp=-1.,
/
EOF_SP

#243	LON	31.156768	LAT	32.058631
#219	LON	30.106029	LAT	31.417649
#227	LON	30.45778	LAT	31.064325
#242	LON	31.102551	LAT	31.29895
#250	LON	31.461231	LAT	30.952482
#214	LON	29.72162	LAT	30.122474
#226	LON	30.386568	LAT	30.345565
#235	LON	30.7384	LAT	29.973466
#246	LON	31.395053	LAT	30.201613



# run the programs:
#---
$ap_command  ${progdir}/${binary_tclim}
#---

# update namelist for consistency check ...
cat > INPUT_TCLIM << EOF_tclim
&t_clim_raw_data
  raw_data_t_clim_path='${data_dir}',
  raw_data_t_clim_filename='${raw_data_tclim_fine}',
  raw_data_t_id = 1
/  

&t_clim_io_extpar
  t_clim_buffer_file='${buffer_tclim}',
  t_clim_output_file='${output_tclim}'
/  
EOF_tclim

if [ "${urban_mode}" = 1 ]
 then 
     #---
     $ap_command ${progdir}/${binary_ahf}
     #---
     $ap_command ${progdir}/${binary_isa}
     #---
fi
#---
$ap_command  ${progdir}/${binary_topo}
#---
$ap_command  ${progdir}/${binary_soil}
#---
$ap_command  ${progdir}/${binary_lu}
#---
$ap_command  ${progdir}/${binary_tclim}
#---
# do not use the extpar preprocessor for the albedo as it is terrible
# slow and the algorithm is insufficient for ICON.
# $ap_command ${progdir}/${binary_alb}
# the following procedure has been proven to be much faster for high resolution:
rm -f weights.nc alb-dis.nc alnid-dis.nc aluvd.nc
time cdo -f nc4 -P $OMP_NUM_THREADS gendis,$icon_grid_dir/$icon_grid_file ${data_dir}${raw_data_alb} weights.nc
time cdo -f nc4 -P $OMP_NUM_THREADS remap,$icon_grid_dir/$icon_grid_file,weights.nc ${data_dir}${raw_data_alb} alb-dis.nc
time cdo -f nc4 -P $OMP_NUM_THREADS remap,$icon_grid_dir/$icon_grid_file,weights.nc ${data_dir}${raw_data_alnid} alnid-dis.nc
time cdo -f nc4 -P $OMP_NUM_THREADS remap,$icon_grid_dir/$icon_grid_file,weights.nc ${data_dir}${raw_data_aluvd} aluvd-dis.nc
time ${progdir}/cdo2alb-buffer.py
mv alb-dis_BUFFER.nc ${buffer_alb} 
#---
$ap_command  ${progdir}/${binary_aot}
#---
# do not use the extpar preprocessor for NDVI as it is terrible
# slow and a bilinear interpolation only which can be done much much faster for high resolution
# ====> $ap_command  ${progdir}/${binary_ndvi}
rm -f weights.nc ndvi-ycon.nc
time cdo -f nc4 -P $OMP_NUM_THREADS genycon,$icon_grid_dir/$icon_grid_file ${data_dir}${raw_data_ndvi} weights.nc
time cdo -f nc4 -P $OMP_NUM_THREADS settaxis,1111-01-01,0,1mo -remap,$icon_grid_dir/$icon_grid_file,weights.nc ${data_dir}${raw_data_ndvi} ndvi-ycon.nc
time ${progdir}/cdo2ndvi-buffer.py
mv ndvi-ycon_BUFFER.nc ${buffer_ndvi}
#---
$ap_command  ${progdir}/${binary_flake}
#---
exit
# the consistency check requires the output of ${binary_aot},
#   ${binary_tclim}, ${binary_lu}, ${binary_globe}, ${binary_ndvi},
#   ${binary_soil}, and ${binary_flake}

printenv | grep GRIB_

cat > INPUT_TCLIM_FINAL << EOF_tclim
&t_clim_raw_data
  raw_data_t_clim_path='${data_dir}',
  raw_data_t_clim_filename='${raw_data_tclim_fine}',
  raw_data_t_id = 1
/  

&t_clim_io_extpar
  t_clim_buffer_file='crutemp_climF_extpar_BUFFER.nc',
  t_clim_output_file='crutemp_climC_extpar_BUFFER.nc'
/  
EOF_tclim

#$ap_command ${progdir}/${binary_extpar_consistency_check}  | tee log_consistency_check.txt

# this is necessary to prevent size mismatch errors for nested domains
#rm -f glcc_landuse_BUFFER.nc
