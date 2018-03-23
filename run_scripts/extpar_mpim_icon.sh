#! /bin/bash
#_____________________________________________________________________________
#SBATCH --job-name=extpar
#SBATCH --partition=gpu
#SBATCH --exclusive
#SBATCH --nodelist=mg204
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=18
#SBATCH --mem=864G
#SBATCH --time=04:00:00
#SBATCH --mail-type=FAIL
#SBATCH --account=mh0287
#SBATCH --output=extpar.o%j
#SBATCH --error=extpar.e%j
#_____________________________________________________________________________

icon_grid_file=icon_grid_0005_R02B04_G.nc
#icon_grid_file=icon_grid_0010_R02B09_G.nc
#icon_grid_file=icon_grid_0011_R02B10_G.nc

icon_grid=$(echo $icon_grid_file | awk -F'[_.]' '{print $3 "_" $4 "_" $5}')

tile_mode=1
urban_mode=0

bindir=/home/mpim/m214089/extpar5/bin

workdir=/scratch/m/m214089/${SLURM_JOBID:-$$}
workdir=/scratch/m/m214089

griddir=/pool/data/ICON/grids/private/mpim/icon_preprocessing/source/preliminary_grids
datadir=/pool/data/ICON/grids/private/mpim/icon_preprocessing/source/extpar_input.2016/

export OMP_NUM_THREADS=1

today=$(date  +%Y%m%d)

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

#nest_TAG=$(basename ${icon_grid_file} .nc | awk -F '_' '{print $5}' | cut -c 1)
#nest_ID=$(basename ${icon_grid_file} .nc | awk -F '_' '{print $5}' | cut -c 3)
#echo "Check to recognize nested domain with nest_TAG=${nest_TAG} and nest_ID=${nest_ID}"

if [[ ! -d ${workdir} ]] ; then
  mkdir -p ${workdir} 
fi
cd ${workdir}
pwd

binary_extpar_consistency_check=extpar_consistency_check.exe
binary_aot=extpar_aot_to_buffer.exe
binary_tclim=extpar_cru_to_buffer.exe
binary_lu=extpar_landuse_to_buffer.exe
binary_topo=extpar_topo_to_buffer.exe
binary_ndvi=extpar_ndvi_to_buffer.exe
binary_soil=extpar_soil_to_buffer.exe
binary_flake=extpar_flake_to_buffer.exe
binary_alb=extpar_alb_to_buffer.exe

cat > INPUT_grid_org << EOF
&GRID_DEF 
 igrid_type = 1,
 domain_def_namelist='INPUT_ICON_GRID'
/ 
EOF

cat > INPUT_ICON_GRID << EOF
&icon_grid_info
  icon_grid_dir='${griddir}'
  icon_grid_nc_file ='${icon_grid_file}'
/
EOF

if [ "${tile_mode}" = 1 ]
then 
    grib_output_filename="icon_extpar_${filename}_${today}_tiles.g2"
    netcdf_output_filename="icon_extpar_${filename}_${today}_tiles.nc"
    echo "Output file is $netcdf_output_filename"
else
    grib_output_filename="icon_extpar_${filename}_${today}.g2"
    netcdf_output_filename="icon_extpar_${filename}_${today}.nc"
    echo "Output file is $netcdf_output_filename"
fi

echo $netcdf_output_filename
echo $grib_output_filename

raw_data_alb='month_alb_new.nc'
raw_data_alnid='month_alnid_new.nc'
raw_data_aluvd='month_aluvd_new.nc'
buffer_alb='month_alb_BUFFER.nc'
output_alb='month_alb_extpar_ICON.nc'

raw_data_aot='aerosol_optical_thickness.nc'
buffer_aot='extpar_aot_BUFFER.nc'
output_aot='aot_extpar_ICON.nc'

raw_data_tclim_coarse='absolute_hadcrut3.nc'
raw_data_tclim_fine='CRU_T2M_SURF_clim.nc'
buffer_tclim='crutemp_climF_extpar_BUFFER.nc'
output_tclim='crutemp_climF_extpar_ICON.nc'

raw_data_glc2000='glc2000_byte.nc'
buffer_glc2000='extpar_landuse_BUFFER.nc'
output_glc2000='extpar_landuse_ICON.nc'
raw_data_glcc='glcc_usgs_class_byte.nc'
buffer_glcc='glcc_landuse_BUFFER.nc'
output_glcc='glcc_landuse_ICON.nc'

# raw_data_globcover='GLOBCOVER_L4_200901_200912_V2.3_int16.nc'
raw_data_globcover_0='GLOBCOVER_0_16bit.nc'
raw_data_globcover_1='GLOBCOVER_1_16bit.nc'
raw_data_globcover_2='GLOBCOVER_2_16bit.nc'
raw_data_globcover_3='GLOBCOVER_3_16bit.nc'
raw_data_globcover_4='GLOBCOVER_4_16bit.nc'
raw_data_globcover_5='GLOBCOVER_5_16bit.nc'
buffer_lu='extpar_landuse_BUFFER.nc'
output_lu='extpar_landuse_ICON.nc'

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

raw_data_ndvi='NDVI_1998_2003.nc'
buffer_ndvi='NDVI_BUFFER.nc'
output_ndvi='ndvi_extpar_ICON.nc'

raw_data_soil_FAO='FAO_DSMW_DP.nc'
raw_data_soil_HWSD='HWSD0_30_texture_2.nc'
raw_data_deep_soil='HWSD30_100_texture_2.nc'
buffer_soil='SOIL_BUFFER.nc'
output_soil='SOIL_ICON.nc'

raw_lookup_table_HWSD='LU_TAB_HWSD_UF.data'
raw_HWSD_data='HWSD_DATA_COSMO.data'
raw_HWSD_data_deep='HWSD_DATA_COSMO_S.data'
raw_HWSD_data_extpar='HWSD_DATA_COSMO_EXTPAR.asc'

raw_data_flake='lakedepth.nc'
buffer_flake='flake_BUFFER.nc'
output_flake='ext_par_flake_ICON.nc'

# NOAA imperviou surface area dataset (default)
raw_data_isa_0='NOAA_ISA_16bit.nc'

# # EEA impervious surface area dataset
# raw_data_isa_0='EEA_ISA_4_16bit.nc'


buffer_isa='ISA_BUFFER.nc'
output_isa='ISA_extpar_ICON.nc'

# this file is adapted from Flanner (2009)
raw_data_ahf='AHF_2006_2.5min_latreverse.nc' 

# # AHF is redistributed at 25km scales according to NOAA ISA.
# raw_data_ahf='AHF_2006_NOAAISAredistr.nc'

buffer_ahf='AHF_BUFFER.nc'
output_ahf='AHF_extpar_ICON.nc'

# create input namelists 
cat > INPUT_AOT << EOF
&aerosol_raw_data
  raw_data_aot_path='${datadir}',
  raw_data_aot_filename='${raw_data_aot}'
/  
&aerosol_io_extpar
  aot_buffer_file='${buffer_aot}',
  aot_output_file='${output_aot}'
/
EOF
#---
cat > INPUT_TCLIM << EOF
&t_clim_raw_data
  raw_data_t_clim_path='${datadir}',
  raw_data_t_clim_filename='${raw_data_tclim_coarse}'
  raw_data_t_id = 1
/  

&t_clim_io_extpar
  t_clim_buffer_file='crutemp_climC_extpar_BUFFER.nc',
  t_clim_output_file='crutemp_climC_extpar_ICON.nc'
/  
EOF
#---
cat > INPUT_LU << EOF
&lu_raw_data
   raw_data_lu_path='',
   raw_data_lu_filename='${datadir}${raw_data_globcover_0}' '${datadir}${raw_data_globcover_1}' '${datadir}${raw_data_globcover_2}' '${datadir}${raw_data_globcover_3}' '${datadir}${raw_data_globcover_4}' '${datadir}${raw_data_globcover_5}',
   i_landuse_data=1,
   ilookup_table_lu=1,
   ntiles_globcover=6
/
&lu_io_extpar
   lu_buffer_file='${buffer_lu}',
   lu_output_file='${output_lu}'
/
&glcc_raw_data
   raw_data_glcc_path='${datadir}',
   raw_data_glcc_filename='${raw_data_glcc}'
/
&glcc_io_extpar
   glcc_buffer_file='${buffer_glcc}',
   glcc_output_file='${output_glcc}'
/
EOF
#---
cat > INPUT_ORO << EOF
&orography_io_extpar
  orography_buffer_file='${buffer_topo}',
  orography_output_file='${output_topo}'
/
&orography_raw_data
 itopo_type = 1
 lsso_param = .TRUE., 
 lfilter_topo = .FALSE.,
 lsubtract_mean_slope = .FALSE., 
 raw_data_orography_path='',
 ntiles_column = 4,
 ntiles_row = 4,
 topo_files =  '${datadir}${raw_data_globe_A10}' '${datadir}${raw_data_globe_B10}'  '${datadir}${raw_data_globe_C10}'  '${datadir}${raw_data_globe_D10}'  '${datadir}${raw_data_globe_E10}'  '${datadir}${raw_data_globe_F10}'  '${datadir}${raw_data_globe_G10}'  '${datadir}${raw_data_globe_H10}'  '${datadir}${raw_data_globe_I10}'  '${datadir}${raw_data_globe_J10}'  '${datadir}${raw_data_globe_K10}'  '${datadir}${raw_data_globe_L10}'  '${datadir}${raw_data_globe_M10}'  '${datadir}${raw_data_globe_N10}'  '${datadir}${raw_data_globe_O10}'  '${datadir}${raw_data_globe_P10}'   
/
EOF
cat > INPUT_OROSMOOTH << EOF
&orography_smoothing
  lfilter_oro=.TRUE.
/
EOF
#---
cat > INPUT_RADTOPO << EOF
&radtopo
  lradtopo=.FALSE.,
  nhori=24,
/
EOF
#---
#---
cat > INPUT_NDVI << EOF
&ndvi_raw_data
  raw_data_ndvi_path='${datadir}',
  raw_data_ndvi_filename='${raw_data_ndvi}'
/  
&ndvi_io_extpar
 ndvi_buffer_file='${buffer_ndvi}',
 ndvi_output_file='${output_ndvi}'
/
EOF
#---
cat > INPUT_SOIL << EOF
&soil_raw_data
 isoil_data = 1,
 ldeep_soil = .false.,
 raw_data_soil_path='${datadir}',
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
EOF
#---
if [ "${urban_mode}" = 1 ]
then 
    cat > INPUT_ISA << EOF
&isa_raw_data
   raw_data_isa_path='',
   raw_data_isa_filename='${datadir}${raw_data_isa_0}'
   ntiles_isa=1
/
&isa_io_extpar
   isa_buffer_file='${buffer_isa}',
   isa_output_file='${output_isa}'
/
EOF
    #---
    cat > INPUT_AHF << EOF
&ahf_raw_data
  raw_data_ahf_path='${datadir}',
  raw_data_ahf_filename='${raw_data_ahf}'
/
&ahf_io_extpar
 ahf_buffer_file='${buffer_ahf}',
 ahf_output_file='${output_ahf}'
/
EOF
fi
#---
cat > INPUT_FLAKE << EOF
&flake_raw_data
   raw_data_flake_path='${datadir}',
   raw_data_flake_filename='${raw_data_flake}'
/
&flake_io_extpar
   flake_buffer_file='${buffer_flake}'
   flake_output_file='${output_flake}'
/
EOF
#---
cat > INPUT_ALB << EOF
&alb_raw_data
  raw_data_alb_path='${datadir}',
  raw_data_alb_filename='${raw_data_alb}'
/
&alnid_raw_data
  raw_data_alb_path='${datadir}',
  raw_data_alnid_filename='${raw_data_alnid}'
/
&aluvd_raw_data
  raw_data_alb_path='${datadir}',
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
EOF
#_______________________________________________________________________________
# consistency check
cat > INPUT_CHECK << EOF
&extpar_consistency_check_io
  grib_output_filename='${grib_output_filename}'
  netcdf_output_filename='${netcdf_output_filename}'
  orography_buffer_file='${buffer_topo}'
  soil_buffer_file='${buffer_soil}'
  lu_buffer_file='${buffer_lu}'
  glcc_buffer_file='${buffer_glcc}'
  flake_buffer_file='${buffer_flake}'
  ndvi_buffer_file='${buffer_ndvi}'
  t_clim_buffer_file='${buffer_tclim}'
  aot_buffer_file='${buffer_aot}'
  alb_buffer_file='${buffer_alb}'
  i_lsm_data=1
  land_sea_mask_file=""
  number_special_points=3
  tile_mode=${tile_mode}
/  
EOF

# Modifications for Falkenberg -> k=646652 in 0026_R3B7 ! Caution in Fortran k=1,...N, in IDL k=0,...N-1!
# R3B07 global: 
#    lon_geo_sp=14.1077,
#    lat_geo_sp=52.1215,
# R3B08 Nest: 
#    lon_geo_sp=14.1222,
#    lat_geo_sp=52.1665,
#
#cat > INPUT_SP_1 << EOF_SP
#&special_points
#    lon_geo_sp=14.1222,
#    lat_geo_sp=52.1665,
#    soiltype_sp=3.,
#    z0_sp=0.03,
#    rootdp_sp=0.6,
#    plcovmn_sp=0.55,
#    plcovmx_sp=0.8,
#    laimn_sp=0.5,
#    laimx_sp=2.5,
#    for_d_sp=-1.,
#    for_e_sp=-1.,
#    fr_land_sp=-1.,
#/
#EOF_SP
# Modifications for Waldstation -> k=646649 in 0026_R3B7
# R3B07 global:
#    lon_geo_sp=13.9894,
#    lat_geo_sp=52.1873,
# R3B08 Nest: 
#    lon_geo_sp=13.9525
#    lat_geo_sp=52.1817,
#
#cat > INPUT_SP_2 << EOF_SP
#&special_points
#    lon_geo_sp=13.9525,
#    lat_geo_sp=52.1817,
#    soiltype_sp=3.,
#    z0_sp=0.91,
#    rootdp_sp=0.6,
#    plcovmn_sp=0.79,
#    plcovmx_sp=0.81,
#    laimn_sp=3.0,
#    laimx_sp=4.0,
#    for_d_sp=-1.,
#    for_e_sp=-1.,
#    fr_land_sp=-1.,
#/
#EOF_SP
# Modifications for MOL -> k=646665 in 0026_R3B7
# R3B07 global:
#    lon_geo_sp=14.313,
#    lat_geo_sp=52.260,
# R3B08 Nest: 
#    lon_geo_sp=14.1199,
#    lat_geo_sp=52.2096,
#
#cat > INPUT_SP_3 << EOF_SP
#&special_points
#    lon_geo_sp=14.1199,
#    lat_geo_sp=52.2096,
#    soiltype_sp=5.,
#    z0_sp=0.13,
#    rootdp_sp=0.6,
#    plcovmn_sp=0.55,
#    plcovmx_sp=0.88,
#    laimn_sp=0.78,
#    laimx_sp=3.18,
#    for_d_sp=-1.,
#    for_e_sp=-1.,
#    fr_land_sp=-1.,
#/
#EOF_SP

#---
cat > INPUT_TCLIM << EOF
&t_clim_raw_data
  raw_data_t_clim_path='${datadir}',
  raw_data_t_clim_filename='${raw_data_tclim_fine}',
  raw_data_t_id = 1
/  

&t_clim_io_extpar
  t_clim_buffer_file='${buffer_tclim}',
  t_clim_output_file='${output_tclim}'
/  
EOF
if [ "${urban_mode}" = 1 ]
then 
    ${bindir}/${binary_ahf}
    ${bindir}/${binary_isa}
fi

${bindir}/${binary_topo}

exit

${bindir}/${binary_soil}
${bindir}/${binary_lu}
${bindir}/${binary_tclim}
${bindir}/${binary_alb}
${bindir}/${binary_aot}
${bindir}/${binary_ndvi}
${bindir}/${binary_flake}

# the consistency check requires the output of 
#   ${binary_aot}, ${binary_tclim}, ${binary_lu}, ${binary_globe}, 
#   ${binary_ndvi}, ${binary_soil} and ${binary_flake}

cat > INPUT_TCLIM_FINAL << EOF
&t_clim_raw_data
  raw_data_t_clim_path='${datadir}',
  raw_data_t_clim_filename='${raw_data_tclim_fine}',
  raw_data_t_id = 1
/  

&t_clim_io_extpar
  t_clim_buffer_file='crutemp_climF_extpar_BUFFER.nc',
  t_clim_output_file='crutemp_climC_extpar_BUFFER.nc'
/  
EOF

${bindir}/${binary_extpar_consistency_check}

#echo 'Generation of global NetCDF attributes' 
#
#sed -n '/Code information/,/code information/p' log_consistency_check.txt > attrib_comment.txt
#library_name=$(fgrep  "Library name" attrib_comment.txt | awk -F ":"  '{print $2}')
#tag_name=$(fgrep  "Tag name" attrib_comment.txt | awk -F ":"  '{print $2}')
#revision_number=$(fgrep  "Revision number" attrib_comment.txt | awk -F ":"  '{print $2}')
#checkin_date=$(fgrep  "Checkin-Date" attrib_comment.txt | awk -F ":"  '{print $2}')
#code_modified=$(fgrep  "Code is modified" attrib_comment.txt | awk -F ":"  '{print $2}')
#compile_date=$(fgrep  "Compile-Date" attrib_comment.txt | awk -F ":"  '{print $2}')
#compiled_by=$(fgrep  "Compiled by" attrib_comment.txt | awk -F ":"  '{print $2}')
#compiled_on=$(fgrep  "Compiled on" attrib_comment.txt | awk -F ":"  '{print $2}')
#start_time=$(fgrep  "Current start time" attrib_comment.txt | awk -F ":"  '{print $2}')
#grib_md5sum=$(md5sum -b ${grib_output_filename})
#
#echo "MD5SUM of the GRIB file is $grib_md5sum"
#
#if [ "${nest_TAG}" = "N" ]
#then
#    if [ "${nest_ID}" = "2" ]
#    then
#        /usr/local/pkg/grib_api/1.12.3/CRAY/bin/grib_set -sgeneratingProcessIdentifier=2 ${grib_output_filename} ${grib_output_filename}_GPI2
#        echo "For nested domain nest_TAG = ${nest_TAG} with generatingProcessIdentifier=${nest_ID} the result file in GRIB2 is ${grib_output_filename}_GPI2"
#    fi
#    if [ "${nest_ID}" = "3" ]
#    then
#        /usr/local/pkg/grib_api/1.12.3/CRAY/bin/grib_set -sgeneratingProcessIdentifier=3 ${grib_output_filename} ${grib_output_filename}_GPI3
#        echo "For nested domain nest_TAG = ${nest_TAG} with generatingProcessIdentifier=${nest_ID} the result file in GRIB2 is ${grib_output_filename}_GPI3"
#    fi
#fi
#
#/e/uhome/jhelmert/bin/nco-4.4.7/src/nco/ncatted -h -a library_name,global,o,c,"$library_name"\
# -a tag_name,global,o,c,"$tag_name"\
# -a comment,global,d,, \
# -a history,global,d,, \
# -a bin_dir,global,o,c,"$bindir"\
# -a revision_number,global,o,c,"$revision_number"\
# -a checkin_date,global,o,c,"$checkin_date"\
# -a code_modified,global,o,c,"$code_modified"\
# -a compile_date,global,o,c,"$compile_date"\
# -a compiled_by,global,o,c,"$compiled_by"\
# -a compiled_on,global,o,c,"$compiled_on"\
# -a start_time,global,o,c,"$start_time"\
# -a tile_mode,global,o,c,"$tile_mode"\
# -a grib_md5sum,global,o,c,"$grib_md5sum" ${netcdf_output_filename}
#
#/e/uhome/jhelmert/bin/nco-4.4.7/src/nco/ncbo -O $workdir/${netcdf_output_filename} ${testfile}  diff.nc; /e/uhome/jhelmert/bin/cdo_xce infov diff.nc
#
echo "Test of external parameters against version ${testfile} has been performed"
echo "Difference to former version using cdo_xce infov diff.nc should be min,max,avg=0 for every parameter!"

echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "External parameters for refinement grid of ICON model generated in $workdir"
echo "$testfile has been used as reference. Please check $workdir/diff.nc"

#echo "Check consistency to former versions with ncdiff -O ${netcdf_output_filename} /e/uscratch/jhelmert/ICON_EXTPAR_20141202/icon_extpar_0026_R03B07_G_20141202.nc diff.nc; cdo infov diff.nc"

#echo "For plotting results with bplot use: export ICON_COORDINATE_FILE=$griddir/$icon_grid_file"
#echo "or modify ${HOME}/.profile_bplot  "
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

