#! /bin/bash
#
#  Interpolate ERA interim monthly means of SST and W_SNOW to ICON grid
#  CAUTION!!! Runs ONLY on LCE not on XCE !!!
#
#____________________________________________________________________________________
#
set -eu

set +u
. $MODULESHOME/init/bash
set -u
#_____________________________________________________________________2_______________

eradir=/work/mh0287/m214089/extpar_T_input
if [[ ! -d ${eradir} ]]
then
    echo "ERROR: Input directory ${eradir} does not exist..."
    exit 1
fi

workdir=/scratch/m/m214089/preprocessing/extpar_160km
workdir=/work/mh0287/users/uwe/luis/nextpar_160km
if [[ ! -d ${workdir} ]]
then
    mkdir -p ${workdir}
fi
cd ${workdir}

export ICON_GRID_DIR=/pool/data/ICON/grids/public/mpim/0013
icon_grid_file=icon_grid_0013_R02B04_G.nc
ires=0013_R02B04_G
ln -sf $ICON_GRID_DIR/$icon_grid_file $icon_grid_file

first_year=1986   # first year for climatology
last_year=2015   # last year for climatology

ifs_file="ei_an${first_year}-${last_year}"

#____________________________________________________________________________________
# get input files
cp ${eradir}/${ifs_file} ${ifs_file}

# calculate monthly mean over all years and adapt to DWD convention on static boundary conditions
cdo -setyear,1111 -setday,11 -ymonmean ${ifs_file} ${ifs_file}.mean

# interpolate from IFS to ICON grid

number_of_threads=2

interpolation_method=nn               # nearest-neighbor scalar

icon_file=$(basename ${ifs_file})_${ires}
rm -f ${icon_file}.nc

# parameter table
cat > parameterTable <<EOF
&parameter
  param=235.128
  name=SKT
  out_name=T_S
  long_name="Skin temperature"
  units="K" 
/
&parameter
  param=34.128
  name=SST
  out_name=T_SEA
  long_name="Sea surface temperature"
  units="K"
/
&parameter
  param=141.128
  name=SD
  out_name=W_SNOW
  long_name="Snow depth"
  units="m of water equivalent"
  standard_name=lwe_thickness_of_surface_snow_amount
/
EOF

ls -l
    
cdo -r -f nc -P ${number_of_threads} setpartabp,parameterTable -remap${interpolation_method},${icon_grid_file} ${ifs_file}.mean ${icon_file}.nc

ls -l ${icon_file}.nc

# beautify the resulting netcdf to use it in EXTPAR consistency check
# that expects a special buffer format time,ie,je,ke rename first
# dimension into ie
ncrename -d ncells,ie ei_an1986-2015_${ires}.nc ei_an1986-2015_${ires}_ie.nc
# add second dimension je -> Caution! ordering of dimensions in
# opposite in NetCDF and Fortran
ncap2 -s 'defdim("je",1);T_SEA_ieje[$time,$je,$ie]=T_SEA' -s 'W_SNOW_ieje[$time,$je,$ie]=W_SNOW' -O ei_an1986-2015_${ires}_ie.nc ei_an1986-2015_${ires}_ieje.nc
#add third dimension ke
ncap2 -s 'defdim("ke",1);T_SEA_iejeke[$time,$ke,$je,$ie]=T_SEA_ieje' -s 'W_SNOW_iejeke[$time,$ke,$je,$ie]=W_SNOW_ieje' -O ei_an1986-2015_${ires}_ieje.nc ei_an1986-2015_${ires}_iejeke.nc
# delete old variables
ncks -C -O -x -v T_SEA,T_SEA_ieje,W_SNOW,W_SNOW_ieje  ei_an1986-2015_${ires}_iejeke.nc ei_an1986-2015_${ires}_iejeke_tmp.nc
# rename variables back
ncrename -h -O -v T_SEA_iejeke,T_SEA -v W_SNOW_iejeke,W_SNOW ei_an1986-2015_${ires}_iejeke_tmp.nc ei_sst_an1986-2015_${ires}_BUFFER.nc
