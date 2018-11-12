#! /bin/bash
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

H_name='HSURF_eiCL'
T_name='T_2M_eiCL'

first_year=1986   # first year for climatology
last_year=2015   # last year for climatology

ifs_file="ei_2t_an${first_year}-${last_year}"
ifs_oro_file="ei_oro_${first_year}"

#____________________________________________________________________________________
#
# get input files
cp ${eradir}/${ifs_file} ${ifs_file}
cp ${eradir}/${ifs_oro_file} ${ifs_oro_file}.tmp

# calculate monthly mean over all years and adapt to DWD convention on static boundary conditions
cdo -setyear,1111 -setday,11 -ymonmean ${ifs_file} ${ifs_file}.mean

cdo infon ${ifs_file}.mean

# convert from geopotential to height
c=1  # to convert from geopotential to heights
cdo -setparam,8.2 -mulc,${c} ${ifs_oro_file}.tmp ${ifs_oro_file}
rm ${ifs_oro_file}.tmp

cdo infon ${ifs_oro_file}

# interpolate from IFS to ICON grid

number_of_threads=2

interpolation_method=nn                # nearest-neighbor scalar

file_basename="eit2icon"

icon_file=$(basename ${ifs_file})_${ires}
rm -f ${icon_file}.nv
        
# parameter table
cat > parameterTable <<EOF
&parameter
  param=167.128
  name=2T
  out_name=${T_name}
  long_name="2 metre temperature"
  units="K"
/
&parameter
  param=8.2
  name=HSURF
  out_name=${H_name}
  long_name="Geometrical height"
  units="m"
/
EOF
    
ls -l

cdo -r -f nc -P ${number_of_threads} setpartabp,parameterTable -remap${interpolation_method},${icon_grid_file} ${ifs_oro_file}  ${icon_file}.00
cdo -r -f nc -P ${number_of_threads} setpartabp,parameterTable -remap${interpolation_method},${icon_grid_file} ${ifs_file}.mean ${icon_file}.nc

ls -l ${icon_file}.nc

#           remove time dimension from orography netCDF file
h_tmp=h.nc
ncwa -a time ${icon_file}.00 ${h_tmp}
#           append orography file to temperature file
ncks -C -A -v HSURF_eiCL ${h_tmp} ${icon_file}.nc
rm ${h_tmp}

# beautify the resulting netcdf to use it in EXTPAR consistency check
# that expects a special buffer format time,ie,je,ke rename first
# dimension into ie
ncrename -d ncells,ie ${icon_file}.nc ${icon_file}_ie.nc
# add second dimension je -> Caution! ordering of dimensions in
# opposite in NetCDF and Fortran
ncap2 -s 'defdim("je",1);T_2M_eiCL_ieje[$time,$je,$ie]=T_2M_eiCL' -s 'HSURF_eiCL_ieje[$je,$ie]=HSURF_eiCL' -O ${icon_file}_ie.nc ${icon_file}_ieje.nc
# add third dimension ke
ncap2 -s 'defdim("ke",1);T_2M_eiCL_iejeke[$time,$ke,$je,$ie]=T_2M_eiCL_ieje' -s 'HSURF_eiCL_iejeke[$ke,$je,$ie]=HSURF_eiCL_ieje' -O ${icon_file}_ieje.nc ${icon_file}_iejeke.nc
# delete old variables
ncks -C -O -x -v T_2M_eiCL,T_2M_eiCL_ieje,HSURF_eiCL,HSURF_eiCL_ieje  ${icon_file}_iejeke.nc ${icon_file}_iejeke_tmp.nc
# rename variables back
ncrename -h -O -v T_2M_eiCL_iejeke,T_2M_CLIM -v HSURF_eiCL_iejeke,TOPO_CLIM ${icon_file}_iejeke_tmp.nc  ei_t2m_an1986-2015_${ires}_BUFFER.nc

