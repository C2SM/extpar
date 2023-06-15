#!/bin/ksh -l

 codes_info

 data_dir=/hpc/rhome/routfor/routfox/extpar/rawdata/

 filename=$1
 progdir=$2
 workdir=$3
 icon_grid_dir=$4

#print -- "run_extpar_ICON_hwsdART: OMP_NUM_THREADS=${OMP_NUM_THREADS}"

 file_date=$(date  +%Y%m%d)

 # set target grid definition 
 icon_grid_file=icon_grid_${filename}.nc

 if [[ ! -d ${workdir} ]] ; then
    mkdir -p ${workdir} 
 fi
 cd ${workdir}
 pwd

 binary_hwsdART=extpar_hwsdART_to_buffer.exe

#---

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

 grib_output_filename="icon_extpar_${filename}_${file_date}.g2"
 netcdf_output_filename="icon_extpar_${filename}_${file_date}.nc"
 grib_sample='GRIB2'

 print -- $netcdf_output_filename
 print -- $grib_output_filename

#----

 raw_data_hwsdART='HWSD0_USDA.nc'
 buffer_hwsdART='hwsdART_extpar_BUFFER.nc'
 output_hwsdART='hwsdART_extpar_ICON.nc'

#---

 cat > INPUT_hwsdART << EOF_hwsdART
&hwsdART_nml
  raw_data_hwsdART_path='${data_dir}',
  raw_data_hwsdART_filename='${raw_data_hwsdART}'
  hwsdART_buffer_file='${buffer_hwsdART}'
  hwsdART_output_file='${output_hwsdART}'
/
EOF_hwsdART

#---

# run the programs

 ${progdir}/${binary_hwsdART}

 print -- 'External parameters for refinement grid of ICON model generated in'
 print --  ${workdir}
 print -- "For plotting results with bplot use: export ICON_COORDINATE_FILE=$icon_grid_dir/$icon_grid_file"

