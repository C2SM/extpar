#! /bin/bash
#_________________________________________________________________________________________________
#SBATCH --job-name=ifs2icon
#SBATCH --partition=gpu
#SBATCH --exclusive
##SBATCH --nodelist=mg204
#SBATCH --ntasks-per-node=72
#SBATCH --nodes=1
#SBATCH --mem=0
#SBATCH --time=12:00:00
#SBATCH --mail-type=FAIL
#SBATCH --account=mh0287
#SBATCH --output=LOG.ifs2icon.run.%j
#SBATCH --error=LOG.ifs2icon.run.%j
#_________________________________________________________________________________________________
#
set -eu

ulimit -s unlimited
ulimit -c 0
#_________________________________________________________________________________________________


#declare -A grid_ids=( ["0019_R02B05"]="80km" )
declare -A grid_ids=( ["0015_R02B09"]="5km" )
#                      ["0021_R02B06"]="40km"
#                      ["0023_R02B07"]="20km" )
#declare -A grid_ids=( ["0025_R02B08"]="10km"
#                      ["0015_R02B09"]="5km" )

ifs_grid_file=/work/mh0287/m214089/dyamond_ecmwf_initial/ifs_oper_T1279_2016080100.grb
ifs_data_file=/work/mh0287/m214089/dyamond_ecmwf_initial/ifs_oper_T1279_2016080100.grb
    
for grid_id in ${!grid_ids[@]}
do

work_dir=/scratch/m/m214089/preprocessing/dyamond/${grid_ids[$grid_id]}
output_dir=/home/mpim/m214089/experiments/input/${grid_ids[$grid_id]}

# for the icontools
rootdir=/mnt/lustre01/sw/rhel6-x64/icontools/dwd_icon_tools
progdir=${rootdir}/icontools

# for mpirun
bindir=/mnt/lustre01/sw/rhel6-x64/icontools/gcc/local.gcc/bin

export PATH=${progdir}:${bindir}:$PATH

cd $work_dir

icon_grid_file=icon_grid_${grid_id}_G.nc
grid_dir=${grid_id%%_*}

if [[ ! -f "$icon_grid_file"  ]]
then
    if [[ -d /pool/data/ICON/grids/public/mpim/$grid_dir ]]
    then
        cp /pool/data/ICON/grids/public/mpim/$grid_dir/$icon_grid_file .
    else
        wget http://icon-downloads.mpimet.mpg.de/grids/public/mpim/$grid_dir/$icon_grid_file
    fi
fi

icon_grid_dir=${work_dir}

#_________________________________________________________________________________________________
#
export OMP_NUM_THREADS=8
export OMP_SCHEDULE="static"
export OMP_DYNAMIC="false"
export OMP_STACKSIZE=4G

export GRIB_DEFINITION_PATH=$HOME/local.gcc/share/eccodes/dwd/definitions:$HOME/local.gcc/share/eccodes/definitions
#_________________________________________________________________________________________________

icon_data_file="${work_dir}/icon_oper_analysis.nc"

if [[ -e "${icon_data_file}" ]]
then
    rm ${icon_data_file}
fi

cat > ${rootdir}/example/tmp.nml << REMAP_NML_EOF
! REMAPPING NAMELIST FILE
!
&remap_nml
 in_grid_filename   = "${ifs_grid_file}"
 in_filename        = "${ifs_data_file}"
 in_type            = 1
 !
 out_grid_filename  = "${icon_grid_file}"
 out_filename       = "${icon_data_file}"
 out_type           = 2
 !
 out_filetype       = 5
 l_have3dbuffer     = .false.
 !
/
! DEFINITION FOR INPUT DATA FIELD
!
&input_field_nml  ! sea surface temperature
 inputname      = "SST"
 outputname     = "SST" 
 code           = 34
 intp_method    = 3    
/
&input_field_nml  ! temperature
 inputname      = "T"         
 outputname     = "T"          
 code           = 130  
 intp_method    = 3            
/
!&input_field_nml  ! horiz. wind comp. u
! inputname      = "U","V"         
! outputname     = "VN"
! code           = 131 
! intp_method    = 3             
!/
&input_field_nml  
 inputname      = "U" 
 outputname     = "U"
 code           = 131          
 intp_method    = 3    
/
&input_field_nml  
 inputname      = "V" 
 outputname     = "V"
 code           = 132          
 intp_method    = 3    
/
&input_field_nml  ! vertical velocity
 inputname      = "W" 
 outputname     = "W"
 code           = 135
 intp_method    = 3    
/
&input_field_nml  ! surface pressure
 inputname      = "LNSP"
 outputname     = "LNPS" 
 code           = 152
 intp_method    = 3    
/
&input_field_nml  ! geopotential
 inputname      = "Z"
 outputname     = "GEOP_SFC"
 code           = 129
 intp_method    = 3    
/
&input_field_nml  ! geopotential
 inputname      = "Z"
 outputname     = "GEOP_ML"
 code           = 129
 intp_method    = 3    
/ 
&input_field_nml  ! specific humidity
 inputname      = "Q"
 outputname     = "QV" 
 code           = 133
 intp_method    = 3    
/
&input_field_nml  ! cloud liquid water content
 inputname      = "CLWC"
 outputname     = "QC" 
 code           = 246
 intp_method    = 3    
/
&input_field_nml  ! cloud ice content
 inputname      = "CIWC"
 outputname     = "QI" 
 code           = 247 
 intp_method    = 3    
/
&input_field_nml  ! rain water content
 inputname      = "CRWC"
 outputname     = "QR" 
 code           = 75 
 intp_method    = 3    
/
&input_field_nml  ! snow water content
 inputname      = "CSWC"
 outputname     = "QS" 
 code           = 76
 intp_method    = 3    
/
&input_field_nml  ! snow temperature
 inputname      = "TSN"
 outputname     = "T_SNOW" 
 code           = 238
 intp_method    = 3    
/
&input_field_nml  ! water content of snow
 inputname      = "SD"
 outputname     = "W_SNOW" 
 code           = 141
 intp_method    = 3    
/
&input_field_nml  ! density of snow
 inputname      = "RSN"
 outputname     = "RHO_SNOW" 
 code           = 33
 intp_method    = 3    
/
&input_field_nml  ! snow albedo
 inputname      = "ASN"
 outputname     = "ALB_SNOW" 
 code           = 32
 intp_method    = 3    
/
&input_field_nml  ! skin temperature
 inputname      = "SKT"
 outputname     = "SKT" 
 code           = 235
 intp_method    = 3    
/
&input_field_nml  ! soil temperature level 1
 inputname      = "STL1"
 outputname     = "STL1" 
 code           = 139
 intp_method    = 3    
/
&input_field_nml  ! soil temperature level 2
 inputname      = "STL2"
 outputname     = "STL2" 
 code           = 170
 intp_method    = 3    
/
&input_field_nml  ! soil temperature level 3
 inputname      = "STL3"
 outputname     = "STL3" 
 code           = 183
 intp_method    = 3    
/
&input_field_nml  ! soil temperature level 4
 inputname      = "STL4"
 outputname     = "STL4" 
 code           = 236
 intp_method    = 3    
/
&input_field_nml  ! sea-ice cover
 inputname      = "CI"
 outputname     = "CI" 
 code           = 31
 intp_method    = 3    
/
&input_field_nml  ! water cont. of interception storage
 inputname      = "SRC"
 outputname     = "W_I" 
 code           = 198 
 intp_method    = 3    
/
&input_field_nml  ! surface roughness
 inputname      = "SR"
 outputname     = "Z0" 
 code           = 173
 intp_method    = 3    
/
&input_field_nml  ! Land/sea mask
 inputname      = "LSM"
 outputname     = "LSM" 
 code           = 172 
 intp_method    = 3    
/
&input_field_nml  ! soil moisture index layer 1
 inputname      = "SWVL1"
 outputname     = "SMIL1" 
 code           = 39
 intp_method    = 3    
 var_in_mask    = "LSM"
 code_in_mask   = 172
 in_mask_threshold = 0.5
 in_mask_below  = .TRUE.
/
&input_field_nml  ! soil moisture index layer 2
 inputname      = "SWVL2"
 outputname     = "SMIL2" 
 code           = 40
 intp_method    = 3    
 var_in_mask    = "LSM"
 code_in_mask   = 172
 in_mask_threshold = 0.5
 in_mask_below  = .TRUE.
/
&input_field_nml  ! soil moisture index layer 3
 inputname      = "SWVL3"
 outputname     = "SMIL3" 
 code           = 41
 intp_method    = 3    
 var_in_mask    = "LSM"
 code_in_mask   = 172
 in_mask_threshold = 0.5
 in_mask_below  = .TRUE.
/
&input_field_nml  ! soil moisture index layer 4
 inputname      = "SWVL4"
 outputname     = "SMIL4" 
 code           = 42
 intp_method    = 3    
 var_in_mask    = "LSM"
 code_in_mask   = 172
 in_mask_threshold = 0.5
 in_mask_below  = .TRUE.
/
REMAP_NML_EOF
#_________________________________________________________________________________________________

mpiexec -np 1 ${rootdir}/icontools/iconremap_mpi --remap_nml ${rootdir}/example/tmp.nml -v

target="ifs2icon_2016080100_${grid_id}_G.nc"
mv ${icon_data_file} ${output_dir}/${target}

#_________________________________________________________________________________________________
#
done
