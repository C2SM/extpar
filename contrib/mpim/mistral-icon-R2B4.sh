#! /bin/bash
#________________________________________________________________________________________________
#
set -eu
#
export OMP_NUM_THREADS=8
#________________________________________________________________________________________________
#
basedir=$(pwd)

icon_grid_dir=/pool/data/ICON/grids/private/mpim/icon_preprocessing/source/grids

grids="icon_grid_0005_R02B04_G.nc"

year=2016

input_dir=/pool/data/ICON/grids/private/mpim/icon_preprocessing/source/extpar_input.$year

output_dir=./

# path to binaries

progdir=$basedir/bin

binary_topo=extpar_topo_to_buffer.exe
binary_tclim=extpar_cru_to_buffer.exe
binary_aot=extpar_aot_to_buffer.exe
binary_lu=extpar_landuse_to_buffer.exe
binary_ndvi=extpar_ndvi_to_buffer.exe
binary_soil=extpar_soil_to_buffer.exe
binary_flake=extpar_flake_to_buffer.exe
binary_alb=extpar_alb_to_buffer.exe

binary_extpar_consistency_check=extpar_consistency_check.exe

#________________________________________________________________________________________________
#
# path and filename settings

workdir=$basedir/work

if [[ ! -d ${workdir} ]]
then
  mkdir -p ${workdir} 
fi
cd ${workdir}

for grid in ${grids}
do
    input_file=${grid}
    output_file=${grid/icon/icon-extpar}

    #____
    #
    
    grib_output_filename=${output_file/nc/grb}
    netcdf_output_filename=$output_file
    grib_sample='GRIB2'
    
    echo $netcdf_output_filename
    echo $grib_output_filename
    
    #____________________________________________________________________________________________
    #
    # set target grid definition 
    #
    cat > INPUT_grid_org << EOF
&GRID_DEF 
 igrid_type          = 1
 domain_def_namelist = 'INPUT_ICON_GRID'
/ 
EOF
    #____
    #
    cat > INPUT_ICON_GRID << EOF
&icon_grid_info
  icon_grid_dir     = '${icon_grid_dir}'
  icon_grid_nc_file = '${input_file}'
/
EOF

    #____
    #
    
    raw_data_aot=aerosol_optical_thickness.nc
    buffer_aot=aot_buffer.nc
    output_aot=extpar_aot_icon.nc
    
    raw_data_tclim_coarse=absolute_hadcrut3.nc
    raw_data_tclim_fine=CRU_T2M_SURF_clim.nc
    buffer_tclim=crutemp_clim_extpar_buffer.nc
    output_tclim=extpar_crutemp_clim_icon.nc
    
    raw_data_glc2000=glc2000_byte.nc
    buffer_glc2000=extpar_landuse_buffer.nc
    output_glc2000=extpar_landuse_icon.nc
    
    raw_data_glcc=glcc_usgs_class_byte.nc
    buffer_glcc=glcc_landuse_buffer.nc
    output_glcc=extpar_glcc_landuse_icon.nc
    
    raw_data_globcover_0=GLOBCOVER_0_16bit.nc
    raw_data_globcover_1=GLOBCOVER_1_16bit.nc
    raw_data_globcover_2=GLOBCOVER_2_16bit.nc
    raw_data_globcover_3=GLOBCOVER_3_16bit.nc
    raw_data_globcover_4=GLOBCOVER_4_16bit.nc
    raw_data_globcover_5=GLOBCOVER_5_16bit.nc
    buffer_lu=extpar_landuse_buffer.nc
    output_lu=extpar_landuse_icon.nc

    case $year in
        2011)
            raw_data_globe_A10=GLOBE_A10g_corr.nc
            raw_data_globe_B10=GLOBE_B10g_corr.nc
            raw_data_globe_C10=GLOBE_C10g_corr.nc
            raw_data_globe_D10=GLOBE_D10g_corr.nc
            raw_data_globe_E10=GLOBE_E10g_corr.nc
            raw_data_globe_F10=GLOBE_F10g_corr.nc
            raw_data_globe_G10=GLOBE_G10g_corr.nc
            raw_data_globe_H10=GLOBE_H10g_corr.nc
            raw_data_globe_I10=GLOBE_I10g_corr.nc
            raw_data_globe_J10=GLOBE_J10g_corr.nc
            raw_data_globe_K10=GLOBE_K10g_corr.nc
            raw_data_globe_L10=GLOBE_L10g_corr.nc
            raw_data_globe_M10=GLOBE_M10g_corr.nc
            raw_data_globe_N10=GLOBE_N10g_corr.nc
            raw_data_globe_O10=GLOBE_O10g_corr.nc
            raw_data_globe_P10=GLOBE_P10g_corr.nc
            ;;
        2016)
            raw_data_globe_A10=ECMWF_A10_short.nc
            raw_data_globe_B10=ECMWF_B10_short.nc
            raw_data_globe_C10=ECMWF_C10_short.nc
            raw_data_globe_D10=ECMWF_D10_short.nc
            raw_data_globe_E10=ECMWF_E10_short.nc
            raw_data_globe_F10=ECMWF_F10_short.nc
            raw_data_globe_G10=ECMWF_G10_short.nc
            raw_data_globe_H10=ECMWF_H10_short.nc
            raw_data_globe_I10=ECMWF_I10_short.nc
            raw_data_globe_J10=ECMWF_J10_short.nc
            raw_data_globe_K10=ECMWF_K10_short.nc
            raw_data_globe_L10=ECMWF_L10_short.nc
            raw_data_globe_M10=ECMWF_M10_short.nc
            raw_data_globe_N10=ECMWF_N10_short.nc
            raw_data_globe_O10=ECMWF_O10_short.nc
            raw_data_globe_P10=ECMWF_P10_short.nc
            ;;
        *)
            echo "Input data only available from 2011 and 2016. Selected is $year ..."
            exit 1
            ;;
    esac
    buffer_globe=GLOBE_buffer.nc
    output_globe=extpar_GLOBE_icon.nc
    
    buffer_topo=topography_buffer.nc
    output_topo=extpar_topography_icon.nc
    
    raw_data_ndvi=NDVI_1998_2003.nc
    buffer_ndvi=NDVI_buffer.nc
    output_ndvi=extpar_NDVI_icon.nc
    
    raw_data_soil=FAO_DSMW_DP.nc
    buffer_soil=FAO_DSMW_buffer.nc
    output_soil=extpar_FAO_DSMW_icon.nc
    
    raw_data_flake=lakedepth.nc
    buffer_flake=flake_buffer.nc
    output_flake=extpar_flake_icon.nc
    
    raw_data_alb=month_alb.nc
    raw_data_alnid=month_alnid.nc
    raw_data_aluvd=month_aluvd.nc
    buffer_alb=alb_buffer.nc
    output_alb=extpar_alb_icon.nc
    
    # create input namelists 
    
    #____
    #
    cat > INPUT_AOT << EOF
&aerosol_raw_data
  raw_data_aot_path=''
  raw_data_aot_filename='${raw_data_aot}'
/  
&aerosol_io_extpar
  aot_buffer_file='${buffer_aot}'
  aot_output_file='${output_aot}'
/
EOF
    #____
    #
    cat > INPUT_TCLIM << EOF
&t_clim_raw_data
  raw_data_t_clim_path     = ''
  raw_data_t_clim_filename = '${raw_data_tclim_coarse}'
  raw_data_t_id            = 2
/
&t_clim_io_extpar
  t_clim_buffer_file       = '${buffer_tclim}'
  t_clim_output_file       = '${output_tclim}'
/
EOF
    #____
    #
    cat > INPUT_LU << EOF
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
EOF
    #____
    #
    cat > INPUT_ORO << EOF
&orography_io_extpar
  orography_buffer_file='${buffer_topo}',
  orography_output_file='${output_topo}'
/
&orography_raw_data
 itopo_type = 1
 lsso_param = .TRUE., 
 raw_data_orography_path='',
 ntiles_column = 4,
 ntiles_row = 4,
 topo_files = '${raw_data_globe_A10}', '${raw_data_globe_B10}', '${raw_data_globe_C10}',
              '${raw_data_globe_D10}', '${raw_data_globe_E10}', '${raw_data_globe_F10}',
              '${raw_data_globe_G10}', '${raw_data_globe_H10}', '${raw_data_globe_I10}',  
              '${raw_data_globe_J10}', '${raw_data_globe_K10}', '${raw_data_globe_L10}',  
              '${raw_data_globe_M10}', '${raw_data_globe_N10}', '${raw_data_globe_O10}',  
              '${raw_data_globe_P10}'  
/
EOF
    #____
    #
    cat > INPUT_OROSMOOTH << EOF
&orography_smoothing
  lfilter_oro=.TRUE.
/
EOF
    #____
    #
    cat > INPUT_RADTOPO << EOF
&radtopo
  lradtopo=.FALSE.,
  nhori=24,
/
EOF
    #____
    #
    cat > INPUT_NDVI << EOF
&ndvi_raw_data
  raw_data_ndvi_path='',
  raw_data_ndvi_filename='${raw_data_ndvi}'
/  
&ndvi_io_extpar
 ndvi_buffer_file='${buffer_ndvi}',
 ndvi_output_file='${output_ndvi}'
/
EOF
    #____
    #
    cat > INPUT_SOIL << EOF
&soil_raw_data
 isoil_data = 1,
 raw_data_soil_path='',
 raw_data_soil_filename='${raw_data_soil}'
/
&soil_io_extpar
  soil_buffer_file='${buffer_soil}',
  soil_output_file='${output_soil}'
/ 
EOF
    #____
    #
    cat > INPUT_FLAKE << EOF
&flake_raw_data
   raw_data_flake_path='',
   raw_data_flake_filename='${raw_data_flake}'
/
&flake_io_extpar
   flake_buffer_file='${buffer_flake}'
   flake_output_file='${output_flake}'
/
EOF
    #____
    #
    cat > INPUT_ALB << EOF
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
  alb_source='al'
  alnid_source='alnid'
  aluvd_source='aluvd'
/
EOF
    #____
    #
    # consistency check
    cat > INPUT_CHECK << EOF
&extpar_consistency_check_io
  grib_output_filename='${grib_output_filename}',
  netcdf_output_filename='${netcdf_output_filename}',
  orography_buffer_file='${buffer_topo}',
  soil_buffer_file='${buffer_soil}',
  lu_buffer_file='${buffer_lu}',
  glcc_buffer_file='${buffer_glcc}',
  flake_buffer_file='${buffer_flake}',
  ndvi_buffer_file='${buffer_ndvi}',
  t_clim_buffer_file='${buffer_tclim}',
  aot_buffer_file='${buffer_aot}',
  alb_buffer_file='${buffer_alb}',
  i_lsm_data=1,
  land_sea_mask_file='',
  number_special_points=0
/  
EOF
    #_____________________________________________________________________________________________
    #
    # link raw data files to local workdir
    
    #ln -sf ${raw_data_lu_path}/${raw_data_globcover}
    
    ln -sf ${input_dir}/${raw_data_aot}
    
    ln -sf ${input_dir}/${raw_data_tclim_coarse}
    ln -sf ${input_dir}/${raw_data_tclim_fine}
    
    ln -sf ${input_dir}/${raw_data_globcover_0}
    ln -sf ${input_dir}/${raw_data_globcover_1}
    ln -sf ${input_dir}/${raw_data_globcover_2}
    ln -sf ${input_dir}/${raw_data_globcover_3}
    ln -sf ${input_dir}/${raw_data_globcover_4}
    ln -sf ${input_dir}/${raw_data_globcover_5}
    
    ln -sf ${input_dir}/${raw_data_glc2000}
    ln -sf ${input_dir}/${raw_data_glcc}
    
    ln -sf ${input_dir}/${raw_data_globe_A10} 
    ln -sf ${input_dir}/${raw_data_globe_B10}
    ln -sf ${input_dir}/${raw_data_globe_C10} 
    ln -sf ${input_dir}/${raw_data_globe_D10} 
    ln -sf ${input_dir}/${raw_data_globe_E10}
    ln -sf ${input_dir}/${raw_data_globe_F10} 
    ln -sf ${input_dir}/${raw_data_globe_G10} 
    ln -sf ${input_dir}/${raw_data_globe_H10} 
    ln -sf ${input_dir}/${raw_data_globe_I10} 
    ln -sf ${input_dir}/${raw_data_globe_J10} 
    ln -sf ${input_dir}/${raw_data_globe_K10} 
    ln -sf ${input_dir}/${raw_data_globe_L10} 
    ln -sf ${input_dir}/${raw_data_globe_M10} 
    ln -sf ${input_dir}/${raw_data_globe_N10} 
    ln -sf ${input_dir}/${raw_data_globe_O10} 
    ln -sf ${input_dir}/${raw_data_globe_P10} 
    
    ln -sf ${input_dir}/${raw_data_ndvi}
    
    ln -sf ${input_dir}/${raw_data_soil}

    ln -sf ${input_dir}/${raw_data_flake}
    
    ln -sf ${input_dir}/${raw_data_alb}
    ln -sf ${input_dir}/${raw_data_alnid}
    ln -sf ${input_dir}/${raw_data_aluvd}
    
    #_____________________________________________________________________________________________
    #
    # run the programs:
    #
    # the next seven programs can run independent of each other
    
    set -o monitor
    trap add_next_job CHLD
    
    todo[0]=${progdir}/${binary_topo}
    todo[1]=${progdir}/${binary_alb} 
    todo[2]=${progdir}/${binary_aot}  
    todo[3]=${progdir}/${binary_tclim}
    todo[4]=${progdir}/${binary_lu}  
    todo[5]=${progdir}/${binary_ndvi}  
    todo[6]=${progdir}/${binary_soil}
    todo[7]=${progdir}/${binary_flake} 
    
    index=0
    max_jobs=1
   
    function add_next_job {
        if [[ $index -lt ${#todo[*]} ]]
        then
            ${todo[$index]} &
            (( index += 1 ))
        fi
    }
   
    while [[ $index -lt $max_jobs ]]
    do
        add_next_job
    done
    wait
    set +o monitor
    echo done

    # the consistency check requires the output of all previous processing step

    # ${progdir}/${binary_extpar_consistency_check}
    
    #_____________________________________________________________________________________________
    #
    # output renaming and mv to output directory    
    
    # if [[ ! -d ${output_dir} ]]
    # then
    #     mkdir -p ${output_dir} 
    # fi

    # processing_date=$(date +"%Y%m%d")
    # extpar=${grid/icon/icon_extpar}
    # extpar=${extpar/\.nc/_${processing_date}.nc}

    # mv $output_topo $output_dir/$extpar
    
    #_____________________________________________________________________________________________
    #
done




