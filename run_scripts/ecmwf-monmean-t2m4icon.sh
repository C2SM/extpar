#!/bin/ksh
module unload cdo
module load nco
module load cdo
module load grib_api

#
#  Interpolate ERA interim monthly means of SST and W_SNOW to ICON grid
#
export ICON_GRID_DIR=$2
 extpar_dir=$3   # working directory
 era_dir=$SCRATCH/extpar_data

integer retrieve=0
integer interpolate=1

# output file type GRIB or netCDF
#typeset oft ext H_name T_name
#oft=${1:-nc}

#case ${oft} in
#    g*|G*) oft=2
#        ext=g2
#        H_name='HSURF'
#        T_name='T_2M'
#        ;;
#    n*) oft=4
#        ext=nc
#        H_name='HSURF_eiCL'
#        T_name='T_2M_eiCL'
#        ;;
#     *) print -- "Unknown output file type ${oft}. Possible are GRIB2 or netCDF" >&2
#        exit 5
#        ;;
#esac

set -e
        ext=nc
        H_name='HSURF_eiCL'
        T_name='T_2M_eiCL'




typeset extpar_date=20171108
typeset validi_date=${extpar_date}00         # localValidityDate*

typeset year0=1986   # first year for climatology
typeset yearn=2015   # last year for climatology

if [[ ! -d ${extpar_dir} ]]; then
    mkdir -p ${extpar_dir}
fi
if [[ ! -d ${era_dir} ]]; then
    mkdir -p ${era_dir}
fi

typeset ifs_file era_oro
ifs_file="${era_dir}/ei_2t_an${year0}-${yearn}"
era_oro="${era_dir}/ei_oro_${year0}"

cd ${extpar_dir}

typeset modul_iconremap=/home/ms/de/dfr/routfox/abs/iconremap_mpi
integer year month day hour
typeset -Z2 mo

if (( retrieve == 1 )); then
#
# Retrieve monthly means of ERA Interim data from mars
#
    typeset dates

    year=${year0}
    while (( year <= ${yearn} )); do
        month=1
        while (( month <= 12 )); do
            mo=${month}
            dates="${dates}/${year}${mo}01"
            month=month+1
        done
        year=year+1
    done
    dates=${dates:1}

    mars <<EOF
retrieve, class=ei, stream=moda, type=an, expver=1,
  date=${dates},
  repress=gg, resol=av, grid=128, gaussian=regular,levtype=sfc, param=2t,
  target="${ifs_file}"
EOF
    mars <<EOF
retrieve, class=ei, stream=moda, type=an, expver=1,
  date=${year0}0101,
  repress=gg, resol=av, grid=128, gaussian=regular,levtype=sfc, param=z,
  target="${era_oro}.z"
EOF

#   Calculate mean over all years
    cdo ymonmean ${ifs_file} ${ifs_file}.tmp
    typeset climdate="centuryOfReferenceTimeOfData=12,yearOfCentury=11,day=11,hour=0"
    grib_set -s ${climdate} ${ifs_file}.tmp ${ifs_file}.mean
    rm ${ifs_file}.tmp

#   Convert from geopotential to height
    float ginv=1./9.80665  # to convert from geopotential to heights
    grib_set -s indicatorOfParameter=8,table2Version=2,scaleValuesBy=${ginv},${climdate} \
        ${era_oro}.z ${era_oro}

fi

if (( interpolate == 1 )); then
#
#  Interpolate from IFS to ICON grid
#
    typeset id='eiT2icon'
    integer nt=20   # number of parallel threads
    typeset intp_method
#   intp_method=2            # area-weighted formula
#    intp_method=3            # RBF scalar
    intp_method=4            # nearest-neighbor scalar
#   intp_method=8            # barycentric scalar

    typeset remap_nml=${id}.nml$$

    set -A locCreateDate $(datconv -C $(date +%Y%m%d%H))
    set -A valdate $(datconv -C ${validi_date})
    integer gridNr
#
#  Write rule files for grib_filter
#

#   ICON_GRID_DIR and ICON_XML_GRID_TABLE are needed for icon_grid_get
    if [[ -z "${ICON_GRID_DIR}" ]]; then
        export ICON_GRID_DIR=$PERM/grids
    fi
    if [[ -z "${ICON_XML_GRID_TABLE}" ]]; then
        export ICON_XML_GRID_TABLE=/hpc/rhome/routfor/routfox/icon/xml/dwd_grids.xml
    fi

#   iconremap sets one date for all GRIB records. Therefore, the interpolation must
#   be done separately for each month. Split up the original file
    grib_copy ${ifs_file}.mean ${id}.[month]

#    for gridNr in 23 24 25 26 27 28
#   do
# 0026_R03B07_G
    icon_res=$1 
    icon_grid_file=icon_grid_${icon_res}.nc
        if [[ ${icon_grid_file} == *N02.nc ]]; then
            gpid=2
        else
            gpid=1
        fi
    icon_grid=${ICON_GRID_DIR}/${icon_grid_file}


#        icon_grid=$(icon_grid_get -q ${gridNr} 78 255)
#        if [[ ${icon_grid} == *N02.nc ]]; then
#            gpid=2
#        else
#            gpid=1
#        fi
        ires=$(basename ${icon_grid})
        ires=${ires/icon_grid_/}
        ires=${ires/.nc/}
        ncstorage_file=eiT2icon${ires}_storage.nc
        icon_file=$(basename ${ifs_file})_${ires}
        rm -f ${icon_file}

#       iconremap sets one date for all GRIB records. Therefore, the interpolation must
#       be done separately for each month. Month 0 is the orography
        for month in $(seq 0 12)
        do
            mo=${month}

            if (( month == 0 )); then
#             interpolate orography
              year=1
              day=1
              hour=0
              cat > ${remap_nml} <<-REMAP_EOF0
&remap_nml
 in_grid_filename="${era_oro}",
 in_filename="${era_oro}",
	REMAP_EOF0
            else
              year=1111
              day=11
              hour=11
              cat > ${remap_nml} <<-REMAP_EOF0
&remap_nml
 in_grid_filename="${id}.${month}",
 in_filename="${id}.${month}",
	REMAP_EOF0
            fi

            cat >> ${remap_nml} <<-REMAP_EOF
 out_filename="${icon_file}.${mo}",
 out_grid_filename="${icon_grid}",
 in_type  = 1,
 out_type = 2,                    ! ICON triangular grid
 out_filetype=${oft},             ! output format GRIB2 or netCDF
! ncstorage_file="${ncstorage_file}"
 extra_grib_keys_int = "centre",78, "subCentre",255, 
                       "year",${year}, "day",${day}, "hour",${hour}, "minute",0,
                       "grib2LocalSectionPresent",1, "localDefinitionNumber",254, 
                       "localCreationDateYear", ${locCreateDate[1]}, "localValidityDateYear", ${valdate[1]},
                       "localCreationDateMonth",${locCreateDate[2]}, "localValidityDateMonth",${valdate[2]},
                       "localCreationDateDay",  ${locCreateDate[3]}, "localValidityDateDay",  ${valdate[3]},
                       "localCreationDateHour", ${locCreateDate[4]}, "localValidityDateHour", ${valdate[4]},
                       "shapeOfTheEarth",6, "typeOfGeneratingProcess",9, "backgroundProcess",0,
                       "generatingProcessIdentifier", ${gpid},
                       "significanceOfReferenceTime",0, "typeOfProcessedData",0,
                       "bitsPerValue",16,
	REMAP_EOF

            if (( month == 0 )); then
              cat >> ${remap_nml} <<-REMAP_FIELD
              "typeOfGeneratingProcess",9,"productDefinitionTemplateNumber",8,"typeOfStatisticalProcessing",0
/
&input_field_nml
!inputname='Z',
!code = 129,
 inputname='HSURF',
 code = 8,
 outputname='${H_name}',
 intp_method=${intp_method}
/
	REMAP_FIELD
            else
              cat >> ${remap_nml} <<-REMAP_FIELD
             "typeOfFirstFixedSurface",103,"scaledValueOfFirstFixedSurface",2,
/
&input_field_nml
 inputname='2T',
 outputname='${T_name}',
 code = 167,
 intp_method=${intp_method}
/
	REMAP_FIELD
            fi

            cat ${remap_nml}

            aprun -n $EC_total_tasks  -N $EC_tasks_per_node \
              -j $EC_hyperthreads -d $EC_threads_per_task \
            ${modul_iconremap} -vv --remap_nml=${remap_nml}
            ls -l ${icon_file}.${mo}
            cat nml.log
#           rm -f nml.log

        done

        if [[ ${ext} == 'g2' ]]; then
            cat ${icon_file}.[01][0-9] > ${icon_file}.${ext}
        elif [[ ${ext} == 'nc' ]]; then
#           paste T_2M_eiCL together
            cdo copy ${icon_file}.0[1-9] ${icon_file}.1[0-2] ${icon_file}.${ext}
#           remove time dimension from orography netCDF file
            typeset h_tmp=h.nc$$
        /hpc/rhome/software/nco/4.6.2/bin/ncwa -a time ${icon_file}.00 ${h_tmp}
#           append orography file to temperature file
        /hpc/rhome/software/nco/4.6.2/bin/ncks -C -A -v HSURF_eiCL ${h_tmp} ${icon_file}.${ext}
            rm ${h_tmp}
        else
            print -- "Unknown extension ${ext}! You should not get here!" >& 2
            exit 6
        fi
#        rm ${icon_file}.[01][0-9]
#        ls -l ${icon_file}.${ext}


#    rm ${remap_nml}
#    rm ${id}.*

fi
#Beautify the resulting netcdf to use it in EXTPAR consistency check that expects a special buffer format time,ie,je,ke
ncrename -d ncells,ie ${icon_file}.nc ${icon_file}_ie.nc
#add second dimension je -> Caution! ordering of dimensions in opposite in NetCDF and Fortran
ncap2 -s 'defdim("je",1);T_2M_eiCL_ieje[$time,$je,$ie]=T_2M_eiCL' -s 'HSURF_eiCL_ieje[$je,$ie]=HSURF_eiCL' -O ${icon_file}_ie.nc ${icon_file}_ieje.nc
#add third dimension ke
ncap2 -s 'defdim("ke",1);T_2M_eiCL_iejeke[$time,$ke,$je,$ie]=T_2M_eiCL_ieje' -s 'HSURF_eiCL_iejeke[$ke,$je,$ie]=HSURF_eiCL_ieje' -O ${icon_file}_ieje.nc ${icon_file}_iejeke.nc
# delete old variables
ncks -C -O -x -v T_2M_eiCL,T_2M_eiCL_ieje,HSURF_eiCL,HSURF_eiCL_ieje  ${icon_file}_iejeke.nc ${icon_file}_iejeke_tmp.nc
# rename variables back
ncrename  -O -v T_2M_eiCL_iejeke,T_2M_CLIM -v HSURF_eiCL_iejeke,TOPO_CLIM ${icon_file}_iejeke_tmp.nc ${icon_file}_BUFFER.nc
rm -rf ${icon_file}_ie.nc ${icon_file}_ieje.nc ${icon_file}_iejeke.nc ${icon_file}_iejeke_tmp.nc
cdo -f copy grib2  ${icon_file}_BUFFER.nc ${icon_file}_BUFFER.g2

