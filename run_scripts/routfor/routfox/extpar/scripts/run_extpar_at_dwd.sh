#!/bin/ksh

# ---  run_extpar_at_dwd.sh --- #

########################################################################
#                                                                      #
#  generate EXTernal PARameters for ICON                               #
#                                                                      #
#  call:         run_extpar_at_dwd.sh            \                     #
#                  -dir_grid=${dir_grid}         \                     #
#                  -info_grid=${info_grid}                             #
#                                                                      #
#                                                                      #
#  parameters:   -dir_grid   - directory of the ICON grid file         #
#                              default is  ROUTINE_DIR_IGLO_EXT        #
#  routine_config | fgrep -i ROUTINE_DIR_IGLO_EXT                      #
#                -info_grid  - ICON grid base (e.g. 0026_R03B07_G)     #
#                -work_dir   - working directory (optional)            #
#  -m Use MERIT/REMA orography (default)                               #
#  -x NoArt-Flag: hwsdART is processed                                 #
#  -q Queue-Version: Use BIG class with suffix _rh8                    #
#                                                                      #
#  example: run_extpar_at_dwd.sh -dir_grid=${filename}          \      #
#                                 -info_grid=${info_grid} -m -x -q     #
#                                                                      #
#  modified by:  J. Helmert from extpar_cosmo_routi.sh by Ch. Koziar   #
#                                                                      #
#  last change: 2019/06/13  by  J. Helmert                             #
#                                                                      #
########################################################################
#                                                                      #
#  Modifications:                                                      #
#  --------------                                                      #
#                                                                      #
#  2019/06/13  -  First version                                        #
#  2019/10/14  -  Add notes for stand-alone version in user env        #
#  2021/05/25  -  Adaptations for extpar V5.5                          #
#  2024/03/01  -  Add switches -x for NoArt and -q for queue version   #
#                                                                      #
########################################################################

# To create an own copy of the EXTPAR environment:
# 1. copy the content of ~routfor/routfox/extpar to an own directory, e.g. ${HOME}/OWN_EXTPAR
# 2. Adapt below in the script EXTPAR_HOME=${HOME}/OWN_EXTPAR
# 3. For using own binaries of EXTPAR point below in the script PROGDIR_LC to the location of the EXTPAR binaries
# 4. Ensure that the correct names of the binaries are used, i.e.
#       modify the content of FileNames in the ADM workbench (remove prefix ->tst-<, add suffix ->.new<-):
# BINARIES    = $(STDROOT)/extpar_aot_to_buffer.new \
#               $(STDROOT)/extpar_consistency_check.new \
#               $(STDROOT)/extpar_cru_to_buffer.new \
#               $(STDROOT)/extpar_flake_to_buffer.new \
#               $(STDROOT)/extpar_landuse_to_buffer.new \
#               $(STDROOT)/extpar_ndvi_to_buffer.new \
#               $(STDROOT)/extpar_soil_to_buffer.new \
#               $(STDROOT)/extpar_alb_to_buffer.new \
#               $(STDROOT)/extpar_topo_to_buffer.new \
#               $(STDROOT)/extpar_ahf_to_buffer.new \
#               $(STDROOT)/extpar_isa_to_buffer.new \
#               $(STDROOT)/extpar_emiss_to_buffer.new \
#               $(STDROOT)/extpar_sgsl_to_buffer.new

set -x

# exit on error:
set -e

# ----------------------------------------------------------- #
# ---  modify/add: PATH, FPATH, ROUTINE_PERL, path_batch  --- #
# ----------------------------------------------------------- #
 ( . inc.define_path 2>/dev/null ) \
    && . inc.define_path           \
    || . ~routfor/routfox/bin/inc.define_path

# only until Thomas has synchronized his bins:
#export PATH=~for1han/bin:~for1han/ebin:$PATH

module list
module clear <<< 'y'
module list

error_count=0
run_command()
{
    set +e
    echo ">> Run ${1%% *} ..."
    start=$(date +%s.%N)
    eval qsubw  $1 > ${logfile} 2>&1
    rc=$?
    set -e
    printf "   Return code: %i\n" $rc
    end=$(date +%s.%N)
    runtime=$(bc -l <<< "$end - $start")
    if (( rc > 0 ))
    then
        (( error_count++ ))
    fi
    case $rc in
        0)
            echo "   SUCCESS ${1%% *}"
            echo "   LOGFILE is $2"
            ;;
        127)
            echo "   ERROR ${1%% *}: command not found"
            echo "   LOGFILE is $2"
            ;;
        *)
            echo "   ERROR ${1%% *}: fatal error - return code $rc"
            echo "   LOGFILE is $2"
            ;;
    esac
    echo "   execution time: $runtime s"
    [[ -f ${logfile} ]] && cat ${logfile}
    [[ -f $2 ]] && cat $2

    # @ckoziar, 20190926:
    if (( rc != 0 )) ; then
        exit ${rc}
    fi
}

scriptpath=$0
scriptname=${scriptpath##*/}
logfile=${scriptname%.*}_$(date +%Y%m%d%H%M%S).log


eval $(routine_config)

# ------------------------------------------------ #
# ---  help variables to print error masseges  --- #
# ------------------------------------------------ #
typeset          name=${0##*/} ; name=${name%.ne[uw]}
typeset          cn="${name}:"
typeset -L${#cn} bn=''
typeset -u       NAME=${name}

# ------------------------ #
# ---  usage function  --- #
# ------------------------ #
function usage { print "usage: ${name} -h" ; }

# ---------------------------------------------- #
# ---  parameter test - simple               --- #
# ---  (OPTIND: next parameter to evaluate)  --- #
# ---------------------------------------------- #
typeset    argstring="$@"

integer    errflag=0
typeset    errparm=''
typeset    ECHO=''

typeset    xopt=''

typeset    fflag=''
typeset    Dflag=''

typeset    logdir=${TMPDIR}

if (( errflag == 1 )) ; then
  print "${cn} parameter error (${errparm})  -  usage: ${name} -h" >&2
  exit 201
fi

# ---------------------------------------------- #
# ---  parameter test - extended             --- #
# ---  (OPTIND: next parameter to evaluate)  --- #
# ---------------------------------------------- #
typeset    argstring="$@"
typeset    arg_array
set -A     arg_array -- "$@"

integer    errflag=0
typeset    errparm=''
typeset    ECHO=''

typeset    dir_grid=$ROUTINE_DIR_IGLO_EXT
typeset    gpid=''
typeset    info_grid=''
typeset    work_dir=''
typeset    nml_dir=''
typeset    queue_version='rc_big'
typeset    globeFlag=''
typeset    meritFlag='true'
# @ckoziar, 20200326: enable emissFlag for temporary tests only !!!
# typeset    emissFlag=''
typeset    emissFlag='true'
typeset    NoArtFlag=''

while getopts d:exqghi:mn:w: option ; do
  case ${option} in
     d) # dir_grid
        optarg=${OPTARG}
        case ${optarg%%=*} in
          ir_grid) if [[ "${arg_array[$((${OPTIND}-2))]%%=*}" == "-${option}${optarg%%=*}" ]] ; then
                       if [[ ":${optarg%%=*}" == ":${optarg}" ]] ; then
                         dir_grid="${arg_array[$((${OPTIND}-1))]}"
                         OPTIND=$((${OPTIND} + 1))
                       else
                         dir_grid="${OPTARG#*=}"
                       fi
                     fi
                     ;;
           *)        errflag=1
                     ;;
        esac
        ;;
     e) emissFlag='true'
        ;;
     x) NoArtFlag='false'
        ;;
     q) queue_version='rc_big_rh8'
        ;;    
     g) optarg=${OPTARG}
        case ${optarg%%=*} in
          pid) # generating process identifier:
               if [[ "${arg_array[$((${OPTIND}-2))]%%=*}" == "-${option}${optarg%%=*}" ]] ; then
                  if [[ ":${optarg%%=*}" == ":${optarg}" ]] ; then
                    gpid="${arg_array[$((${OPTIND}-1))]}"
                    OPTIND=$((${OPTIND} + 1))
                  else
                    gpid="${OPTARG#*=}"
                  fi
               fi
               ;;
          *)   # disable-ASTER / enable-GLOBE flag
               globeFlag='true'
               ;;
        esac
        ;;
     h) scr_info $0
        exit 0
        ;;
     i) # info_grid
        optarg=${OPTARG}
        case ${optarg%%=*} in
          nfo_grid) if [[ "${arg_array[$((${OPTIND}-2))]%%=*}" == "-${option}${optarg%%=*}" ]] ; then
                       if [[ ":${optarg%%=*}" == ":${optarg}" ]] ; then
                         info_grid="${arg_array[$((${OPTIND}-1))]}"
                         OPTIND=$((${OPTIND} + 1))
                       else
                         info_grid="${OPTARG#*=}"
                       fi
                    fi
                    ;;
          *)        errflag=1
                    ;;
        esac
        ;;
     m) meritFlag='true'
        ;;
     n) # namelist directory
        optarg=${OPTARG}
        case ${optarg%%=*} in
          ml_dir) if [[ "${arg_array[$((${OPTIND}-2))]%%=*}" == "-${option}${optarg%%=*}" ]] ; then
                     if [[ ":${optarg%%=*}" == ":${optarg}" ]] ; then
                       nml_dir="${arg_array[$((${OPTIND}-1))]}"
                       OPTIND=$((${OPTIND} + 1))
                     else
                       nml_dir="${OPTARG#*=}"
                     fi
                  fi
                  ;;
          *)      errflag=1
                   ;;
        esac
        ;;
     w) # work_dir
        optarg=${OPTARG}
        case ${optarg%%=*} in
          ork_dir) if [[ "${arg_array[$((${OPTIND}-2))]%%=*}" == "-${option}${optarg%%=*}" ]] ; then
                       if [[ ":${optarg%%=*}" == ":${optarg}" ]] ; then
                         work_dir="${arg_array[$((${OPTIND}-1))]}"
                         OPTIND=$((${OPTIND} + 1))
                       else
                         work_dir="${OPTARG#*=}"
                       fi
                     fi
                     ;;
           *)        errflag=1
                     ;;
        esac
        ;;
    \?) errflag=1
        ;;
  esac
  if (( errflag == 1 )) ; then
    eval errparm=\"\${$((${OPTIND}-1))}\"
    break
  fi
done
shift $((${OPTIND} - 1))

if [[ -z "${info_grid}" ]] ; then
  print "${cn} parameter error - usage: ${name} -h" >&2
  exit 201
fi

if [[ "${dir_grid}" == "${ROUTINE_DIR_IGLO_EXT}" && \
      ! -f ${dir_grid}/${info_grid} ]] ; then
   info_grid=icon_grid_${info_grid}.nc
fi

if [[ ! -f ${dir_grid}/${info_grid} ]] ; then
  print "${cn} parameter error - no corresponding grid file found (${dir_grid}/${info_grid})." >&2
  exit 201
fi

if [[ -n "${nml_dir}" ]] ; then
  if [[ ! -d ${nml_dir} ]] ; then
    print "${cn} parameter error - no valid directory specified by argument '-nml_dir' ('${nml_dir}')." >&2
    exit 201
  fi
fi

if [[ -n "${meritFlag}" ]] ; then
  # disable-ASTER / enable-GLOBE flag in case of MERIT:
  globeFlag='true'
fi

logdir=${work_dir:-${TMPDIR}}

###################################################################################################


# @ckoziar, 20190806:
# TODAY=$(date  +%Y%m%d)
TODAY=$(date '+%s')

ICON_GRID=${info_grid}
# PROGDIR_LC=~jhelmert/EXTPAR_GIT/release5.0/bin
# PROGDIR_LC=/usr/local/pkg/for0adm/abs/
# PROGDIR_LC=~routfor/routfox/patch-rc
# PROGDIR_LC=~routfor/routfox/patch-rc:/hpc/rhome/for0adm/x86/abs
# PROGDIR_LC=~routfor/routfox/patch-rc:/hpc/rhome/for0adm/X86/abs
# PROGDIR_LC=~routfor/routfox/patch-rc:/hpc/rhome/for0adm/nwp/x86/util/bin
PROGDIR_LC=~routfor/routfox/patch-rc:/hpc/rhome/for0adm/nwp/x86/util/bin/
GRIDDIR=${dir_grid}
# EXTPAR_HOME=~routfor/routfox/extpar
 EXTPAR_HOME=/hpc/rhome/for0adm/nwp/x86/util/run_scripts/routfor/routfox/extpar
# EXTPAR_HOME=/hpc/uhome/jhelmert/runscript
TMP=${TMPDIR}
# @ckoziar, 20190806:
# WORKDIR=$TMP/${ICON_GRID}_${TODAY}
WORKDIR=${work_dir:-${TMP}/${ICON_GRID}_${TODAY}_$$}
print ">>>  Settings:"
if [[ -n "${globeFlag}" ]] ; then print ">>>  use_GLOBE=${globeFlag}" ; fi
if [[ -n "${emissFlag}" ]] ; then print ">>>  emissFlag=${emissFlag}" ; fi
if [[ -n "${NoArtFlag}" ]] ; then print ">>>  NoArtFlag=${NoArtFlag}" ; fi
if [[ -n "${meritFlag}" ]] ; then print ">>>  meritFlag=${meritFlag}" ; fi
print ">>>  dir_grid=${dir_grid}"
print ">>>  info_grid=${info_grid}"
print ">>>  progdir_lc=${PROGDIR_LC}"
print ">>>  extpar_home=${EXTPAR_HOME}"
print ">>>  tmp=${TMP}"
print ">>>  workdir=${WORKDIR}"

if (( errflag == 1 )) ; then
  print "${cn} parameter error (${errparm})  -  usage: ${name} -h" >&2
  exit 201
fi

# @ckoziar, 20210210: uncomment following 2 lines:
# ncks='/hpc/sw/nco/4.8.1/bin/ncks'
# ncwa='/hpc/sw/nco/4.8.1/bin/ncwa'

echo "cd $TMP"
cd $TMP

# @jhelmert/ckoziar, 20191022:
# check for minimum and maximum latitude in case of ASTER orography:
if [[ -s ${GRIDDIR}/${ICON_GRID} ]] ; then
  meshsize_km=$(ncks -M ${GRIDDIR}/${ICON_GRID} | \
                awk '/grid_root/  { GRIDROOT=$(NF-1)  }
                     /grid_level/ { GRIDLEVEL=$(NF-1) }
                     END { print(5050./(GRIDROOT*2**GRIDLEVEL))}
                    ' -
               )

  # @fprill/ckoziar, 20191204: cdo infov doesn't work in every case --> replace by ncwa/ncks
# min_lat_vert=$(cdo infov ${GRIDDIR}/${ICON_GRID} | \
#                grep -i latitude_vertices | awk '{ print $9 * 180. / atan2(0, -1) }')
  ncwa -y min -v latitude_vertices -O ${GRIDDIR}/${ICON_GRID} foo.$$.nc && \
  min_lat_vert=$(ncks -H -s "%f" foo.$$.nc | awk '/^[-+0-9.,]+$/ { print $1 * 180. / atan2(0, -1) }')

# max_lat_vert=$(cdo infov ${GRIDDIR}/${ICON_GRID} | \
#                grep -i latitude_vertices | awk '{ print $11 * 180. / atan2(0, -1) }')
  ncwa -y max -v latitude_vertices -O ${GRIDDIR}/${ICON_GRID} foo.$$.nc && \
  max_lat_vert=$(ncks -H -s "%f" foo.$$.nc | awk '/^[-+0-9.,]+$/ { print $1 * 180. / atan2(0, -1) }')

  if [[ "${meshsize_km}"  != +([-+0-9.,]) ]] ; then
    print -- "ERROR: Cannot determine horizontal resolution of target grid file!\n" >&2
    exit 202
  fi

  if [[ "${min_lat_vert}" != +([-+0-9.,]) || \
        "${max_lat_vert}" != +([-+0-9.,])    \
     ]] ; then
    print -- "ERROR: Cannot determine minimum/maximum latitude of target grid file!\n" >&2
    exit 203
  fi

  rm -f foo.$$.nc

  grid_number=$(ncks -M ${GRIDDIR}/${ICON_GRID} | \
                awk '/number_of_grid_used/ { print $(NF-1) }')

  if [[ -z "${grid_number}" && -w ${GRIDDIR}/${ICON_GRID} ]] ; then
    print -- "try to set grid number of '${GRIDDIR}/${ICON_GRID}' to 9999 now..."
    ncatted -h -a number_of_grid_used,global,o,l,9999 ${GRIDDIR}/${ICON_GRID} && grid_number=9999
  fi

  print -- "mesh size of '${GRIDDIR}/${ICON_GRID}': ${meshsize_km} km"
  print -- "grid number of '${GRIDDIR}/${ICON_GRID}': '${grid_number}'"
  print -- "minimum latitude of vertices in '${GRIDDIR}/${ICON_GRID}': ${min_lat_vert} deg"
  print -- "maximum latitude of vertices in '${GRIDDIR}/${ICON_GRID}': ${max_lat_vert} deg"

# @ckoziar, 20220711: bugfix (globeFlag AND meritFlag have to be empty!)
# if [[ -z "${globeFlag}" || -z "${meritFlag}" ]] ; then
# if [[ -z "${globeFlag}" && -z "${meritFlag}" ]] ; then
      
   if (( meshsize_km < 3.0 )) ; then
      # switch from GLOBE to ASTER orography for resolutions higher than 3.0 km:
      print -- " Try to switch to ASTER orography because target horizontal resolution is higher than 3.0 km"

       if (( max_lat_vert < 60.0 && min_lat_vert > -60.0 )) ; then
        echo "Switch to ASTER orography successful because target domain is within 60degN and 60degS"
        meritFlag=''
        globeFlag=''
       fi
#  fi
  fi
fi

# cat > igrid <<EOF
# ${ICON_GRID}
# EOF
# cat > igrid_dir <<EOF
# ${GRIDDIR}
# EOF
# cat > iwork_dir <<EOF
# ${WORKDIR}
# EOF
# cat > progdir_lc <<EOF
# ${PROGDIR_LC}
# EOF
# cat > progdir_xc <<EOF
# ${PROGDIR_XC}
# EOF

###############################################################################################
echo "Step 1/3 for grid ${ICON_GRID} ... Prepare EXTPAR environment in $WORKDIR"
sed -e "s:%igrid%:${ICON_GRID}:g
        s:%igrid_dir%:${GRIDDIR}:g
        s:%iwork_dir%:${WORKDIR}:g
        s@%progdir_lc%@${PROGDIR_LC}@g
        s:%queue_version%:${queue_version}:g
        s:%loggingfile%:${logdir}/o.extpar_icon.extpar_prep.$$.log:g
       " ${EXTPAR_HOME}/job/extpar_prep.sh > extpar_prep.$$.sh && \
run_command extpar_prep.$$.sh ${logdir}/o.extpar_icon.extpar_prep.$$.log
echo "Step 1/3 for grid ${ICON_GRID} ... Prepare EXTPAR environment in $WORKDIR ... finished"
###############################################################################################




###############################################################################################
logfile4zip=external_parameter_icon_${ICON_GRID%.nc}_tiles.extpar.log
echo "Step 2/3 for grid ${ICON_GRID} ... Start EXTPAR modules in $WORKDIR"
sed -e "s:%igrid%:${ICON_GRID}:g
        s:%igrid_dir%:${GRIDDIR}:g
        s:%iwork_dir%:${WORKDIR}:g
        s@%progdir_lc%@${PROGDIR_LC}@g
        s:%nml_dir%:${nml_dir}:g
        s:%queue_version%:${queue_version}:g
        s:%globeFlag%:${globeFlag}:g
        s:%meritFlag%:${meritFlag}:g
        s:%emissFlag%:${emissFlag}:g
        s:%NoArtFlag%:${NoArtFlag}:g
        s:%loggingfile%:${logdir}/o.extpar_icon.extpar_modules.$$.log:g
        s:%logfile4zip%:${logdir}/${logfile4zip}:g
       " ${EXTPAR_HOME}/job/extpar_modules.sh > extpar_modules.$$.sh && \
run_command extpar_modules.$$.sh ${logdir}/o.extpar_icon.extpar_modules.$$.log
echo "Step 2/3 for grid ${ICON_GRID} ... EXTPAR modules ... finished"
###############################################################################################


###############################################################################################
echo "Step 3/3 for grid ${ICON_GRID} ... Checking EXTPAR file and converting to GRIB2"
[[ -s ${WORKDIR}/external_parameter_icon_${ICON_GRID%.nc}_tiles.nc ]] || echo "EXTPAR file is external_parameter_icon_${ICON_GRID%.nc}_tiles.nc in ${WORKDIR}"
if [[ ! -s ${WORKDIR}/external_parameter_icon_${ICON_GRID%.nc}_tiles.nc ]]; then
   icon_grid=${ICON_GRID%.nc} ; icon_grid=${icon_grid#icon_grid_} ;
   if [[ ! -s ${WORKDIR}/external_parameter_icon_${icon_grid}_tiles.nc ]]; then
      echo "$0 error: ${WORKDIR}/external_parameter_icon_${icon_grid}_tiles.nc is missing, something was going wrong before ..."
      exit
   fi
fi
sed -e "s:%igrid%:${ICON_GRID}:g
        s:%iwork_dir%:${WORKDIR}:g
        s:%queue_version%:${queue_version}:g
        s:%loggingfile%:${logdir}/o.extpar_icon.extpar2grib2.$$.log:g
       " ${EXTPAR_HOME}/job/extpar2grib2.sh > extpar2grib2.$$.sh && \
run_command extpar2grib2.$$.sh ${logdir}/o.extpar_icon.extpar2grib2.$$.log
[[ -s ${WORKDIR}/external_parameter_icon_${ICON_GRID%.nc}_tiles.g2 ]] || echo "EXTPAR file in GRIB2 is external_parameter_icon_${ICON_GRID%.nc}_tiles.g2 in ${WORKDIR}"
if [[ ! -s ${WORKDIR}/external_parameter_icon_${ICON_GRID%.nc}_tiles.g2 ]]; then
   icon_grid=${ICON_GRID%.nc} ; icon_grid=${icon_grid#icon_grid_} ;
   if [[ ! -s ${WORKDIR}/external_parameter_icon_${icon_grid}_tiles.g2 ]]; then
      echo "$0 error: ${WORKDIR}/external_parameter_icon_${icon_grid}_tiles.g2 is missing, something was going wrong before ..."
      exit
   fi
fi
echo "Step 3/3 for grid ${ICON_GRID} ... Checking EXTPAR file and converting to GRIB2 ... finished"
###############################################################################################

# revise generating process identifier:
if [[ -z "${gpid}" ]] ; then
  if [[ "${ICON_GRID%.nc}" == icon_grid_[0-9][0-9][0-9][0-9]_R[0-9][0-9]B[0-9][0-9]_L ]] ; then
    gpid=11
  elif [[ "${ICON_GRID%.nc}" == icon_grid_[0-9][0-9][0-9][0-9]_R[0-9][0-9]B[0-9][0-9]_LN+([0-9]) ]] ; then
    gpid=12
  fi
fi

if [[ -n "${gpid}" ]] ; then
  cd ${WORKDIR}
  ls -l external_parameter_icon_${ICON_GRID%.nc}_tiles.g2
  echo "grib_set -v -s generatingProcessIdentifier=${gpid} external_parameter_icon_${ICON_GRID%.nc}_tiles.g2 external_parameter_icon_${ICON_GRID%.nc}_tiles.g2.tmp"
  grib_set -v -s generatingProcessIdentifier=${gpid} external_parameter_icon_${ICON_GRID%.nc}_tiles.g2 external_parameter_icon_${ICON_GRID%.nc}_tiles.g2.$$.tmp
  echo "mv external_parameter_icon_${ICON_GRID%.nc}_tiles.g2.tmp external_parameter_icon_${ICON_GRID%.nc}_tiles.g2"
  mv external_parameter_icon_${ICON_GRID%.nc}_tiles.g2.$$.tmp external_parameter_icon_${ICON_GRID%.nc}_tiles.g2
  ls -l external_parameter_icon_${ICON_GRID%.nc}_tiles.g2
  cd -
fi
