#!/bin/ksh -l

#PBS -w NAME=extpar_modules
#PBS -N ${NAME}
#PBS -S /bin/ksh
#PBS -q %queue_version%
#PBS -m n
#PBS -r n
#PBS -l cpunum_job=32
##PBS -l memsz_prc=192gb
#PBS -l vmemsz_prc=192gb
#PBS -l memsz_job=384gb
#PBS -l vmemsz_job=384gb
#PBS -l elapstim_req=172800
#PBS -o %loggingfile%
#PBS -j o

set -x

# exit on error:
set -e

set +e
type module 2>&1 >/dev/null
rcmod=$?
set -e
if (( rcmod != 0 )) ; then
  . /usr/share/Modules/init/ksh
fi

module list
module clear <<< 'y'
module list

module load cdo
# @ckoziar, 20220518: uncommented
#module load python
export PYTHONPATH="Hurra"

module list

ulimit -s unlimited
ulimit -c 0

#________________________________________________________________________________
error_count=0
run_command()
{
    progdir="${PATH}" # :default
    argno=$(echo "$*" | awk '{ print NF-1 }')
    if (( argno == 0 )) ; then
      # Fortran binary (arguments overgiven implicit by namelist):
      progdir="${PROGDIR}"
    elif (( argno > 0 )) ; then
      # Shell Skript (arguments overgiven explicit by skript parameters):
      progdir="${EXTPAR_HOME}/job"
    fi
    modul=$(echo "$*" | awk '{ print $1 }')
    args=$(echo "$*" | awk '{ sep=""; for(i=2; i<=NF; i++) { printf "%s%s", sep, $i; sep=" " } }')
    set +e
    modul=$(PATH=${progdir} whence ${modul})
    rc=$?
    set -e
    if (( rc != 0 )) ; then
      modul=modul_not_set
    fi
    if [[ ! -e ${modul} || ! -x ${modul} ]] ; then
      print -- "   ERROR: command '${modul}' not found or not executable."
      exit 201
    fi
    echo ">> Run $* ..."
    start=$(date +%s.%N)
    set +e
#   eval ./$* > ${logfile} 2>&1
    eval ${modul} ${args} > ${logfile} 2>&1
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
            ;;
        127)
            echo "   ERROR ${1%% *}: command not found"
            ;;
        130)
            echo "   ERROR ${1%% *}: script terminated by Ctrl-C"
            ;;
        *)
            echo "   ERROR ${1%% *}: fatal error - return code $rc"
            ;;
    esac
    echo "   execution time: $runtime s"
    [[ -f ${logfile} ]] && ls -l ${logfile} && cat ${logfile}
    [[ -f ${logfile} ]] && cat ${logfile} >> ${logfileforzip}

    # @ckoziar, 20210510: special logfiles of extpar v5.5
    logfile55=${modul##*/} ; logfile55=${logfile55%.exe}.log
    [[ -f ${logfile55} ]] && ls -l ${logfile55} && cat ${logfile55}
    [[ -f ${logfile55} ]] && cat ${logfile55} >> ${logfileforzip}

    # @ckoziar, 20190926:
    if (( rc != 0 )) ; then
        exit ${rc}
    fi
}

run_python()
{
    progdir="${PATH}" # :default
    argno=$(echo "$*" | awk '{ print NF-1 }')
    if (( argno == 0 )) ; then
      # Fortran binary (arguments overgiven implicit by namelist):
      progdir="${WORKDIR}"
    elif (( argno > 0 )) ; then
      # Shell Skript (arguments overgiven explicit by skript parameters):
      progdir="${EXTPAR_HOME}/job"
    fi
    modul=$(echo "$*" | awk '{ print $1 }')
    args=$(echo "$*" | awk '{ sep=""; for(i=2; i<=NF; i++) { printf "%s%s", sep, $i; sep=" " } }')
    set +e
    modul=$(PATH=${progdir} whence ${modul})
    rc=$?
    set -e
    if (( rc != 0 )) ; then
      modul=modul_not_set
    fi
    if [[ ! -e ${modul} || ! -x ${modul} ]] ; then
      print -- "   ERROR: command '${modul}' not found or not executable."
      exit 201
    fi
    echo ">> Run $* ..."
    start=$(date +%s.%N)
    set +e
#   eval ./$* > ${logfile} 2>&1
    eval ${modul} ${args} > ${logfile} 2>&1
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
            ;;
        127)
            echo "   ERROR ${1%% *}: command not found"
            ;;
        130)
            echo "   ERROR ${1%% *}: script terminated by Ctrl-C"
            ;;
        *)
            echo "   ERROR ${1%% *}: fatal error - return code $rc"
            ;;
    esac
    echo "   execution time: $runtime s"
    [[ -f ${logfile} ]] && ls -l ${logfile} && cat ${logfile}
    [[ -f ${logfile} ]] && cat ${logfile} >> ${logfileforzip}

    # @ckoziar, 20210510: special logfiles of extpar v5.5
    logfile55=${modul##*/} ; logfile55=${logfile55%.py}.log
    [[ -f ${logfile55} ]] && ls -l ${logfile55} && cat ${logfile55}
    [[ -f ${logfile55} ]] && cat ${logfile55} >> ${logfileforzip}

    # @ckoziar, 20190926:
    if (( rc != 0 )) ; then
        exit ${rc}
    fi
}



eval $(routine_config)

TMP=${TMPDIR}
cd $TMP

logfile=extpar_modules_$(date +%Y%m%d%H%M%S).log
logfileforzip=%logfile4zip%
globeFlag=%globeFlag%
meritFlag=
nml_dir=%nml_dir%

ICON_GRID=%igrid%
GRIDDIR=%igrid_dir%
WORKDIR=%iwork_dir%
PROGDIR=%progdir_lc%
EXTPAR_HOME=~routfor/routfox/extpar

emissFlag=%emissFlag%
NoArtFlag=%NoArtFlag%

typeset merit_merge_flag='true'
typeset add_globe_flag='true'
typeset asterFlag=

#________________________________________________________________________________
# Names of executables

# Python
binary_era=extpar_era_to_buffer.py
binary_alb=extpar_alb_to_buffer.py
binary_tclim=extpar_cru_to_buffer.py
binary_ndvi=extpar_ndvi_to_buffer.py
binary_emiss=extpar_emiss_to_buffer.py
# Fortran
binary_aot=extpar_aot_to_buffer.exe.new
binary_soil=extpar_soil_to_buffer.exe.new
binary_flake=extpar_flake_to_buffer.exe.new
binary_lu=extpar_landuse_to_buffer.exe.new
binary_topo=extpar_topo_to_buffer.exe.new
binary_consistency_check=extpar_consistency_check.exe.new
binary_hwsdART=extpar_hwsdART_to_buffer.exe.new
#________________________________________________________________________________

# # cp $PROGDIR/$binary_alb $WORKDIR/
# cp ${EXTPAR_HOME}/job/$binary_alb $WORKDIR/
# cp $PROGDIR/$binary_ndvi $WORKDIR/
# cp $PROGDIR/$binary_emiss $WORKDIR/
# cp $PROGDIR/$binary_tclim $WORKDIR/
# cp $PROGDIR/$binary_lu $WORKDIR/
# cp $PROGDIR/$binary_topo $WORKDIR/
# cp $PROGDIR/$binary_aot $WORKDIR/
# cp $PROGDIR/$binary_soil $WORKDIR/
# cp $PROGDIR/$binary_flake $WORKDIR/
# cp $PROGDIR/$binary_consistency_check $WORKDIR/
# for file in $PROGDIR/*.py ; do
#   if [[ -s ${file} ]] ; then
#     cp -p ${file} $WORKDIR/
#   fi
# done

# 20220512: corresponding to "PBS -l cpunum_job=32"
#export OMP_NUM_THREADS=8
export OMP_NUM_THREADS=32

export NETCDF_OUTPUT_FILETYPE=NETCDF4

cd $WORKDIR

# @ckoziar, 20210211: uncomment following 2 lines:
# ncks='/hpc/sw/nco/4.8.1/bin/ncks'
# ncwa='/hpc/sw/nco/4.8.1/bin/ncwa'

# Assume GLOBE orography as default
cp INPUT_ORO_GLOBE INPUT_ORO
# @fprill/ckoziar, 20190924:
# Use ASTER orography for resolutions higher than 3.0 km
if [[ -s ${GRIDDIR}/${ICON_GRID} ]] ; then
  meshsize_km=$(ncks -M ${GRIDDIR}/${ICON_GRID} | \
                awk '/grid_root/  { GRIDROOT=$(NF-1)  }
                     /grid_level/ { GRIDLEVEL=$(NF-1) }
                     END { print(5050./(GRIDROOT*2**GRIDLEVEL))}
                    ' -
               )

  # @fprill/ckoziar, 20191204: cdo infov doesn't work in every case --> replace by ncwa/ncks
# min_lat_dom=$(cdo infov ${GRIDDIR}/${ICON_GRID} | \
#                grep -i latitude_vertices | awk '{ print $9 * 180. / atan2(0, -1) }')
  ncwa -y min -v latitude_vertices -O ${GRIDDIR}/${ICON_GRID} foo.$$.nc && \
  min_lat_dom=$(ncks -H -s "%f" foo.$$.nc | awk '/^[-+0-9.,]+$/ { printf "%f\n", $1 * 180. / atan2(0, -1) }')

# max_lat_dom=$(cdo infov ${GRIDDIR}/${ICON_GRID} | \
#                grep -i latitude_vertices | awk '{ print $11 * 180. / atan2(0, -1) }')
  ncwa -y max -v latitude_vertices -O ${GRIDDIR}/${ICON_GRID} foo.$$.nc && \
  max_lat_dom=$(ncks -H -s "%f" foo.$$.nc | awk '/^[-+0-9.,]+$/ { printf "%f\n", $1 * 180. / atan2(0, -1) }')

# min_lon_dom=$(cdo infov ${GRIDDIR}/${ICON_GRID} | \
#                grep -i longitude_vertices | awk '{ print $9 * 180. / atan2(0, -1) }')
  ncwa -y min -v longitude_vertices -O ${GRIDDIR}/${ICON_GRID} foo.$$.nc && \
  min_lon_dom=$(ncks -H -s "%f" foo.$$.nc | awk '/^[-+0-9.,]+$/ { printf "%f\n", $1 * 180. / atan2(0, -1) }')
# min_lon_dom=$(ncks -H -s "%f" foo.$$.nc | awk '/^[-+0-9.,]+$/ { lon = $1 * 180. / atan2(0, -1) ; while ( lon <= -180.0 ) { lon += 360.0 } printf "%f\n", lon }')

# max_lon_dom=$(cdo infov ${GRIDDIR}/${ICON_GRID} | \
#                grep -i longitude_vertices | awk '{ print $11 * 180. / atan2(0, -1) }')
  ncwa -y max -v longitude_vertices -O ${GRIDDIR}/${ICON_GRID} foo.$$.nc && \
  max_lon_dom=$(ncks -H -s "%f" foo.$$.nc | awk '/^[-+0-9.,]+$/ { printf "%f\n", $1 * 180. / atan2(0, -1) }')
# max_lon_dom=$(ncks -H -s "%f" foo.$$.nc | awk '/^[-+0-9.,]+$/ { lon = $1 * 180. / atan2(0, -1) ; while ( lon <= -180.0 ) { lon += 360.0 } printf "%f\n", lon }')

  if [[ "${meshsize_km}"  != +([-+0-9.,]) ]] ; then
    print -- "ERROR: Cannot determine horizontal resolution of target grid file!\n" >&2
    exit 202
  fi

  if [[ "${min_lat_dom}" != +([-+0-9.,]) || \
        "${max_lat_dom}" != +([-+0-9.,]) || \
        "${min_lon_dom}" != +([-+0-9.,]) || \
        "${max_lon_dom}" != +([-+0-9.,])    \
     ]] ; then
    print -- "ERROR: Cannot determine minimum/maximum latitude/longitude of target grid file!\n" >&2
    exit 203
  fi

  rm -f foo.$$.nc

  print -- "        grid file information for '${GRIDDIR}/${ICON_GRID}': "
  print -- "        --->  mesh size: ${meshsize_km} km"
  print -- "        --->  north latitude of vertices: ${max_lat_dom} deg"
  print -- "        --->  south latitude of vertices: ${min_lat_dom} deg"
  print -- "        --->  west longitude of vertices: ${min_lon_dom} deg"
  print -- "        --->  east longitude of vertices: ${max_lon_dom} deg"

  merit_merge_flag='true'
  if (( max_lat_dom > -60.0 && min_lat_dom > -60.0 ||
        max_lat_dom < -62.0 && min_lat_dom < -62.0 )) ; then
    merit_merge_flag=''
  fi

# @ckoziar, 20220711: bugfix (globeFlag AND meritFlag have to be empty!)
# if [[ -z "${globeFlag}" || -z "${meritFlag}" ]] ; then
#  if [[ -z "${globeFlag}" && -z "${meritFlag}" ]] ; then
    if (( meshsize_km < 3.0 &&  max_lat_dom < 60.0 && min_lat_dom > -60.0 )) ; then
	print -- "switch to ASTER orography because target horizontal resolution is higher than 3.0 km and "
        print -- "target domain is is within 60degN and 60degS"
      cp INPUT_ORO_ASTER INPUT_ORO
      asterFlag='true'
      add_globe_flag=

      aster_dir=$(awk '$1 == "raw_data_orography_path" { print substr($NF,2,length($NF)-2) }' INPUT_ORO)
      set -A arr_tilenames
      set -A arr_tilerows
      set -A arr_tilecols
      typeset -i row
      typeset -i col
      typeset -i xint
      typeset -i datelineFlag=0

      if (( min_lon_dom < -179.95 &&
            max_lon_dom >  179.95 ))
      then
        datelineFlag=1
        ncdump -v longitude_vertices ${GRIDDIR}/${ICON_GRID} > lulu.$$.txt
        linestot=$(wc -l lulu.$$.txt | awk '{ print $1 }')
        (( linestot += 1 ))
        awk -v linestot=${linestot} \
            '/^ *longitude_vertices *= */ \
             {
               linestot = NR
             }
             NR >= linestot \
             {
               for(i=1; i<=NF; i++)
               {
                 if($i ~ /[0-9\.,]+/)
                 {
                   lon = $i
                   gsub(/,/, "", lon)
                   lon = lon * 180. / atan2(0, -1)
                   while ( lon < 0. )
                   {
                     lon += 360.
                   }
                   printf "%.14f\n", lon
                 }
               }
             }' lulu.$$.txt > longitudes.$$.txt

        eval $(awk 'BEGIN \
                    {
                      max = -180.0
                      min = 180.0
                    }
                    {
                      if($1 > max)
                      {
                        max = $1
                      }
                      if($1 < min)
                      {
                        min = $1
                      }
                    }
                    END \
                    {
                      while(min > 180.)
                      {
                        min -= 360.
                      }
                      while(max > 180.)
                      {
                        max -= 360.
                      }
                      printf "min_lon_dom2=%.14f\n", min
                      printf "max_lon_dom2=%.14f\n", max
                    }' longitudes.$$.txt
              )

        print -- "        ---> REAL west longitude of vertices: ${min_lon_dom2} deg"
        print -- "        ---> REAL east longitude of vertices: ${max_lon_dom2} deg"
        min_lon_dom=${min_lon_dom2}
        max_lon_dom=${max_lon_dom2}
      fi

      for tile in $(eval_numlist -f %03d 1-240) ; do
        tile_name="ASTER_orig_T${tile}.nc"
        if [[ -s ${aster_dir}/${tile_name}.xz ]] ; then

          # prepare tile info file before by:
          # for i in ASTER_orig_T0??.nc ; do
          #   min_lon_tile=$(cdo sinfov $i 2>/dev/null | awk '$NF == "degrees_east" { print $3 }')
          #   max_lon_tile=$(cdo sinfov $i 2>/dev/null | awk '$NF == "degrees_east" { print $5 }')
          #   min_lat_tile=$(cdo sinfov $i 2>/dev/null | awk '$NF == "degrees_north" { print $5 }')
          #   max_lat_tile=$(cdo sinfov $i 2>/dev/null | awk '$NF == "degrees_north" { print $3 }')
          #   echo "${i} ${min_lon_tile} ${max_lon_tile} ${min_lat_tile} ${max_lat_tile}" >> \
          #        ${aster_dir}/ASTER_tile_domains.txt
          # done

          # Scanning mode for lon is west to east
          min_lon_tile=$(awk '$1 == "'"${tile_name}"'" { print $2 }' ${aster_dir}/ASTER_tile_domains.txt)
          max_lon_tile=$(awk '$1 == "'"${tile_name}"'" { print $3 }' ${aster_dir}/ASTER_tile_domains.txt)
          # Scanning mode for lat is north to south
          min_lat_tile=$(awk '$1 == "'"${tile_name}"'" { print $4 }' ${aster_dir}/ASTER_tile_domains.txt)
          max_lat_tile=$(awk '$1 == "'"${tile_name}"'" { print $5 }' ${aster_dir}/ASTER_tile_domains.txt)

          # echo "${tile}: longitude: ( ${min_lon_tile} ... ${max_lon_tile} )"
          # echo "${tile}: latitude: ( ${min_lat_tile} ... ${max_lat_tile} )"

          if (( ( datelineFlag > 0 &&
                  ( ( max_lat_dom > min_lat_tile && max_lat_dom < max_lat_tile ) || \
                    ( max_lat_dom > max_lat_tile && min_lat_dom < min_lat_tile ) || \
                    ( min_lat_dom > min_lat_tile && min_lat_dom < max_lat_tile ) ) && \
                  ( ( min_lon_dom > min_lon_tile && min_lon_dom < max_lon_tile ) || \
                    ( min_lon_dom < min_lon_tile && max_lon_tile <  180.0      ) || \
                    ( max_lon_dom > max_lon_tile && min_lon_tile > -180.0      ) || \
                    ( max_lon_dom > min_lon_tile && max_lon_dom < max_lon_tile ) ) ) || \
                ( ( ( max_lat_dom > min_lat_tile && max_lat_dom < max_lat_tile ) || \
                    ( max_lat_dom > max_lat_tile && min_lat_dom < min_lat_tile ) || \
                    ( min_lat_dom > min_lat_tile && min_lat_dom < max_lat_tile ) ) && \
                  ( ( min_lon_dom > min_lon_tile && min_lon_dom < max_lon_tile ) || \
                    ( min_lon_dom < min_lon_tile && max_lon_dom > max_lon_tile ) || \
                    ( max_lon_dom > min_lon_tile && max_lon_dom < max_lon_tile ) ) ) )) ; then
              arr_tilenames[${#arr_tilenames[*]}]=${tile_name}

              # ASTER consist of 240 tiles, 12 per row, 20 per column;
              # one tile is 30 (latitudinal) x 6 (longitudinal):
              (( xint = ( tile - 1 ) / 12 ))
              (( col = tile - ( 12 * xint ) ))
              (( row = ( tile + 11 ) / 12 ))
              arr_tilerows[${row}]=${row}
              arr_tilecols[${col}]=${col}
          fi
        fi
      done

      echo "needed ASTER tile names ${arr_tilenames[*]}"
      echo "needed ASTER rows ${arr_tilerows[*]}"
      echo "needed number of ASTER rows ${#arr_tilerows[*]}"
      echo "needed ASTER columns ${arr_tilecols[*]}"
      echo "needed number of ASTER columns ${#arr_tilecols[*]}"

      # prepare strings for pshell/namelist:
      typeset cmdfile_pshell="cmd_xz.pshell.$$"
      rm -f ${cmdfile_pshell}
      typeset parallel_pshell=8
      str=''
      blank=''
      for i in ${arr_tilenames[*]} ; do
        # copy ASTER tiles to $WORKDIR und uncompress them (do later on by pshell):
        print -- "[[ -s ${WORKDIR}/${i} ]] || cp ${aster_dir}/${i}.xz ${WORKDIR}/ && xz -dv ${WORKDIR}/${i}.xz" >> ${cmdfile_pshell}

        str=$(printf -- "%s%s'%s'" "${str}" "${blank}" ${i})
        blank=' '
      done

      # modify namelist for ASTER:
      sed "s/ topo_files = .*$/ topo_files = ${str}/" INPUT_ORO > foo1.$$
      sed "s/ ntiles_row = .*$/ ntiles_row = ${#arr_tilerows[*]}/" foo1.$$ > foo2.$$
      sed "s/ ntiles_column = .*$/ ntiles_column = ${#arr_tilecols[*]}/" foo2.$$ > foo3.$$
      sed "s: raw_data_orography_path = .*$: raw_data_orography_path = '${WORKDIR}/':" foo3.$$ > INPUT_ORO

      rm -f foo?.$$
      cat INPUT_ORO

      # uncompress ASTER tiles in parallel:
      if [[ -s ${cmdfile_pshell} ]] ; then
        cat ${cmdfile_pshell}

        date
        set +e
        pshell -p ${parallel_pshell} -taback -taback_fast_sleep -L -T -f ${cmdfile_pshell}
        rc=$?
        set -e
        date

        if (( rc != 0 )) ; then
          print -- "ERROR: Cannot uncompress ASTER tiles by pshell!\n" >&2
          exit ${rc}
        else
          set -A arr_rc_parallel -- $(<rc_parallel)
          rc=$(print -- ${arr_rc_parallel[*]} | tr ' ' '+' | bc)
          if (( rc != 0 )) ; then
            print -- "ERROR: Cannot uncompress all ASTER tiles by pshell!\n" >&2
            ls -l ${WORKDIR}/ASTER_orig_T*
            exit ${rc}
          else
            ls -l ${WORKDIR}/ASTER_orig_T*
          fi
        fi
        rm -f ${cmdfile_pshell}
      fi
    fi
#  else
#    if (( meshsize_km < 1.5 )) ; then
#      # @ckoziar, 20191108:
#      # don't generate SSO parameters for GLOBE
#      # in case of horizontal resolution < 1.5 km:
#      sed "s/ lsso_param = .*$/ lsso_param = .FALSE./" INPUT_ORO > foo.$$ && mv foo.$$ INPUT_ORO
#    fi
#  fi
fi

# final namelists modifications:
if [[ -n "${nml_dir}" && -d ${nml_dir} ]] ; then
  for file in ${nml_dir}/INPUT_* ; do
    if [[ -s ${file} ]] ; then
      if [[ "${file}" == *.mod ]] ; then
        ifile=${file##*/} ; ifile=${ifile%.mod}
        if [[ -s ${ifile} ]] ; then
          print -- "modify namelist '${ifile}' by following instructions:"
          cat ${file}
          echo
          print -- "namelist before modification:"
          cat ${ifile}
          set +e
          # @ckoziar, 20191108:
          # convert AWK to Perl because named arrays are not supported by gawk v3:
          # ~fprill/software/gawk-4.1.1-build/bin/awk -f ${EXTPAR_HOME}/job/modify_nml.awk -v modfile=${file} ${ifile} > ${ifile}.$$.new
          perl > ${ifile}.$$.new <<__HERE__
\$modfile = "${file}";
\$nmlfile = "${ifile}";
\$active_nml = "";
if(open(INPUT, "<\$modfile"))
{
  printf("!  read modification table from '%s'.\n!\n", \$modfile);
  while(<INPUT>)
  {
    @a = split('\|', \$_);
    for(\$i = 0; \$i < scalar(@a); \$i++)
    {
      \$a[\$i] =~ s/^\s+//; \$a[\$i] =~ s/\s+$//; # trim
    }
    if(scalar(@a) == 3)
    {
      \$mod_dict{"&" . \$a[0]}{\$a[1]} = \$a[2];
    }
  }
  close(INPUT);
  foreach \$nml ( sort keys(%mod_dict) )
  {
    foreach \$key ( sort keys %{ \$mod_dict{\$nml} } )
    {
      printf("!  %s %s %s\n", \$nml, \$key, \$mod_dict{\$nml}{\$key});
    }
  }
  printf("!\n!  original namelist file: '%s'.\n\n", \$nmlfile);
}
else
{
  print STDERR "Error: Cannot open namelist modification file '\$modfile' for reading!\n";
  exit 204;
}

if(open(INPUT, "<\$nmlfile"))
{
  while(\$line = <INPUT>)
  {
    chomp(\$line);
    if(\$line =~ /^&[a-zA-Z_]+/)
    { # parse namelist start:
      \$active_nml = \$line;
      \$active_nml =~ s/^\s+//; \$active_nml =~ s/\s+$//; # trim
      printf("%s\n", \$line);
      %visited_arr = ();
      next;
    }
    elsif(\$line =~ /\/\s*$/)
    { # parse namelist close:
      if(defined(\$mod_dict{\$active_nml}))
      {
        foreach \$key ( sort keys %{ \$mod_dict{\$active_nml} } )
        {
          if(!defined(\$visited_arr{\$key}))
          {
            printf("%s = %s\n", \$key, \$mod_dict{\$active_nml}{\$key});
          }
        }
      }
      \$active_nml = "";
      printf("%s\n", \$line);
      next;
    }
    else
    { # search for namelist keys to replace:
      @a = split('=', \$line);
      \$a[0] =~ s/^\s+//; \$a[0] =~ s/\s+$//; # trim
      if(defined(\$mod_dict{\$active_nml}) && defined(\$mod_dict{\$active_nml}{\$a[0]}))
      {
        printf("%s = %s\n", \$a[0], \$mod_dict{\$active_nml}{\$a[0]});
        \$visited_arr{\$a[0]} = 1;
      }
      else
      {
        printf("%s\n", \$line);
      }
    }
  }
}
else
{
  print STDERR "Error: Cannot open namelist input file '\$nmlfile' for reading!\n";
  exit 205;
}

exit 0;
__HERE__
          rc=$?
          set -e
          if (( rc != 0 )) ; then
            print -- "ERROR: Cannot modify namelist '${ifile}' by '${file}'!\n" >&2
            exit ${rc}
          fi

          mv ${ifile}.$$.new ${ifile}
          print -- "namelist after modification:"
          cat ${ifile}
        fi
      else
        ifile=${file##*/}
        print -- "total namelist replacing by 'cp ${file} ${PWD}/${ifile}'"
        set +e
        ls -l ${file} ${PWD}/${ifile}
        print -- "diff ${file} ${PWD}/${ifile}:"
        diff ${file} ${PWD}/${ifile}
        set -e
        cp ${file} ${PWD}/${ifile}
      fi
    fi
  done
fi

cd ${WORKDIR}

run_python ${binary_era}
cdo copy -selvar,T_SEA,W_SNOW era_buffer.nc era_sst_buffer.nc
cdo copy -selvar,TOPO_CLIM,T_2M_CLIM era_buffer.nc era_2t_buffer.nc

run_python ${binary_alb}
run_python ${binary_ndvi}
run_python ${binary_emiss}

# Fortran
run_command ${binary_aot}
run_command ${binary_soil}
run_command ${binary_flake}

run_command ${binary_lu}
if [[ -z "${meritFlag}" || ( -n "${meritFlag}" && -n "${merit_merge_flag}" ) || ( -n "${meritFlag}" && -n "${add_globe_flag}" ) ]] ; then
  run_command ${binary_topo}
fi
if [[ -n "${meritFlag}" ||  -n "${add_globe_flag}" || ( -z  "${asterFlag}" ) ]] ; then
  # input and output fields
  typeset file_globe="${WORKDIR}/topography_BUFFER_GLOBE.nc" # :GLOBE-DEM-File
  typeset file_merit="${WORKDIR}/topography_BUFFER_MERIT.nc" # :MERIT/REMA-File
  typeset file_final="${WORKDIR}/topography_BUFFER.nc"       # :merged output file
  typeset ftmp_globe="${TMPDIR}/topography_BUFFER_GLOBE.nc" # :temporary GLOBE-DEM-File
  typeset ftmp_merit="${TMPDIR}/topography_BUFFER_MERIT.nc" # :temporary MERIT/REMA-File
  typeset f_temp="newtopo.nc"

  typeset south='-62'
  typeset north='-60'
  typeset lat="lat"
  # @ckoziar, 20220504: 
  # don't process SSO parameters for GLOBE in case of horizontal resolution < 1.5 km:
# typeset ee="HSURF,SSO_STDH,SSO_SIGMA,SSO_THETA,SSO_GAMMA,Z0_TOPO,FR_LAND_TOPO"
  typeset ee
  if (( meshsize_km < 1.5 )) ; then
    ee="HSURF,Z0_TOPO,FR_LAND_TOPO"
  else
    ee="HSURF,SSO_STDH,SSO_SIGMA,SSO_THETA,SSO_GAMMA,Z0_TOPO,FR_LAND_TOPO"
  fi

  if [[ -n "${merit_merge_flag}" ||  -n "${add_globe_flag}" ]] ; then
    # rename GLOBE file because MERIT will result in the same filename:
    if [[ -s ${file_final} ]] ; then
      echo "mv ${file_final} ${file_globe}"
      mv ${file_final} ${file_globe}
    else
      print -- "ERROR: GLOBE orography file '${file_final}' is empty or non-existent!\n" >&2
      exit 97
    fi
  fi

  print -- "processing MERIT orography..."
  cp INPUT_ORO_MERIT INPUT_ORO

  # prepare strings for pshell/namelist:
  typeset cmdfile_pshellm="cmd_merit.pshell.$$"
  rm -f ${cmdfile_pshellm}
  typeset parallel_pshellm=8

  merit_dir=$(awk '$1 == "raw_data_orography_path" { print substr($NF,2,length($NF)-2) }' INPUT_ORO)
  for i in ${merit_dir}/REMA_BKG_S60-S90_*.nc4 ; do
    # copy REMA files to $WORKDIR (do later on by pshell):
    print -- "cp ${i} ${WORKDIR}/" >> ${cmdfile_pshellm}
  done
  for i in ${merit_dir}/MERIT_*.nc4.xz ; do
    # copy MERIT files to $WORKDIR und decompress (do later on by pshell):
    file=$(basename ${i})
    print -- "cp ${i} ${WORKDIR}/ && xz -dv ${WORKDIR}/${file}" >> ${cmdfile_pshellm}
  done

  # copy MERIT tiles in parallel:
  if [[ -s ${cmdfile_pshellm} ]] ; then
    cat ${cmdfile_pshellm}

    date
    set +e
    pshell -p ${parallel_pshellm} -taback -taback_fast_sleep -L -T -f ${cmdfile_pshellm}
    rc=$?
    set -e
    date

    if (( rc != 0 )) ; then
      print -- "ERROR: Cannot copy MERIT tiles by pshell!\n" >&2
      exit ${rc}
    else
      set -A arr_rc_parallel -- $(<rc_parallel)
      rc=$(print -- ${arr_rc_parallel[*]} | tr ' ' '+' | bc)
      if (( rc != 0 )) ; then
        print -- "ERROR: Cannot copy all MERIT tiles by pshell!\n" >&2
        ls -l ${WORKDIR}/REMA_BKG_S60-S90_*
        exit ${rc}
      else
        ls -l ${WORKDIR}/REMA_BKG_S60-S90_*
      fi
    fi
    rm -f ${cmdfile_pshellm}
  fi

  sed "s: raw_data_orography_path = .*$: raw_data_orography_path = '${WORKDIR}/':" INPUT_ORO > INPUT_ORO.lulu.$$ && \
  mv INPUT_ORO.lulu.$$ INPUT_ORO

  run_command ${binary_topo}

  if [[ -n "${merit_merge_flag}" ]] ; then
    # merge topography_BUFFER_GLOBE.nc and topography_BUFFER_MERIT.nc
    # to topography_BUFFER.nc
    print -- "merge GLOBE and MERIT orography in the region of -60...-62 degree south..."

    if [[ -s ${file_final} ]] ; then
      echo "mv ${file_final} ${file_merit}"
      mv ${file_final} ${file_merit}
    else
      print -- "ERROR: file '${file_final}' is empty or non-existent!\n" >&2
      exit 97
    fi

    # Extract latitude and ${ee} from the original files
    cdo selvar,${lat} ${file_merit} ${ftmp_merit}.${lat}
    cdo selvar,${ee}  ${file_merit} ${ftmp_merit}.topo
    cdo selvar,${ee}  ${file_globe} ${ftmp_globe}.topo

    # Create the masks for the files
    cdo add -lec,${south} ${ftmp_merit}.${lat} \
            -gec,${north} ${ftmp_merit}.${lat} ${ftmp_merit}.mask
    cdo nec,1 ${ftmp_merit}.mask ${ftmp_globe}.mask

    # Create masked vaules of the files
    cdo mul ${ftmp_merit}.topo ${ftmp_merit}.mask ${ftmp_merit}.masked
    cdo mul ${ftmp_globe}.topo ${ftmp_globe}.mask ${ftmp_globe}.masked

    # Add both masked files to new file
    cdo add ${ftmp_merit}.masked ${ftmp_globe}.masked ${f_temp}

    # remove orography fields from input file
    cdo delname,${ee} ${file_merit} tempfile.nc

    ncks -A ${f_temp} tempfile.nc

    mv tempfile.nc ${file_final}
  fi

  if [[ -s ${file_final} ]] ; then
    ls -l ${file_final}
    if [[ -n "${merit_merge_flag}" ]] ; then
      ls -l ${file_globe} ${file_merit}
    fi
  else
    print -- "ERROR: file '${file_final}' is empty or non-existent!\n" >&2
    ls -l ${file_globe} ${file_merit} ${file_final}
    exit 97
  fi

# exchange FR_LAND_TOPO in MERIT-BUFFER file with GLOBE for correct FR_LAKE determination

  mv topography_BUFFER.nc topography_BUFFER_MERIT.nc
  cdo selvar,FR_LAND_TOPO topography_BUFFER_GLOBE.nc FR_LAND_TOPO_GLOBE.nc
  cdo selvar,FR_LAND_TOPO topography_BUFFER_MERIT.nc FR_LAND_TOPO_MERIT.nc

  cdo delname,FR_LAND_TOPO topography_BUFFER_MERIT.nc topography_BUFFER_tmp1.nc

  cdo merge topography_BUFFER_tmp1.nc FR_LAND_TOPO_GLOBE.nc topography_BUFFER.nc

  # cleanup:
  rm -f ${ftmp_merit}.mask ${ftmp_merit}.masked
  rm -f ${ftmp_globe}.mask ${ftmp_globe}.masked

  for i in ${WORKDIR}/REMA_BKG_S60-S90_*.nc4 ; do
    echo "rm -f ${i}"
    rm -f ${i}
  done

  for i in ${WORKDIR}/MERIT_*.nc4 ; do
    echo "rm -f ${i}"
    rm -f ${i}
  done
fi
# T_CLIM needs the orography output file - call after TOPO
run_python ${binary_tclim}

# @jhelmert/ckoziar, 20191212:
# avoid 'SIGSEGV, segmentation fault' in consistency_check
# by switching to intel module version 2018.2.199:
# module_count=$(module list 2>&1 | grep -c 'intel/')
# if [[ -z "${module_count}" ]] || (( module_count == 0 )) ; then
#   module load intel/2018.2.199
# else
#   module swap intel/2018.2.199
# fi

echo
echo "checking existence of INPUT_EMISS before consistency check:"
ls -l INPUT_*
echo

run_command ${binary_consistency_check}

# runhwsdART only if NoArtFlag is empty, i.e. not activated by -x on command line
if [[ -n "${NoArtFlag}" ]] ; then
      if (( meshsize_km > 1.5 )) ; then
      print -- "Horizontal resolution of target grid is coarser than 1.5 km: Computation of HWSD-ART soil data"
      run_command ${binary_hwsdART}
      else
      print -- "Horizontal resolution of target grid is higher than 1.5 km: Computation of HWSD-ART soil data not feasible"
      fi
fi
