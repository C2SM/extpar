#!/bin/ksh -l

#PBS -w NAME=extpar_2grib2
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
#PBS -l elapstim_req=1800
#PBS -o %loggingfile%
#PBS -j o

set -e
set -x

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

module list

# @ckoziar, 20210507:
ncatted='/hpc/sw/nco/4.8.1/bin/ncatted'

eval $(routine_config)

TMP=${TMPDIR}
cd $TMP

for grid_id in %igrid%
do
work_dir=%iwork_dir%
cd $work_dir

# @fprill/ckoziar, 20190903: revise grid_type - as long as extpar itself don't set it:
# (cdo revise it for grid sizes > 999 points only)
#ncatted -O -a CDI_grid_type,,c,c,"unstructured" external_parameter_icon_${grid_id%.nc}_tiles.nc
${ncatted} -O -a CDI_grid_type,,c,c,"unstructured" \
              -a standard_name,time,c,c,time \
              -a units,time,c,c,"days as %Y%m%d.%f"  external_parameter_icon_${grid_id%.nc}_tiles.nc

# extract the fields needed: for GRIB LON and LAT are needed and mapped to CLON and CLAT
cdo copy -selvar,SOILTYP,FR_LAND,PLCOV_MX,LAI_MX,RSMIN,FOR_D,FOR_E,EMIS_RAD,ROOTDP,Z0,NDVI_MAX,topography_c,SSO_STDH,SSO_THETA,SSO_GAMMA,SSO_SIGMA,T_CL,FR_LAKE,DEPTH_LK,LU_CLASS_FRACTION,NDVI,NDVI_MRAT,AER_BC,AER_DUST,AER_ORG,AER_SO4,AER_SS,ALB,ALNID,ALUVD,T_SEA,W_SNOW,lon,lat external_parameter_icon_${grid_id%.nc}_tiles.nc test_selvar.nc
# change names -- write lon -> CLON, lat -> CLAT to have unit degN, degE for coordinated in GRIB
cdo copy -chname,topography_c,HSURF -chname,T_CL,T_2M_CL -chname,LU_CLASS_FRACTION,FR_LUC -chname,AER_BC,AER_BC12 -chname,AER_DUST,AER_DUST12 -chname,AER_ORG,AER_ORG12 -chname,AER_SO4,AER_SO412 -chname,AER_SS,AER_SS12 -chname,ALB,ALB_DIF12 -chname,ALNID,ALB_NI12 -chname,ALUVD,ALB_UV12 -chname,lon,CLON -chname,lat,CLAT test_selvar.nc new-names-select-external_parameter_icon_${grid_id%.nc}_tiles.nc

# transform to GRIB2, field by field and use grib_filter for corrections
for var in SOILTYP FR_LAND PLCOV_MX LAI_MX RSMIN FOR_D FOR_E EMIS_RAD ROOTDP Z0 NDVI_MAX HSURF SSO_STDH SSO_THETA SSO_GAMMA SSO_SIGMA T_2M_CL FR_LAKE DEPTH_LK FR_LUC NDVI NDVI_MRAT AER_BC12 AER_DUST12 AER_ORG12 AER_SO412 AER_SS12 ALB_DIF12 ALB_NI12 ALB_UV12 T_SEA W_SNOW CLON CLAT ; do

var3=$( echo ${var} | head -c3)

if [[ "${var}" == "CLON" || "${var}" == "CLAT" ]] ; then
  cdo -f grb2 -b 24 copy -selvar,${var} new-names-select-external_parameter_icon_${grid_id%.nc}_tiles.nc external_parameter_icon_${grid_id%.nc}_tiles_${var}.g2
else
  cdo -f grb2 -b 16 copy -selvar,${var} new-names-select-external_parameter_icon_${grid_id%.nc}_tiles.nc external_parameter_icon_${grid_id%.nc}_tiles_${var}.g2
fi

#Generate variables for time stamps in local sections

YEAR=$(date  +%Y)
MONTH=$(date  +%m)
DAY=$(date  +%d)
HOUR=$(date  +%H)
MINUTE=$(date  +%M)
SECOND=$(date  +%S)


cat > script_extpar_filter_rules << end_input
# Aufruf: grib_filter  -o <grib2-file> -f -v filter_rules_fr_luc  <grib2-file-first>
#

transient sname = "${var}" ;
transient sname3 = "${var3}" ;
print "*************** [sname] **** [sname3] ****************************** " ;

print "******************************************************************** " ;
print "************ Variable   [shortName] [stepType] [typeOfLevel] [level]";

#------------------------------------------------------------------------------------
# Gitterdefinitionen
   set shapeOfTheEarth=6;
   set numberOfGridInReference=1;

# lokale Sektion besetzen
#if(centre==78 && (localDefinitionNumber==254 || localDefinitionNumber==253)) {
 set grib2LocalSectionPresent =1;
 set localDefinitionNumber = 254;
 set localHostIdentifier = 255;

 set localCreationDateYear = ${YEAR};
 set localCreationDateMonth = ${MONTH};
 set localCreationDateDay = ${DAY};
 set localCreationDateHour = ${HOUR};
 set localCreationDateMinute = ${MINUTE};
 set localCreationDateSecond = ${SECOND};
 set localValidityDateYear = 0;
 set localValidityDateMonth = 0;
 set localValidityDateDay = 0;
 set localValidityDateHour = 0;
 set localValidityDateMinute = 0;
 set localValidityDateSecond = 0;
 set localInformationNumber = 255;
 set localVersionNumber = 0;
 set localNumberOfExperiment = 255;

# set setLocalDefinition =1;
# set backgroundGeneratingProcessIdentifier = bgp         ;
# set localNumberOfExperiment = nexp ;
#
# print "+++ Lokale Sektion besetzen: bgp=[backgroundGeneratingProcessIdentifier],nexp=[localNumberOfExperiment]";
#  }

# invariante Daten (d=invar)
if(dataDate == 10101) {
  print "+++ Invariant data d=invar" ;
  set productionStatusOfProcessedData = 0;
  set typeOfGeneratingProcess = 196 ;
}
# Klimadaten (d=1111mm1111)
#if(yearOfCentury==11 && day==11 && hour==11 && centuryOfReferenceTimeOfData==12) {
if(year==1111 && day==11) {
  print "+++ Climatic data d=1111mm1111" ;
  set hour = 11;
  set productionStatusOfProcessedData = 0;
  set typeOfGeneratingProcess = 9 ;
}

  set typeOfProcessedData=0;
  set significanceOfReferenceTime=0;
  set generatingProcessIdentifier=1;


#*********************************************************************************
if (sname3 is 'AER') {
  set productDefinitionTemplateNumber=42;
  set typeOfTimeIncrement=3;
  set shortName = sname;
}
#*********************************************************************************
if (sname is 'FR_LUC') {
  set productDefinitionTemplateNumber=53;
  set shortName ='FR_LUC' ;
  set partitionTable=243;
  set numberOfPartitions=23;
# set partitionItems={24,25,26,27,28,2,3,29,30,31,32,33,34,14,35,36,37,38,22,19,22,21,39};
#                                                                                    ##                            ##
  set partitionItems={24,25,26,27,28,2,3,29,30,31,32,33,34,13,35,36,37,38,22,19,20,21,39};

 switch (level) {

  case "1":
    set partitionNumber=24;

  case "2":
    set partitionNumber=25;

  case "3":
    set partitionNumber=26;

  case "4":
    set partitionNumber=27;

  case "5":
    set partitionNumber=28;

  case "6":
    set partitionNumber=2;

  case "7":
    set partitionNumber=3;

  case "8":
    set partitionNumber=29;

  case "9":
    set partitionNumber=30;

  case "10":
    set partitionNumber=31;

  case "11":
    set partitionNumber=32;

  case "12":
    set partitionNumber=33;

  case "13":
    set partitionNumber=34;

  case "14":
    set partitionNumber=13;

  case "15":
    set partitionNumber=35;

  case "16":
    set partitionNumber=36;

  case "17":
    set partitionNumber=37;

  case "18":
    set partitionNumber=38;

  case "19":
    set partitionNumber=22;

  case "20":
    set partitionNumber=19;

  case "21":
    set partitionNumber=20;

  case "22":
    set partitionNumber=21;

  case "23":
    set partitionNumber=39;

  default:
    print " Unexpected level - land use class !";
}

print "+++ Set typeOfLevel to surface and level to missing" ;
set typeOfLevel = "surface";
set scaledValueOfFirstFixedSurface = missing();
set scaleFactorOfFirstFixedSurface = missing();
set scaledValueOfSecondFixedSurface = missing();
set scaleFactorOfSecondFixedSurface = missing();
}
#*********************************************************************************

switch (sname) {

 case "SOILTYP":
  set typeOfOriginalFieldValues=1;
  set shortName=sname;

 case "PLCOV_MX" :
  set productDefinitionTemplateNumber=8;
  set typeOfTimeIncrement=2;
  set shortName=sname;

 case "LAI_MX" :
  set productDefinitionTemplateNumber=8;
  set typeOfTimeIncrement=2;
  set shortName=sname;

 case "RSMIN" :
  set productDefinitionTemplateNumber=8;
  set typeOfStatisticalProcessing=3;
  set typeOfTimeIncrement=2;
  set shortName=sname;

 case "NDVI_MAX" :
  set productDefinitionTemplateNumber=8;
  set typeOfTimeIncrement=2;
  set shortName=sname;

 case "T_2M_CL" :
  set productDefinitionTemplateNumber=8;
  set typeOfStatisticalProcessing=0;
  set typeOfTimeIncrement=3;
  set shortName=sname;

 case "NDVI" :
  set productDefinitionTemplateNumber=8;
  set typeOfStatisticalProcessing=0;
  set typeOfTimeIncrement=3;
  set shortName=sname;

 case "NDVI_MRAT" :
  set productDefinitionTemplateNumber=8;
  set typeOfStatisticalProcessing=0;
  set typeOfTimeIncrement=3;
  set shortName=sname;

 case "ALB_DIF12" :
  set productDefinitionTemplateNumber=8;
  set typeOfStatisticalProcessing=0;
  set typeOfTimeIncrement=3;
  set shortName=sname;

 case "ALB_NI12" :
  set productDefinitionTemplateNumber=8;
  set typeOfStatisticalProcessing=0;
  set typeOfTimeIncrement=3;
  set shortName=sname;

 case "ALB_UV12" :
  set productDefinitionTemplateNumber=8;
  set typeOfTimeIncrement=3;
  set shortName=sname;

 case "T_SEA" :
  set productDefinitionTemplateNumber=8;
  set typeOfStatisticalProcessing=0;
  set typeOfTimeIncrement=3;
  set shortName=sname;

 case "W_SNOW" :
  set productDefinitionTemplateNumber=8;
  set typeOfStatisticalProcessing=0;
  set typeOfTimeIncrement=3;
  set shortName=sname;

 default:
  set shortName=sname;
}
#*********************************************************************************

# set packingType="grid_simple";
#set dataRepresentationTemplateNumber=0;
#set bitsPerValue=16;

print "+++ Surface: Set levels to missing" ;
if (typeOfFirstFixedSurface == 1) {
set scaledValueOfFirstFixedSurface = missing();
set scaleFactorOfFirstFixedSurface = missing();
}
if (typeOfSecondFixedSurface == 255) {
set scaledValueOfSecondFixedSurface = missing();
set scaleFactorOfSecondFixedSurface = missing();
}

print "+++ WRITE";
write ;

end_input


grib_filter script_extpar_filter_rules -o filter_out_${var}.g2 external_parameter_icon_${grid_id%.nc}_tiles_${var}.g2


cat filter_out_*.g2 > external_parameter_icon_${grid_id%.nc}_tiles.g2

 done
done

# check for generatingProcessIdentifier:

gpi_file=external_parameter_icon_${grid_id%.nc}_tiles_filter.g2
cp external_parameter_icon_${grid_id%.nc}_tiles.g2 ${gpi_file}

# Test boundaries of domain for ICON-D2:
#/hpc/rwork0/routfor/routfox/icon/grids/public/edzw$  grib_list_dwd -v -P min,max -w shortName=CLON/CLAT  icon_extpar_0047_R19B07_L_20200527_tiles.g2
#File icon_extpar_0047_R19B07_L_20200527_tiles.g2
#count dataDate dataTime:s step stepUnits:s shortName level:d typeOfLevel min max
#97 10101 0000 0 h CLAT 0 surface 43.0441 58.1647
#98 10101 0000 0 h CLON 0 surface -4.16164 20.5444

# Test boundaries of domain for ICON-N02:
#/hpc/rwork0/routfor/routfox/icon/grids/public/edzw$  grib_list_dwd -v -P min,max -w shortName=CLON/CLAT  icon_extpar_0027_R03B08_N02_20200917_tiles.g2
#File icon_extpar_0027_R03B08_N02_20200917_tiles.g2
#count dataDate dataTime:s step stepUnits:s shortName level:d typeOfLevel min max
#97 10101 0000 0 h CLAT 0 surface 28.9257 71.0702
#98 10101 0000 0 h CLON 0 surface -24.6659 63.6407

set +e
grib_get -P min,max -w shortName=CLON ${gpi_file} | read west east
grib_get -P min,max -w shortName=CLAT ${gpi_file} | read south north
set -e

if [[ ( "${west}" != ?(+|-)+([0-9])?(.)*([0-9]) && \
        "${west}" != ?(+|-).+([0-9]) ) || \
      ( "${east}" != ?(+|-)+([0-9])?(.)*([0-9]) && \
        "${east}" != ?(+|-).+([0-9]) ) || \
      ( "${south}" != ?(+|-)+([0-9])?(.)*([0-9]) && \
        "${south}" != ?(+|-).+([0-9]) ) || \
      ( "${north}" != ?(+|-)+([0-9])?(.)*([0-9]) && \
        "${north}" != ?(+|-).+([0-9]) ) ]] ; then
  print -- "Unvalid domain boundaries detected: W=${west} E=${east} S=${south} N=${north}"
  print -- "Cannot set generatingProcessIdentifier!"
  exit 0
else
  print -- "W=${west} E=${east} S=${south} N=${north}"
fi

if (( south > 42 && north < 60 && west > -5 && east < 22 )) ; then
  print -- "DOMAIN ICON-L (D2): set generatingProcessIdentifier=11"
  gpid=11
elif (( south > 27 && north < 72 && west > -26 && east < 65 )) ; then
  print -- "DOMAIN N02: set generatingProcessIdentifier=2"
  gpid=2
else
  print -- "no N02 or L domain over Europe detected; assume global grid: set generatingProcessIdentifier=1"
  gpid=1
fi

grib_set -s generatingProcessIdentifier=${gpid} ${gpi_file} external_parameter_icon_${grid_id%.nc}_tiles.g2

exit 0
