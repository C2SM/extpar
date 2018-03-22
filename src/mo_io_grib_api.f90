!+ Fortran module containing input/output subroutines with the GRIB_API
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V1_2         2011/03/25 Hermann Asensio
!  Update documentation
! V1_3         2011/04/19 Hermann Asensio
! add support for GRIB1 and GRIB2
! V1_4         2011/04/21 Hermann Asensio, Anne Roche
!  set correct stepType(GRIB1 and GRIB2) and timeRangeIndicator (GRIB1) (ha)
!  clean up of code (roa)
! V1_8         2013-03-12 Frank Brenner
!  small bug fix regarding grib_set "stepType" to "avg" with new GRIB API 
!  versions (above 1.9.9)         
! V1_11        2013/04/16 Juergen Helmert
!  Specify resolutionAndComponentFlags       
! V1_14        2014-07-18 Juergen Helmert
!  Combined COSMO Release
! V2_2         2015-02-23 Juergen Helmert
!  Adaption for usage of standard GRIB-API template         
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module containing input/output subroutines with the GRIB_API
!! and variables for values of the GRIB_API keys (with the same name as the keys)
MODULE mo_io_grib_api
  
  USE grib_api

  !> kind parameters are defined in MODULE data_parameters
  USE mo_kind, ONLY: wp
  USE mo_kind, ONLY: i8
  USE mo_kind, ONLY: i4

  USE mo_io_units, ONLY: filename_max



  IMPLICIT NONE

  PUBLIC

  INTEGER, PARAMETER :: keylen_max = 80 !< maximum length for the length of keys of type character

  CHARACTER (LEN=keylen_max) :: gridType !< type of grid, e.g. rotated_ll

  ! Edition independent keys: geography, rotated_ll
  REAL        :: longitudeOfFirstGridPointInDegrees !< Longitude in degrees of the first grid point
  REAL        :: longitudeOfLastGridPointInDegrees  !< Longitude in degrees of the last grid point
  REAL        :: latitudeOfFirstGridPointInDegrees  !< Latitude in degrees of the first grid point
  REAL        :: latitudeOfLastGridPointInDegrees   !< Latitude in degrees of the last grid point
  INTEGER     :: Ni !< Number of points along a parallel
  INTEGER     :: Nj !< Number of points along a meridian
  REAL        :: iDirectionIncrementInDegrees !< i (longitude) direction increment in degrees
  REAL        :: jDirectionIncrementInDegrees !< j (latitude) direction increment in degrees
  INTEGER     :: iScansNegatively !< i (longitude) is scanned in the negative direction
  INTEGER     :: jScansPositively !< j (latitude) is scanned in the negative direction
  INTEGER     :: jPointsAreConsecutive !< The points in the j direction are scanned consecutively
  REAL        :: longitudeOfSouthernPoleInDegrees !< Longitude in degrees of the southern pole for rotated grid
  REAL        :: latitudeOfSouthernPoleInDegrees  !< Latitude in degrees of the southern pole for rotated grid
  REAL        :: angleOfRotationInDegrees !< Angle of rotation in degrees

  INTEGER     :: numberOfValues !< number of data values in GRIB, read only
  
  ! edition independent keys for parameters
  INTEGER                    :: paramID   !< Unique identifier of the parameter. 
                                          !! A parameter in different formats (grib1/grib2) is identified by the same ID.
  CHARACTER (LEN=keylen_max) :: shortName !< short name of the parameter
  INTEGER                    :: centre !< Identification of originating/generating centre,
  INTEGER, PARAMETER :: dwd_id_grib =78   !< the number 78 is for DWD, Offenbach
!roa mch>
  INTEGER, PARAMETER :: mch_id_grib=215
!roa mch<


  ! keys for represenation of data values
  INTEGER     :: bitsPerValue       !< Number of bits containing each packed value
  INTEGER     :: decimalScaleFactor !< Units decimal scale factor D (10**decimalScaleFactor)

  INTEGER     :: changeDecimalPrecision !< function key used to set the decimal precision
  INTEGER     :: decimalPrecision !< digits of decimal precision after decimal point

  ! local definitions for DWD, Offenbach
  INTEGER     :: localUsePresent       !< use local extension (read only)
  INTEGER     :: localDefinitionNumber !< local defintion number
  INTEGER     :: localElementNumber    !< local element number
  INTEGER (KIND=8)     :: localDecodeDate       !< date of decoding the GRIB record, local extension DWD, read only
  INTEGER     :: localDecodeDateYear   !< year of decoding data the GRIB record, local extension DWD
  INTEGER     :: localDecodeDateMonth  !< month of decoding data the GRIB record, local extension DWD
  INTEGER     :: localDecodeDateDay    !< day of decoding data the GRIB record, local extension DWD
  INTEGER     :: localDecodeDateHour   !< hou rof decoding data the GRIB record, local extension DWD
  INTEGER     :: localDecodeDateMinute !< minute of decoding data the GRIB record, local extension DWD
  INTEGER     :: localVersionNumber    !< version number, local DWD

  !> write COSMO field to grib record
  INTERFACE write_extpar_cosmo_field_grib
    MODULE PROCEDURE write_extpar_cosmo_real_field_grib
    MODULE PROCEDURE write_extpar_cosmo_int_field_grib
  END INTERFACE write_extpar_cosmo_field_grib

  !> write GME field to grib record
  INTERFACE write_extpar_gme_field_grib
    MODULE PROCEDURE write_extpar_gme_real_field_grib
    MODULE PROCEDURE write_extpar_gme_int_field_grib
  END INTERFACE write_extpar_gme_field_grib

  !> write ICON field to grib record
  INTERFACE write_extpar_ICON_field_grib
    MODULE PROCEDURE write_extpar_ICON_real_field_grib
    MODULE PROCEDURE write_extpar_ICON_int_field_grib
  END INTERFACE write_extpar_ICON_field_grib


  CONTAINS

  !> set definitions for rotated longitude latitude grid for the GRIB with GRIB_API
  !! the grib message should have been previously defined, pass the grib_id to this subroutine
  SUBROUTINE set_rotated_ll_grid_gds(gribid, cosmo_grid)
    USE mo_grid_structures, ONLY: rotated_lonlat_grid
    USE mo_utilities_extpar, ONLY: get_rot_spol_coor

    INTEGER, INTENT(IN) :: gribid !< id of grib message (GRIB_API)
    TYPE(rotated_lonlat_grid), INTENT(IN) :: cosmo_grid !< structure which contains the definition of the COSMO grid

    ! local variables
    REAL (KIND=wp) :: spollon
    REAL (KIND=wp) :: spollat
    INTEGER :: errorcode
    INTEGER :: editionNumber !< GRIB edition number

    CALL grib_set(gribid,'gridType','rotated_ll',errorcode)

    CALL get_rot_spol_coor(cosmo_grid%pollon,cosmo_grid%pollat,spollon,spollat)
    longitudeOfSouthernPoleInDegrees = spollon
    latitudeOfSouthernPoleInDegrees = spollat
    angleOfRotationInDegrees = cosmo_grid%polgam

    longitudeOfFirstGridPointInDegrees = cosmo_grid%startlon_rot
    
    latitudeOfFirstGridPointInDegrees = cosmo_grid%startlat_rot

    iDirectionIncrementInDegrees = cosmo_grid%dlon_rot
    jDirectionIncrementInDegrees = cosmo_grid%dlat_rot

    Ni = cosmo_grid%nlon_rot
    Nj = cosmo_grid%nlat_rot

    longitudeOfLastGridPointInDegrees = cosmo_grid%startlon_rot + (cosmo_grid%nlon_rot - 1) * cosmo_grid%dlon_rot
    latitudeOfLastGridPointInDegrees = cosmo_grid%startlat_rot + (cosmo_grid%nlat_rot -1) *  cosmo_grid%dlat_rot

    CALL grib_get(gribid,'editionNumber',editionNumber,errorcode)


    CALL grib_set(gribid,'longitudeOfSouthernPoleInDegrees',longitudeOfSouthernPoleInDegrees,errorcode)
    CALL grib_set(gribid,'latitudeOfSouthernPoleInDegrees',latitudeOfSouthernPoleInDegrees,errorcode)
    CALL grib_set(gribid,'angleOfRotationInDegrees',angleOfRotationInDegrees,errorcode)

    CALL grib_set(gribid,'longitudeOfFirstGridPointInDegrees',longitudeOfFirstGridPointInDegrees,errorcode)
    CALL grib_set(gribid,'latitudeOfFirstGridPointInDegrees',latitudeOfFirstGridPointInDegrees,errorcode)

    CALL grib_set(gribid,'Ni',Ni,errorcode)
    CALL grib_set(gribid,'Nj',Nj,errorcode)

    IF (editionNumber >= 2) THEN 
    CALL grib_set(gribid,'iDirectionIncrementInDegrees',iDirectionIncrementInDegrees,errorcode)
    CALL grib_set(gribid,'jDirectionIncrementInDegrees',jDirectionIncrementInDegrees,errorcode)
    END IF 

    CALL grib_set(gribid,'longitudeOfLastGridPointInDegrees',longitudeOfLastGridPointInDegrees,errorcode)
    CALL grib_set(gribid,'latitudeOfLastGridPointInDegrees',latitudeOfLastGridPointInDegrees,errorcode)


  END  SUBROUTINE set_rotated_ll_grid_gds

  
!  !> set date for an invariant field of the external parameters in a GRIB message with GRIB_API
!  !! the convention at DWD is to set the date for the invariant fields to
!  !! year 1, january 1, 00:00 hour
!  !! the grib message should have been previously defined, pass the grib_id to this subroutine
!  SUBROUTINE get_date_const_field_grib(gribid)
!    INTEGER, INTENT(IN)   :: gribid !< id of grib message (GRIB_API)
!    !local variables
!    INTEGER :: errorcode
!
!    INTEGER  :: cc
!    INTEGER  :: yy
!    INTEGER  :: mm
!    INTEGER  :: dd
!    INTEGER  :: hh
!    INTEGER  :: minute
!    INTEGER (KIND=i8)  :: dataDate
!    INTEGER (KIND=i8)  :: dataTime
!
!    cc = 0
!    yy = 1
!    mm = 1
!    dd = 1
!    hh = 0
!    minute=0
!    ! dataDate format ccyymmdd
!    ! dataTime format hhmm (hour minute)
!    dataDate = dd + (100*mm) + (10000*yy) + (1000000*cc)
!    dataTime = minute + (100*hh)
!    CALL grib_set(gribid,'dataDate ',dataDate)
!    CALL grib_set(gribid,'dataTime ',dataTime)
!
!  END  SUBROUTINE get_date_const_field_grib
!
!  !> set date for an monthly climatology field of the external parameters in a GRIB message with GRIB_API
!  !! the convention at DWD is to set the date for the invariant fields to
!  !! year 1111, month mm, day 11, 00:00 hour
!  !! the grib message should have been previously defined, pass the grib_id to this subroutine
!  SUBROUTINE set_date_mm_extpar_field_grib(gribid,mm)
!    INTEGER, INTENT(IN)   :: gribid !< id of grib message (GRIB_API)
!    INTEGER, INTENT(IN)   :: mm     !< month
!
!    !local variables
!    INTEGER :: errorcode
!
!    INTEGER  :: cc
!    INTEGER  :: yy
!    INTEGER  :: dd
!    INTEGER  :: hh
!    INTEGER  :: minute
!    INTEGER (KIND=i8)  :: dataDate
!    INTEGER (KIND=i8)  :: dataTime
!
!    cc = 11
!    yy = 11
!    dd = 11
!    hh = 0
!    minute=0
!    ! dataDate format ccyymmdd
!    ! dataTime format hhmm (hour minute)
!    dataDate = dd + (100*mm) + (10000*yy) + (1000000*cc)
!    dataTime = minute + (100*hh)
!
!    CALL grib_set(gribid,'dataDate ',dataDate)
!    CALL grib_set(gribid,'dataTime ',dataTime)
!
!  END  SUBROUTINE set_date_mm_extpar_field_grib



  !> set product defintion for a GRIB message with GRIB_API
  !! the grib message should have been previously defined, pass the gribid to this subroutine
  SUBROUTINE set_parameter_grib(gribid,shortName,dataDate,dataTime)

    INTEGER, INTENT(IN)   :: gribid !< id of grib message (GRIB_API)
    CHARACTER (LEN=*), INTENT(IN) :: shortName 
    INTEGER (KIND=i8), INTENT(IN)  :: dataDate  
!< date, for edition independent use of GRIB_API dataDate as Integer in the format ccyymmdd
    INTEGER (KIND=i8), INTENT(IN)  :: dataTime  !< time, for edition independent use GRIB_API dataTime in the format hhmm


    !local variables
    INTEGER :: errorcode
    INTEGER :: centre !< Identification of originating/generating centre,
    INTEGER :: editionNumber !< GRIB edition number


    CHARACTER(len=10) :: ydate
    CHARACTER(len=10) :: ytime
    INTEGER  :: cc
    INTEGER  :: yy
    INTEGER  :: mm
    INTEGER  :: dd
    INTEGER  :: hh
    INTEGER  :: minute

    CALL grib_get(gribid,'centre',centre,errorcode)
    CALL grib_get(gribid,'editionNumber',editionNumber,errorcode)

    ! define precision for the GRIB
    bitsPerValue=16
    CALL grib_set(gribid,"bitsPerValue", bitsPerValue)
    CALL grib_set(gribid,'resolutionAndComponentFlags',0)
     ! set stepType for GRIB1 and GRIB2, as the GRIB sample may not have the correct value
     ! set timeRangeIndicator for GRIB1 explicit
     ! the dataDate values are conventions for external parameter fields 
     ! \TODO maybe in GRIB2 these settings have to be improved
    CALL grib_set(gribid,"stepType","instant")
    IF (editionNumber == 1) THEN
      CALL  grib_set(gribid,"timeRangeIndicator",0) ! GRIB1 setting
      call grib_set (gribid,'stepRange',0)
      call grib_set (gribid,'startStep',0)
      call grib_set (gribid,'endStep',0)
      call grib_set (gribid,'PVPresent',0)
      call grib_set (gribid,'numberOfVerticalCoordinateValues',0)
    ENDIF
    IF ((dataDate >= 11110111).AND.(dataDate <= 11111211)) THEN
!      CALL grib_set(gribid,"stepType","avg")
      IF (editionNumber == 1) THEN
        CALL  grib_set(gribid,"timeRangeIndicator",3) ! GRIB1 setting
      ELSE
!        print *,'FB, debug, mo_io_grib_api grib_set stepType=avg'
        CALL grib_set(gribid,"stepType","avg")
      ENDIF
    ENDIF

    CALL grib_set(gribid,'dataDate ',dataDate)
    CALL grib_set(gribid,'dataTime ',dataTime)

!    PRINT *,'HA debug: shortName: ',TRIM(shortName)
    CALL grib_set(gribid,'shortName',TRIM(shortName),errorcode) 
    ! this is an "edition independent" setting of the parameter values 
    ! in the product defintion section, the GRIB1 or GRIB2 entries are 
    ! taken from the data of the *.def files in $GRIB_DEFINITION_PATH/grib1 
    ! or $GRIB_DEFINITION_PATH/grib1
    
    IF (editionNumber == 1) THEN
      IF (TRIM(shortName)=='DEPTH_LK') THEN
        CALL grib_set(gribid,"indicatorOfTypeOfLevel", 1) ! ground level type
        CALL grib_set(gribid,"indicatorOfParameter", 96)  ! element
        CALL grib_set(gribid,"table2Version", 201)        ! table
      ENDIF
      IF (TRIM(shortName)=='FR_LAKE') THEN
        CALL grib_set(gribid,"indicatorOfTypeOfLevel", 1) ! ground level type
        CALL grib_set(gribid,"indicatorOfParameter", 55)  ! element
        CALL grib_set(gribid,"table2Version", 202)        ! table
      ENDIF
    ELSEIF (editionNumber == 2) THEN
      IF (TRIM(shortName)=='DEPTH_LK') THEN
        CALL grib_set(gribid,"typeOfFirstFixedSurface", 1) ! ground level type
        CALL grib_set(gribid,"discipline", 1)       ! discipline 'Hydrological Products'
        CALL grib_set(gribid,"parameterCategory", 2)! parameterCategory 'inland water and sediment properties'
        CALL grib_set(gribid,"parameterNumber",0)   ! parameterNumber 'water depth'
      ENDIF
      IF (TRIM(shortName)=='FR_LAKE') THEN
        CALL grib_set(gribid,"typeOfFirstFixedSurface", 1) ! ground level type
        CALL grib_set(gribid,"discipline", 1)       ! discipline 'Hydrological Products'
        CALL grib_set(gribid,"parameterCategory", 2)! parameterCategory 'inland water and sediment properties'
        CALL grib_set(gribid,"parameterNumber",2)   ! parameterNumber 'water fraction'
      ENDIF
    ENDIF

     
    IF ((centre == dwd_id_grib ).AND.(editionNumber == 1)) THEN
      ! put data and time of output generation to DWD local section to the grib
      CALL DATE_AND_TIME(ydate,ytime)
      READ(ydate,'(4I2)') cc,yy,mm,dd
      READ(ytime,'(2I2)') hh, minute
      IF (cc == 20) THEN
        yy = yy + 100 ! in the GRIB, the year yy starts with 1900 for the "localDecodeDate"?
      ENDIF

      CALL grib_set(gribid,'localDecodeDateYear ',yy)
      CALL grib_set(gribid,'localDecodeDateMonth ',mm)
      CALL grib_set(gribid,'localDecodeDateDay ',dd)
      CALL grib_set(gribid,'localDecodeDateHour ',hh)
      CALL grib_set(gribid,'localDecodeDateMinute ',minute)
    ENDIF

     CALL grib_set(gribid,'resolutionAndComponentFlags',0)

  END SUBROUTINE set_parameter_grib

  !> write field from EXTPAR to a GRIB file for a COSMO field (rotated lon lat grid)
  !! the file should have been previously opende, pass the outfile_id to this subroutine
  SUBROUTINE write_extpar_cosmo_real_field_grib(outfile_id,grib_sample,cosmo_grid,extpar_buffer,shortName,dataDate,dataTime)
    USE mo_grid_structures, ONLY: rotated_lonlat_grid

    INTEGER, INTENT(IN) :: outfile_id !< id of the GRIB file
    CHARACTER (len=*), INTENT(IN) :: grib_sample  !< name for grib sample  (sample to be found in $GRIB_SAMPLES_PATH)

    TYPE(rotated_lonlat_grid), INTENT(IN)  :: cosmo_grid !< structure which contains the definition of the COSMO grid
    REAL (KIND=wp), INTENT(IN)             :: extpar_buffer(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1:1)
!< field to write out to GRIB file with outfile_id
    CHARACTER (LEN=*), INTENT(IN) :: shortName !< shortName parameter of the field
    INTEGER (KIND=8), INTENT(IN)  :: dataDate 
 !< date, for edition independent use of GRIB_API dataDate as Integer in the format ccyymmdd
    INTEGER (KIND=8), INTENT(IN)  :: dataTime  !< time, for edition independent use GRIB_API dataTime in the format hhmm


    ! local variables
    INTEGER :: gribid_in !< id of grib message (GRIB_API)
    INTEGER :: gribid_dest !< id of grib message (GRIB_API)

    REAL    :: ds(1:cosmo_grid%nlon_rot*cosmo_grid%nlat_rot) !< working array
    INTEGER :: errorcode
    INTEGER :: ind
    INTEGER :: i,j

    ! create a new grib message from sample, (sample to be found in $GRIB_SAMPLES_PATH)
    CALL grib_new_from_samples(gribid_in, TRIM(grib_sample))
    CALL grib_clone(gribid_in, gribid_dest, errorcode) ! copy the grib message to a new one
    CALL grib_release(gribid_in) ! free memory of first message

    CALL set_rotated_ll_grid_gds(gribid_dest,cosmo_grid) ! set gds values for rotated longitude latitude grid
    CALL set_parameter_grib(gribid_dest,shortName,dataDate,dataTime) ! set parameter values for gthe GRIB

    ind=0
    DO i=1,cosmo_grid%nlon_rot
    DO j=1,cosmo_grid%nlat_rot
      ind=(j-1) * cosmo_grid%nlon_rot + i
      ds(ind) = extpar_buffer(i,j,1)  ! put data to 1D array before putting it to the GRIB record
    ENDDO
    ENDDO

    ! ds = RESHAPE(extpar_buffer,(/ SIZE(ds) /))

    CALL grib_set(gribid_dest,'values',ds,errorcode) ! put data to GRIB message
    CALL grib_write(gribid_dest,outfile_id,errorcode) ! write out GRIB message to file
    CALL grib_release(gribid_dest) ! free memory of grib message

  END  SUBROUTINE write_extpar_cosmo_real_field_grib


  !> write field from EXTPAR to a GRIB file for a COSMO field (rotated lon lat grid)
  !! the file should have been previously opende, pass the outfile_id to this subroutine
  SUBROUTINE write_extpar_cosmo_int_field_grib(outfile_id,grib_sample,cosmo_grid,extpar_buffer,shortName,dataDate,dataTime)
    USE mo_grid_structures, ONLY: rotated_lonlat_grid

    INTEGER, INTENT(IN) :: outfile_id !< id of the GRIB file
    CHARACTER (len=*), INTENT(IN) :: grib_sample  !< name for grib sample  (sample to be found in $GRIB_SAMPLES_PATH)

    TYPE(rotated_lonlat_grid), INTENT(IN)  :: cosmo_grid !< structure which contains the definition of the COSMO grid
    INTEGER(KIND=i4), INTENT(IN)             :: extpar_buffer(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1:1) 
!< field to write out to GRIB file with outfile_id
    CHARACTER (LEN=*), INTENT(IN) :: shortName !< shortName parameter of the field
    INTEGER (KIND=8), INTENT(IN)  :: dataDate  
!< date, for edition independent use of GRIB_API dataDate as Integer in the format ccyymmdd
    INTEGER (KIND=8), INTENT(IN)  :: dataTime  !< time, for edition independent use GRIB_API dataTime in the format hhmm



    ! local variables
    INTEGER :: gribid_in !< id of grib message (GRIB_API)
    INTEGER :: gribid_dest !< id of grib message (GRIB_API)

    INTEGER :: ds(1:cosmo_grid%nlon_rot*cosmo_grid%nlat_rot) !< working array
    INTEGER :: errorcode
    INTEGER :: ind
    INTEGER :: i,j

    ! create a new grib message from sample, (sample to be found in $GRIB_SAMPLES_PATH)
    CALL grib_new_from_samples(gribid_in, TRIM(grib_sample))
    CALL grib_clone(gribid_in, gribid_dest, errorcode) ! copy the grib message to a new one
    CALL grib_release(gribid_in) ! free memory of first message

    CALL set_rotated_ll_grid_gds(gribid_dest,cosmo_grid) ! set gds values for rotated longitude latitude grid
    CALL set_parameter_grib(gribid_dest,shortName,dataDate,dataTime) ! set parameter values for gthe GRIB


    ind=0
    DO i=1,cosmo_grid%nlon_rot
    DO j=1,cosmo_grid%nlat_rot
      ind=(j-1) * cosmo_grid%nlon_rot + i
      ds(ind) = extpar_buffer(i,j,1)  ! put data to 1D array before putting it to the GRIB record
    ENDDO
    ENDDO

    ! ds = RESHAPE(extpar_buffer,(/ SIZE(ds) /))


    CALL grib_set(gribid_dest,'values',ds,errorcode) ! put data to GRIB message
    CALL grib_write(gribid_dest,outfile_id,errorcode) ! write out GRIB message to file
    CALL grib_release(gribid_dest) ! free memory of grib message


  END  SUBROUTINE write_extpar_cosmo_int_field_grib

  
  !> write field from EXTPAR to a GRIB file for a COSMO field (rotated lon lat grid)
  !! the file should have been previously opende, pass the outfile_id to this subroutine
  SUBROUTINE write_extpar_gme_real_field_grib(outfile_id,grib_sample,gme_grid,extpar_buffer,shortName,dataDate,dataTime)
    USE mo_grid_structures, ONLY: gme_triangular_grid

    INTEGER, INTENT(IN) :: outfile_id !< id of the GRIB file
    CHARACTER (len=*), INTENT(IN) :: grib_sample  !< name for grib sample  (sample to be found in $GRIB_SAMPLES_PATH)

    TYPE(gme_triangular_grid), INTENT(IN)  :: gme_grid !< structure which contains the definition of the GME grid
    REAL (KIND=wp), INTENT(IN)             :: extpar_buffer(1:gme_grid%nip1,1:gme_grid%nip1,1:gme_grid%nd) 
!< field to write out to GRIB file with outfile_id
    CHARACTER (LEN=*), INTENT(IN) :: shortName !< shortName parameter of the field
    INTEGER (KIND=8), INTENT(IN)  :: dataDate  
!< date, for edition independent use of GRIB_API dataDate as Integer in the format ccyymmdd
    INTEGER (KIND=8), INTENT(IN)  :: dataTime  !< time, for edition independent use GRIB_API dataTime in the format hhmm


    ! local variables
    INTEGER :: gribid_in !< id of grib message (GRIB_API)
    INTEGER :: gribid_dest !< id of grib message (GRIB_API)

    REAL    :: ds(1:gme_grid%nip1*gme_grid%nip1*gme_grid%nd) !< working array
    INTEGER :: errorcode
    INTEGER :: ind
    INTEGER :: ie,je,ke
    INTEGER :: nip1

    INTEGER :: i,j,k

    ! create a new grib message from sample, (sample to be found in $GRIB_SAMPLES_PATH)
    CALL grib_new_from_samples(gribid_in,TRIM(grib_sample))
    CALL grib_clone(gribid_in, gribid_dest, errorcode) ! copy the grib message to a new one
    CALL grib_release(gribid_in) ! free memory of first message

    ! set GME resulution and typ number
    CALL set_gme_grid_def(gribid_dest, gme_grid)
    CALL set_parameter_grib(gribid_dest,shortName,dataDate,dataTime) ! set parameter values for gthe GRIB

    ! ds = RESHAPE(extpar_buffer,(/ SIZE(ds) /))
    ind=0
    nip1=gme_grid%nip1
    DO ke=1,gme_grid%nd
    DO je=1,gme_grid%nip1
    DO ie=1,gme_grid%nip1
      ind=1 + (ie-1) + (je-1)*nip1 + (ke-1)*nip1*nip1
      ds(ind) = extpar_buffer(ie,je,ke)  ! put data to 1D array before putting it to the GRIB record
    ENDDO
    ENDDO
    ENDDO

    CALL grib_set(gribid_dest,'values',ds,errorcode) ! put data to GRIB message
    CALL grib_write(gribid_dest,outfile_id,errorcode) ! write out GRIB message to file
    CALL grib_release(gribid_dest) ! free memory of grib message


  END  SUBROUTINE write_extpar_gme_real_field_grib


  !> write field from EXTPAR to a GRIB file for a COSMO field (rotated lon lat grid)
  !! the file should have been previously opende, pass the outfile_id to this subroutine
  SUBROUTINE write_extpar_gme_int_field_grib(outfile_id,grib_sample,gme_grid,extpar_buffer,shortName,dataDate,dataTime)
    USE mo_grid_structures, ONLY: gme_triangular_grid

    INTEGER, INTENT(IN) :: outfile_id !< id of the GRIB file
    CHARACTER (len=*), INTENT(IN) :: grib_sample  !< name for grib sample  (sample to be found in $GRIB_SAMPLES_PATH)

    TYPE(gme_triangular_grid), INTENT(IN)  :: gme_grid !< structure which contains the definition of the GME grid
    INTEGER (KIND=i4), INTENT(IN)             :: extpar_buffer(1:gme_grid%nip1,1:gme_grid%nip1,1:gme_grid%nd) 
!< field to write out to GRIB file with outfile_id
    CHARACTER (LEN=*), INTENT(IN) :: shortName !< shortName parameter of the field
    INTEGER (KIND=8), INTENT(IN)  :: dataDate  
!< date, for edition independent use of GRIB_API dataDate as Integer in the format ccyymmdd
    INTEGER (KIND=8), INTENT(IN)  :: dataTime  !< time, for edition independent use GRIB_API dataTime in the format hhmm

    ! local variables
    INTEGER :: gribid_in !< id of grib message (GRIB_API)
    INTEGER :: gribid_dest !< id of grib message (GRIB_API)

    INTEGER :: ds(1:gme_grid%nip1*gme_grid%nip1*gme_grid%nd) !< working array
    INTEGER :: errorcode
    INTEGER :: ind
    INTEGER :: ie,je,ke
    INTEGER :: nip1

    ! create a new grib message from sample, (sample to be found in $GRIB_SAMPLES_PATH)
    CALL grib_new_from_samples(gribid_in, TRIM(grib_sample))
    CALL grib_clone(gribid_in, gribid_dest, errorcode) ! copy the grib message to a new one
    CALL grib_release(gribid_in) ! free memory of first message

    ! set GME resulution and typ number
    CALL set_gme_grid_def(gribid_dest, gme_grid)
    CALL set_parameter_grib(gribid_dest,shortName,dataDate,dataTime) ! set parameter values for gthe GRIB

    ! ds = RESHAPE(extpar_buffer,(/ SIZE(ds) /))
    ind=0
    nip1=gme_grid%nip1
    DO ke=1,gme_grid%nd
    DO je=1,gme_grid%nip1
    DO ie=1,gme_grid%nip1
      ind=1 + (ie-1) + (je-1)*nip1 + (ke-1)*nip1*nip1
      ds(ind) = extpar_buffer(ie,je,ke)  ! put data to 1D array before putting it to the GRIB record
    ENDDO
    ENDDO
    ENDDO

    CALL grib_set(gribid_dest,'values',ds,errorcode) ! put data to GRIB message
    CALL grib_write(gribid_dest,outfile_id,errorcode) ! write out GRIB message to file
    CALL grib_release(gribid_dest) ! free memory of grib message


  END  SUBROUTINE write_extpar_gme_int_field_grib

   !> set definitions for GME grid definition with GRIB_API
  !! the grib message should have been previously defined, pass the grib_id to this subroutine
  SUBROUTINE set_gme_grid_def(gribid, gme_grid)
    USE mo_grid_structures, ONLY: gme_triangular_grid

    INTEGER, INTENT(IN) :: gribid !< id of grib message (GRIB_API)
    TYPE(gme_triangular_grid), INTENT(IN) :: gme_grid !< structure which contains the definition of the COSMO grid

    ! local variables
    INTEGER :: errorcode
    INTEGER :: generatingProcessIdentifier ! typparameter
    INTEGER :: typnr !< GME database typ number
    INTEGER :: Ni ! resolution parameter

    CALL grib_set(gribid,'gridType','triangular_grid',errorcode)

    Ni = gme_grid%ni
    CALL get_typnr_gme(Ni,typnr)
    generatingProcessIdentifier = typnr

    CALL grib_set(gribid,'generatingProcessIdentifier',generatingProcessIdentifier,errorcode)
    CALL grib_set(gribid,'Ni',Ni,errorcode)


  END  SUBROUTINE set_gme_grid_def

  !> get typnr of GME database
  SUBROUTINE get_typnr_gme(ni,typnr)
  INTEGER :: Ni    !< resolution parameter
  INTEGER :: typnr !< GME database typ number
       SELECT CASE (Ni)
           CASE (16)
             typnr = 171
           CASE (32)
             typnr = 141
           CASE (48)
             typnr = 143
           CASE (64)
             typnr = 145
           CASE (96)
             typnr = 147
           CASE (128)
             typnr = 149
           CASE (192)
             typnr = 173
           CASE (256)
             typnr = 175
           CASE (384)
             typnr = 206
           CASE DEFAULT
             typnr = 999
         END SELECT

  END SUBROUTINE


  SUBROUTINE write_extpar_ICON_real_field_grib(outfile_id,gribid_dest,icon_grid,&
                                           extpar_buffer,shortName,dataDate,dataTime)
    USE mo_grid_structures, ONLY: icosahedral_triangular_grid

    INTEGER, INTENT(IN) :: outfile_id   !< id of the GRIB file, for writing
    INTEGER, INTENT(IN) :: gribid_dest  !< id of the GRIB2 process

    TYPE(icosahedral_triangular_grid), INTENT(IN)  :: icon_grid !< structure which contains the definition of the ICON grid
    REAL (KIND=wp), INTENT(IN)             :: extpar_buffer(1:icon_grid%ncell,1,1) 
!< field to write out to GRIB file with outfile_id
    CHARACTER (LEN=*), INTENT(IN) :: shortName !< shortName parameter of the field
    INTEGER (KIND=8), INTENT(IN)  :: dataDate  
!< date, for edition independent use of GRIB_API dataDate as Integer in the format ccyymmdd
    INTEGER (KIND=8), INTENT(IN)  :: dataTime  !< time, for edition independent use GRIB_API dataTime in the format hhmm


    ! local variables
    INTEGER :: igribclone !< id of grib message (GRIB_API)


    REAL    :: ds(1:icon_grid%ncell) !< working array
    INTEGER :: errorcode
    INTEGER :: ind
    INTEGER :: ie,je,ke
    INTEGER :: nip1

    INTEGER :: i,j,k

    CALL grib_clone(gribid_dest,igribclone) ! clone sample before modifying it
    CALL grib_set(igribclone,'shortName',shortName,errorcode) 

    ds=0._wp
 
    DO k=1,1
    DO j=1,1
    DO i=1,icon_grid%ncell
       ds(i) =extpar_buffer(i,j,k)   ! put data to 1D array before putting it to the GRIB record
    ENDDO
    ENDDO
    ENDDO

!    CALL grib_set (igribclone,'dataDate',dataDate)
!    CALL grib_set (igribclone,'dataTime',dataTime)
    CALL grib_set(igribclone,'values',ds,errorcode) ! put data to GRIB message
    CALL grib_write(igribclone,outfile_id,errorcode) ! write out GRIB message to file
    CALL grib_release(igribclone) ! free memory of grib message

  END  SUBROUTINE write_extpar_ICON_real_field_grib

  SUBROUTINE write_extpar_ICON_int_field_grib(outfile_id,gribid_dest,icon_grid,&
                                           extpar_buffer,shortName,dataDate,dataTime)
    USE mo_grid_structures, ONLY: icosahedral_triangular_grid

    INTEGER, INTENT(IN) :: outfile_id   !< id of the GRIB file, for writing
    INTEGER, INTENT(IN) :: gribid_dest  !< id of the GRIB2 process

    TYPE(icosahedral_triangular_grid), INTENT(IN)  :: icon_grid !< structure which contains the definition of the ICON grid
    INTEGER (KIND=i4), INTENT(IN)             :: extpar_buffer(1:icon_grid%ncell,1,1) 
!< field to write out to GRIB file with outfile_id
    CHARACTER (LEN=*), INTENT(IN) :: shortName !< shortName parameter of the field
    INTEGER (KIND=8), INTENT(IN)  :: dataDate  
!< date, for edition independent use of GRIB_API dataDate as Integer in the format ccyymmdd
    INTEGER (KIND=8), INTENT(IN)  :: dataTime  !< time, for edition independent use GRIB_API dataTime in the format hhmm


    ! local variables
    INTEGER :: igribclone !< id of grib message (GRIB_API)


    INTEGER :: ds(1:icon_grid%ncell) !< working array
    INTEGER :: errorcode
    INTEGER :: ind
    INTEGER :: ie,je,ke
    INTEGER :: nip1

    INTEGER :: i,j,k

    CALL grib_clone(gribid_dest,igribclone) ! clone sample before modifying it
    CALL grib_set(igribclone,'shortName',shortName,errorcode) 
    CALL grib_set(igribclone,'typeOfLevel','surface')
    CALL grib_set(igribclone,'typeOfOriginalFieldValues',1)

    ds=0
 
    DO k=1,1
    DO j=1,1
    DO i=1,icon_grid%ncell
       ds(i) =extpar_buffer(i,j,k)   ! put data to 1D array before putting it to the GRIB record
    ENDDO
    ENDDO
    ENDDO

!    CALL grib_set (igribclone,'dataDate',dataDate)
!    CALL grib_set (igribclone,'dataTime',dataTime)
    CALL grib_set(igribclone,'values',ds,errorcode) ! put data to GRIB message
    CALL grib_write(igribclone,outfile_id,errorcode) ! write out GRIB message to file
    CALL grib_release(igribclone) ! free memory of grib message

  END  SUBROUTINE write_extpar_ICON_int_field_grib
!------------------------------------------------------------------------------


END MODULE mo_io_grib_api
