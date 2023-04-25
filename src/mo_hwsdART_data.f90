!+ Fortran module with data fields for soil data
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2014/04/24 Daniel Rieger
!  Initial release
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> \author Daniel Rieger
MODULE mo_hwsdART_data

!> kind parameters are defined in MODULE data_parameters
USE mo_kind, ONLY: wp, &
                   i4, &
                   i8

!> abort_extpar defined in MODULE utilities_extpar
USE mo_utilities_extpar, ONLY: abort_extpar

USE mo_GRID_structures, ONLY: reg_lonlat_grid
                           
IMPLICIT NONE

PRIVATE

PUBLIC :: define_hwsdARTtype,            &
          allocate_raw_hwsdART_fields,   &
          hwsdART_soil_unit,             &
          lon_hwsdART,                   &
          lat_hwsdART

PUBLIC :: undef_hwsdARTtype, default_hwsdARTtype, no_data
PUBLIC :: type_clay_heavy, type_silty_clay, type_clay_light, type_silty_clay_loam
PUBLIC :: type_clay_loam, type_silt, type_silt_loam, type_sandy_clay, type_loam
PUBLIC :: type_sandy_clay_loam, type_sandy_loam, type_loamy_sand, type_sand

PUBLIC :: hwsdART_data, hwsdART_grid

INTEGER (KIND=i4), ALLOCATABLE :: hwsdART_soil_unit(:,:) !< 


TYPE(reg_lonlat_grid) :: hwsdART_grid !< structure with defenition of the raw data grid for the hwsd Soil Map of the World

REAL (KIND=wp), ALLOCATABLE  :: lon_hwsdART(:)          !< longitide coordinates of the hwsdART grid in the geographical (lonlat) system, dimension (nlon_reg)
REAL (KIND=wp), ALLOCATABLE  :: lat_hwsdART(:)          !< latitude coordinates of the hwsdART grid in the geographical (lonlat) system, dimension (nlat_reg)


SAVE

INTEGER (KIND=i4) :: undef_hwsdARTtype    !< undefined value for soil type (ocean/no data)
INTEGER (KIND=i4) :: default_hwsdARTtype  !< default soil type loam (9)
! HWSD TYPES:
INTEGER (KIND=i4) :: type_clay_heavy       !< type for heavy clay
INTEGER (KIND=i4) :: type_silty_clay       !< type for silty clay
INTEGER (KIND=i4) :: type_clay_light       !< type for light clay
INTEGER (KIND=i4) :: type_silty_clay_loam  !< type for silty clay loam
INTEGER (KIND=i4) :: type_clay_loam        !< type for clay loam
INTEGER (KIND=i4) :: type_silt             !< type for silt
INTEGER (KIND=i4) :: type_silt_loam        !< type for silt loam
INTEGER (KIND=i4) :: type_sandy_clay       !< type for sandy clay
INTEGER (KIND=i4) :: type_loam             !< type for loam
INTEGER (KIND=i4) :: type_sandy_clay_loam  !< type for sandy clay loam
INTEGER (KIND=i4) :: type_sandy_loam       !< type for sandy loam
INTEGER (KIND=i4) :: type_loamy_sand       !< type for loamy sand
INTEGER (KIND=i4) :: type_sand             !< type for sand



INTEGER (KIND=i4) :: no_data          !< no data flag for FAO and HWSD

INTEGER(KIND=i4)  :: hwsdART_data

CONTAINS

  SUBROUTINE define_hwsdARTtype() 

    IMPLICIT NONE

    undef_hwsdARTtype     = 0
    default_hwsdARTtype   = 9
    type_clay_heavy       = 1
    type_silty_clay       = 2
    type_clay_light       = 3
    type_silty_clay_loam  = 4
    type_clay_loam        = 5
    type_silt             = 6
    type_silt_loam        = 7
    type_sandy_clay       = 8
    type_loam             = 9
    type_sandy_clay_loam  = 10
    type_sandy_loam       = 11
    type_loamy_sand       = 12
    type_sand             = 13

    no_data               = -1

  END SUBROUTINE define_hwsdARTtype

  !---------------------------------------------------------------------------------------------------------------------------!

  !> allocate raw data fields
  SUBROUTINE allocate_raw_hwsdART_fields(ncolumns,nrows)
  IMPLICIT NONE
  INTEGER , INTENT(IN) :: ncolumns !< number of columns
  INTEGER , INTENT(IN) :: nrows    !< number of rows


  INTEGER :: errorcode !< error status variable


   ALLOCATE(hwsdART_soil_unit(1:ncolumns,1:nrows), STAT=errorcode) ! allocate hwsdART_soil_unit
      IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the field hwsdART_soil_unit')
      hwsdART_soil_unit = 0      ! _FillValue

   ALLOCATE(lon_hwsdART(1:ncolumns), STAT=errorcode) 
      IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the field lon_hwsdART')
    lon_hwsdART = 0. !
 
   ALLOCATE(lat_hwsdART(1:nrows), STAT=errorcode) 
      IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the field lat_hwsdART')
    lat_hwsdART = 0. !

  END  SUBROUTINE allocate_raw_hwsdART_fields


END MODULE mo_hwsdART_data
