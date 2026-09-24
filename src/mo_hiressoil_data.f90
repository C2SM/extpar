MODULE mo_hiressoil_data

USE mo_logging
USE mo_kind, ONLY: wp, i4
USE mo_GRID_structures, ONLY: reg_lonlat_grid
                           
IMPLICIT NONE

PRIVATE

PUBLIC :: define_hiressoiltype,            &
     &    allocate_raw_hiressoil_fields,   &
     &    hiressoil_soil_unit,             &
     &    lon_hiressoil,                   &
     &    lat_hiressoil,                   &
     &    undef_hiressoiltype, default_hiressoiltype, no_data, &
     &    type_clay_heavy, type_silty_clay, type_clay_light, type_silty_clay_loam, &
     &    type_clay_loam, type_silt, type_silt_loam, type_sandy_clay, type_loam, &
     &    type_sandy_clay_loam, type_sandy_loam, type_loamy_sand, type_sand, &
     &    hiressoil_data, hiressoil_grid

INTEGER (KIND=i4), ALLOCATABLE :: hiressoil_soil_unit(:,:) !< 


TYPE(reg_lonlat_grid) :: hiressoil_grid !< structure with defenition of the raw data grid for the hwsd Soil Map of the World

REAL (KIND=wp), ALLOCATABLE  :: lon_hiressoil(:), & !< longitide coordinates of the hiressoil grid in the geographical (lonlat) system, dimension (nlon_reg)
                    &           lat_hiressoil(:)    !< latitude coordinates of the hiressoil grid in the geographical (lonlat) system, dimension (nlat_reg)


SAVE

INTEGER (KIND=i4) :: undef_hiressoiltype , &!< undefined value for soil type (ocean/no data)
    &     default_hiressoiltype, &      !< default soil type loam (9)
    &     type_clay_heavy, &          !< type for heavy clay
    &     type_silty_clay, &          !< type for silty clay
    &     type_clay_light, &          !< type for light clay
    &     type_silty_clay_loam, &     !< type for silty clay loam
    &     type_clay_loam, &           !< type for clay loam
    &     type_silt, &                !< type for silt
    &     type_silt_loam, &           !< type for silt loam
    &     type_sandy_clay, &          !< type for sandy clay
    &     type_loam, &                !< type for loam
    &     type_sandy_clay_loam, &     !< type for sandy clay loam
    &     type_sandy_loam, &          !< type for sandy loam
    &     type_loamy_sand, &          !< type for loamy sand
    &     type_sand, &                !< type for sand
    &     no_data, &                  !< no data flag for FAO and HWSD
    &     hiressoil_data

CONTAINS

  SUBROUTINE define_hiressoiltype() 

    IMPLICIT NONE

    undef_hiressoiltype     = 0
    default_hiressoiltype   = 9
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

  END SUBROUTINE define_hiressoiltype

  !---------------------------------------------------------------------------------------------------------------------------!

  !> allocate raw data fields
  SUBROUTINE allocate_raw_hiressoil_fields(ncolumns,nrows)
  IMPLICIT NONE
  INTEGER , INTENT(IN) :: ncolumns, & !< number of columns
       &                  nrows       !< number of rows

  INTEGER :: errorcode !< error status variable


   ALLOCATE(hiressoil_soil_unit(1:ncolumns,1:nrows), STAT=errorcode) ! allocate hiressoil_soil_unit
      IF(errorcode.NE.0)  CALL logging%error('Cant allocate the field hiressoil_soil_unit',__FILE__,__LINE__)
      hiressoil_soil_unit = 0      ! _FillValue

   ALLOCATE(lon_hiressoil(1:ncolumns), STAT=errorcode) 
      IF(errorcode.NE.0)  CALL logging%error('Cant allocate the field lon_hiressoil',__FILE__,__LINE__)
    lon_hiressoil = 0. !
 
   ALLOCATE(lat_hiressoil(1:nrows), STAT=errorcode) 
      IF(errorcode.NE.0)  CALL logging%error('Cant allocate the field lat_hiressoil',__FILE__,__LINE__)
    lat_hiressoil = 0. !

  END  SUBROUTINE allocate_raw_hiressoil_fields


END MODULE mo_hiressoil_data
