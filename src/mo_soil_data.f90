!+ Fortran module with data fields for soil data
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V2_0         1013/06/04 Martina Messmer
!  adaptations such that HWSD data set can be used (top- and subsoil)
!  Code received from Juergen Helmert
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module with data fields for soil data
!> \author Hermann Asensio
MODULE mo_soil_data

  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4

  USE mo_grid_structures,       ONLY: reg_lonlat_grid
                             
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: define_soiltype,               &
       &    allocate_raw_soil_fields,      &
       &    allocate_raw_deep_soil_fields, &
       &    dsmw_legend,                   & 
       &    soil_texslo,                   &
       &    soil_texslo_deep,              &
       &    dsmw_soil_unit,                &
       &    dsmw_deep_soil_unit,           &
       &    n_unit,                        &
       &    dsmw_grid,                     &
       &    lon_soil,                      &
       &    lat_soil

  PUBLIC :: undef_soiltype, default_soiltype, soiltype_ice, soiltype_water, no_data

  PUBLIC :: FAO_data, HWSD_data, HWSD_map
  PUBLIC :: soil_data
  PUBLIC :: deep_soil

  PUBLIC :: lon_full, lat_full


  !> Definition of Data Type to describe the legend for the FAO Digital Soil Map of the World
  TYPE :: dsmw_legend
       INTEGER (KIND=i4)    :: dsmw_code      !< Code number for FAO Digital Soil Map of the World
       ! secondary code legend: 
       ! water = 9000 
       ! glacier, ice = 9001
       ! rock = 9002 
       ! salt = 9003 
       ! histosols (e.g. peat) = 9004 
       ! dunes, shifting sands = 9005 
       ! no data = 9009 
       ! standard classification < 9000
       REAL (KIND=wp)      :: tex_coarse, &    !< Coarse textured part of the soil unit area
            &                 tex_medium, &    !< Medium textured part of the soil unit area
            &                 tex_fine, &      !< Fine textured part of the soil unit area
            &                 part_undefined, &!< Undefined part of the soil unit area
            &                 flat, &          !< Flat slope part of the soil unit area
            &                 hilly, &         !< Hilly slope part of the soil unit area
            &                 steep         !< Steep slope part of the soil unit area

  END TYPE dsmw_legend

  TYPE(dsmw_legend), ALLOCATABLE :: soil_texslo(:), &                  !< legend for DSMW with texture and slope information
       &                            soil_texslo_deep(:)             !< legend for DSMW with texture and slope information

  TYPE(reg_lonlat_grid)          :: dsmw_grid !< structure with defenition of the raw data grid for the FAO Digital Soil Map of the World

  INTEGER(KIND=i4)               :: n_unit   !< number of soil units

  ! FAO Digital Soil Map of the World, the values represent the soil unit number (see for legend in variable soil_texslo)
  INTEGER (KIND=i4), ALLOCATABLE :: dsmw_soil_unit(:,:)
  ! FAO Digital Soil Map of the World, the values represent the soil unit number (see for legend in variable soil_texslo) 
  INTEGER (KIND=i4), ALLOCATABLE :: dsmw_deep_soil_unit(:,:) 

  ! longitide coordinates of the soil grid in the geographical (lonlat) system, dimension (nlon_reg)
  REAL (KIND=wp), ALLOCATABLE    :: lon_soil(:)
  ! latitude coordinates of the soil grid in the geographical (lonlat) system, dimension (nlat_reg)     
  REAL (KIND=wp), ALLOCATABLE    :: lat_soil(:)
  ! longitide coordinates of the soil grid in the geographical (lonlat) system, dimension (nlon_reg)
  REAL (KIND=wp), ALLOCATABLE    :: lon_full(:)
  ! latitude coordinates of the soil grid in the geographical (lonlat) system, dimension (nlat_reg)     
  REAL (KIND=wp), ALLOCATABLE    :: lat_full(:)


  SAVE

  INTEGER (KIND=i4)            :: undef_soiltype, &    !< undefined value for soil type (ocean/no data)
       &                          default_soiltype, & !< default soil type loam (5)
       &                          soiltype_ice, &     !< soiltype for ice
       &                          soiltype_water, &   !< soiltype for water
       &                          no_data, &          !< no data flag for FAO and HWSD
       &                          soil_data

  INTEGER (KIND=i4), PARAMETER :: FAO_data = 1, &
       &                          HWSD_data = 2, &
       &                          HWSD_map = 3

  LOGICAL                      :: deep_soil

  CONTAINS

  SUBROUTINE define_soiltype(isoil_data, ldeep_soil, &
                             undef_soiltype,         &
                             default_soiltype,       &
                             soiltype_ice,           &
                             soiltype_water,         &
                             soil_data) 

    IMPLICIT NONE

    INTEGER(KIND=i4),  INTENT(IN)  :: isoil_data
    LOGICAL,           INTENT(IN)  :: ldeep_soil
    INTEGER (KIND=i4), INTENT(OUT) :: undef_soiltype, &
         &                            default_soiltype, &
         &                            soiltype_ice, &
         &                            soiltype_water, &
         &                            soil_data

    soil_data = isoil_data
    deep_soil = ldeep_soil

    SELECT CASE(isoil_data)
      CASE(FAO_data, HWSD_map)
        undef_soiltype   = 0
        default_soiltype = 5     !< default soil type loam (5)
        soiltype_ice     = 1     !< soiltype for ice 
        soiltype_water   = 9     !< soiltype for water
        no_data          = 9009  !< no data flag of FAO
      CASE(HWSD_data)
        undef_soiltype   = 0
        default_soiltype = 30      !< default soil type loam (30)
        soiltype_ice     = 132     !< soiltype for ice 
        soiltype_water   = 135     !< soiltype for water
        no_data          = -9999   !< no data flag of HWSD
    END SELECT

  END SUBROUTINE define_soiltype

  !---------------------------------------------------------------------------------------------------------------------------!

  !> allocate raw data fields
  SUBROUTINE allocate_raw_soil_fields(ncolumns,nrows,n_units)

    IMPLICIT NONE

    INTEGER(KIND=i4) , INTENT(IN) :: ncolumns, & !< number of columns
         &                           nrows, &    !< number of rows
         &                           n_units   !< number of soil units

    INTEGER(KIND=i4)              :: errorcode !< error status variable

    CALL logging%info('Enter routine: allocate_raw_soil_fields')

    ALLOCATE(dsmw_soil_unit(1:ncolumns,1:nrows), STAT=errorcode) ! allocate dsmw_soil_unit
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the field dsmw_soil_unit',__FILE__,__LINE__)
    dsmw_soil_unit = 0      ! _FillValue of the DSMW

    ALLOCATE(soil_texslo(1:n_units), STAT=errorcode) ! allocate soil_texslo
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array soil_texslo',__FILE__,__LINE__)

    soil_texslo(:)%dsmw_code = no_data ! no data flag of FAO
    soil_texslo(:)%tex_coarse = 0.
    soil_texslo(:)%tex_medium = 0.
    soil_texslo(:)%tex_fine = 0.
    soil_texslo(:)%part_undefined = 0.
    soil_texslo(:)%flat = 0.
    soil_texslo(:)%hilly = 0.
    soil_texslo(:)%steep = 0.

    ALLOCATE(lon_soil(1:ncolumns), STAT=errorcode) 
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the field lon_soil',__FILE__,__LINE__)
    lon_soil = 0. ! _FillValue of the DSMW
 
    ALLOCATE(lat_soil(1:nrows), STAT=errorcode) 
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the field lat_soil',__FILE__,__LINE__)
    lat_soil = 0. ! _FillValue of the DSMW

  END SUBROUTINE allocate_raw_soil_fields

  !------------------------------------------------------------------------------------------------

  SUBROUTINE allocate_raw_deep_soil_fields(ncolumns,nrows,n_units)
  IMPLICIT NONE
  INTEGER , INTENT(IN) :: ncolumns !< number of columns
  INTEGER , INTENT(IN) :: nrows    !< number of rows
  INTEGER , INTENT(IN) :: n_units   !< number of soil units


  INTEGER :: errorcode !< error status variable

  ALLOCATE(dsmw_deep_soil_unit(1:ncolumns,1:nrows), STAT=errorcode) ! allocate dsmw_deep_soil_unit
  IF(errorcode.NE.0) CALL logging%error('Cant allocate the field dsmw_deep_soil_unit')
  
  dsmw_deep_soil_unit = 0 ! _FillValue of the DSMW
  
  ALLOCATE(soil_texslo_deep(1:n_units), STAT=errorcode) ! allocate soil_texslo_deep
  IF(errorcode.NE.0) CALL logging%error('Cant allocate the array soil_texslo_deep')

  soil_texslo_deep(:)%dsmw_code = no_data ! no data flag of FAO
  soil_texslo_deep(:)%tex_coarse = 0.
  soil_texslo_deep(:)%tex_medium = 0.
  soil_texslo_deep(:)%tex_fine = 0.
  soil_texslo_deep(:)%part_undefined = 0.
  soil_texslo_deep(:)%flat = 0.
  soil_texslo_deep(:)%hilly = 0.
  soil_texslo_deep(:)%steep = 0.

    

END SUBROUTINE allocate_raw_deep_soil_fields



END MODULE mo_soil_data
