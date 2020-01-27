!+ Fortran module with data fields for AHF data
!
!
! Description:
! Fortran module with data fields for AHF data
!
! Current Code Owner: DWD, <Name of person responsible for this code>
!    <smail, phone, fax and email>
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V2_4         2015-05-21 Hendrik Wouters
!  Initial release  based on mo_ndvi_data.f90 (V1_14)
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!=======================================================================
!> Fortran module with data fields for AHF data
!> \author Hermann Asensio
MODULE mo_ahf_data

  USE mo_logging
  USE mo_kind,                  ONLY: wp,i4

  USE mo_grid_structures,       ONLY: reg_lonlat_grid
  USE mo_ahf_tg_fields,         ONLY: ahf_field
                           
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: allocate_raw_ahf_fields, &
       &    ahf_raw_data_grid, &
       &    ahf_field_row, &
       &    lon_ahf, &
       &    lat_ahf, &
       &    deallocate_ahf_fields

  PUBLIC :: undef_ahf, minimal_ahf

  PUBLIC :: iahf_type !_br 14.04.16

  TYPE(reg_lonlat_grid)        :: ahf_raw_data_grid
                           
  REAL (KIND=wp), ALLOCATABLE  :: lon_ahf(:), &          !< longitide coordinates of the soil grid in the geographical (lonlat) syste &
       &                          lat_ahf(:), &          !< latitude coordinates of the soil grid in the geographical (lonlat) system &
       &                          ahf_field_row(:)      !< field for one row of ahf data

  REAL (KIND=wp)               :: undef_ahf = 0.0, &  !< undefined value for AHF data
       &                          minimal_ahf = 0.0 !< minimal AHF value bare soil value

  INTEGER (KIND=i4)            :: iahf_type = 1 !_br 14.04.16

  CONTAINS

  !> allocate raw data fields
  SUBROUTINE allocate_raw_ahf_fields(ncolumns,nrows)

    IMPLICIT NONE

    INTEGER , INTENT(IN) :: ncolumns, nrows

    INTEGER(KIND=i4)     :: errorcode !< error status variable

    CALL logging%info('Enter routine: allocate_raw_ahf_fields')

    ALLOCATE(ahf_field_row(1:ncolumns), STAT=errorcode) 
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the field ahf_field_row',__FILE__,__LINE__)
    ahf_field_row = 0. 

    ALLOCATE(lat_ahf(1:nrows), STAT=errorcode) 
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the field lat_ahf',__FILE__,__LINE__)
    lat_ahf = 0. 
    
    ALLOCATE(lon_ahf(1:ncolumns), STAT=errorcode) 
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the field lon_ahf',__FILE__,__LINE__)
    lon_ahf = 0. 

  END  SUBROUTINE allocate_raw_ahf_fields

  SUBROUTINE deallocate_ahf_fields()

    IMPLICIT NONE

    INTEGER :: errorcode !< error status variable

    CALL logging%info('Enter routine: deallocate_ahf_fields')

    DEALLOCATE(lat_ahf, STAT=errorcode) 
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate the field lat_ahf',__FILE__,__LINE__)
    DEALLOCATE(lon_ahf, STAT=errorcode) 
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate the field lon_ahf',__FILE__,__LINE__)
    DEALLOCATE(ahf_field, STAT=errorcode) 
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate the field ahf_field',__FILE__,__LINE__)
    DEALLOCATE(ahf_field_row, STAT=errorcode) 
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate the field ahf_field_row',__FILE__,__LINE__)

  END SUBROUTINE deallocate_ahf_fields

END MODULE mo_ahf_data
