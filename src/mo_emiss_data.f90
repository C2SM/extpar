!+ Fortran module with data fields for EMISS data
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module with data fields for EMISS data
!> \author Hermann Asensio
MODULE mo_emiss_data

  USE mo_logging
  USE mo_kind,                  ONLY: wp,i4

  USE mo_grid_structures,       ONLY: reg_lonlat_grid
                             
  USE mo_emiss_tg_fields,       ONLY: emiss_field,    &
       &                              emiss_max,      &
       &                              emiss_field_mom,&
       &                              emiss_ratio_mom

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: allocate_raw_emiss_fields, &
       &    emiss_raw_data_grid, &
       &    emiss_field_row_mom, &
       &    emiss_field_row, &
       &    lon_emiss, &
       &    lat_emiss, &
       &    ntime_emiss,&
       &    deallocate_emiss_fields

  PUBLIC :: undef_emiss, minimal_emiss

  TYPE(reg_lonlat_grid) :: emiss_raw_data_grid
                           
  !< longitide coordinates of the soil grid in the geographical (lonlat) system, dimension (nlon_reg)
  REAL (KIND=wp), ALLOCATABLE  :: lon_emiss(:), & 
       &                          lat_emiss(:), & 
       &                          emiss_field_row_mom(:,:), &       !< field for one row of emiss data with monthly mean values
       &                          emiss_field_row(:)      !< field for one row of emiss data

  INTEGER (KIND=i4)            :: ntime_emiss = 12 !< number of timesteps (12 for monthly mean values)

  REAL (KIND=wp)               :: undef_emiss = 0.0, &   !< undefined value for EMISS data
       &                          minimal_emiss = 0.0 !< minimal EMISS value 

  CONTAINS

  !> allocate raw data fields
  SUBROUTINE allocate_raw_emiss_fields(ncolumns,nrows,nt)

    IMPLICIT NONE

    INTEGER(KIND=i4) , INTENT(IN) :: ncolumns, &  !< number of columns
         &                           nrows, &     !< number of rows
         &                           nt !< number of timesteps (12 for monthly mean values)

    INTEGER(KIND=i4)              :: errorcode
    
    CALL logging%info('Enter routine: allocate_raw_emiss_fields')

    ALLOCATE(emiss_field_row_mom(1:ncolumns,1:nt), STAT=errorcode) 
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the field emiss_field_row_mom',__FILE__,__LINE__)
    emiss_field_row_mom = 0. 

    ALLOCATE(emiss_field_row(1:ncolumns), STAT=errorcode) 
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the field emiss_field_row',__FILE__,__LINE__)
    emiss_field_row = 0. 

    ALLOCATE(lat_emiss(1:nrows), STAT=errorcode) 
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the field lat_emiss',__FILE__,__LINE__)
    lat_emiss = 0. 

    ALLOCATE(lon_emiss(1:ncolumns), STAT=errorcode) 
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the field lon_emiss',__FILE__,__LINE__)
    lon_emiss = 0. 

  END  SUBROUTINE allocate_raw_emiss_fields

  SUBROUTINE deallocate_emiss_fields()

    IMPLICIT NONE

    INTEGER(KIND=i4) :: errorcode !< error status variable

    CALL logging%info('Enter routine: deallocate_emiss_fields')

    DEALLOCATE(lat_emiss, STAT=errorcode) 
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate the field lat_emiss',__FILE__,__LINE__)
    DEALLOCATE(lon_emiss, STAT=errorcode) 
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate the field lon_emiss',__FILE__,__LINE__)
    DEALLOCATE(emiss_field_row_mom, STAT=errorcode) 
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate the field emiss_field_row_mom',__FILE__,__LINE__)
    DEALLOCATE(emiss_field, STAT=errorcode) 
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate the field emiss_field',__FILE__,__LINE__)
    DEALLOCATE(emiss_max, STAT=errorcode) 
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate the field emiss_max',__FILE__,__LINE__)
    DEALLOCATE(emiss_field_mom, STAT=errorcode) 
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate the field emiss_field_mom',__FILE__,__LINE__)
    DEALLOCATE(emiss_ratio_mom, STAT=errorcode) 
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate the field emiss_ratio_mom',__FILE__,__LINE__)
    DEALLOCATE(emiss_field_row, STAT=errorcode) 
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate the field emiss_field_row',__FILE__,__LINE__)

    END SUBROUTINE deallocate_emiss_fields

END MODULE mo_emiss_data
