!+ Fortran module with data fields for NDVI data
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
!> Fortran module with data fields for NDVI data
!> \author Hermann Asensio
MODULE mo_ndvi_data

  USE mo_logging
  USE mo_kind,                  ONLY: wp,i4

  USE mo_grid_structures,       ONLY: reg_lonlat_grid

  USE mo_ndvi_tg_fields,        ONLY: ndvi_field,    &
       &                              ndvi_max,      &
       &                              ndvi_field_mom,&
       &                              ndvi_ratio_mom
                             
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: allocate_raw_ndvi_fields, &
       &    ndvi_raw_data_grid, &
       &    ndvi_field_row_mom, &
       &    ndvi_field_row, &
       &    lon_ndvi, &
       &    lat_ndvi, &
       &    ntime_ndvi,&
       &    deallocate_ndvi_fields

  PUBLIC :: undef_ndvi, minimal_ndvi

  TYPE(reg_lonlat_grid)        :: ndvi_raw_data_grid
                           
  !< longitide coordinates of the soil grid in the geographical (lonlat) system, dimension (nlon_reg)
  REAL (KIND=wp), ALLOCATABLE  :: lon_ndvi(:), &
       &                          lat_ndvi(:), &
       &                          ndvi_field_row_mom(:,:), &      !< field for one row of ndvi data with monthly mean values
       &                          ndvi_field_row(:)      !< field for one row of ndvi data

  INTEGER (KIND=i4)            :: ntime_ndvi = 12 !< number of timesteps (12 for monthly mean values)

  REAL (KIND=wp)               :: undef_ndvi = 0.0, &  !< undefined value for NDVI data
       &                          minimal_ndvi = 0.09 !< minimal NDVI value bare soil value

  CONTAINS

  !> allocate raw data fields
  SUBROUTINE allocate_raw_ndvi_fields(ncolumns,nrows,nt)

    IMPLICIT NONE

    INTEGER(KINd=i4) , INTENT(IN) :: ncolumns, & !< number of columns
         &                           nrows, &    !< number of rows
         &                           nt !< number of timesteps (12 for monthly mean values)


    INTEGER(KIND=i4)              :: errorcode !< error status variable

    CALL logging%info('Enter routine: allocate_raw_ndvi_fields')

    ALLOCATE(ndvi_field_row_mom(1:ncolumns,1:nt), STAT=errorcode) 
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the field ndvi_field_row_mom',__FILE__,__LINE__)
    ndvi_field_row_mom = 0. 

    ALLOCATE(ndvi_field_row(1:ncolumns), STAT=errorcode) 
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the field ndvi_field_row',__FILE__,__LINE__)
    ndvi_field_row = 0. 

    ALLOCATE(lat_ndvi(1:nrows), STAT=errorcode) 
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the field lat_ndvi',__FILE__,__LINE__)
    lat_ndvi = 0. 

    ALLOCATE(lon_ndvi(1:ncolumns), STAT=errorcode) 
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the field lon_ndvi',__FILE__,__LINE__)
    lon_ndvi = 0. 

  END SUBROUTINE allocate_raw_ndvi_fields

  SUBROUTINE deallocate_ndvi_fields()


    IMPLICIT NONE

    INTEGER(KIND=i4) :: errorcode !< error status variable

    CALL logging%info('Enter routine: deallocate_ndvi_fields')

    DEALLOCATE(lat_ndvi, STAT=errorcode) 
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate the field lat_ndvi',__FILE__,__LINE__)
    DEALLOCATE(lon_ndvi, STAT=errorcode) 
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate the field lon_ndvi',__FILE__,__LINE__)
    DEALLOCATE(ndvi_field_row_mom, STAT=errorcode) 
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate the field ndvi_field_row_mom',__FILE__,__LINE__)
    DEALLOCATE(ndvi_field, STAT=errorcode) 
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate the field ndvi_field',__FILE__,__LINE__)
    DEALLOCATE(ndvi_max, STAT=errorcode) 
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate the field ndvi_max',__FILE__,__LINE__)
    DEALLOCATE(ndvi_field_mom, STAT=errorcode) 
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate the field ndvi_field_mom',__FILE__,__LINE__)
    DEALLOCATE(ndvi_ratio_mom, STAT=errorcode) 
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate the field ndvi_ratio_mom',__FILE__,__LINE__)
    DEALLOCATE(ndvi_field_row, STAT=errorcode) 
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate the field ndvi_field_row',__FILE__,__LINE__)

  END SUBROUTINE deallocate_ndvi_fields

END MODULE mo_ndvi_data
