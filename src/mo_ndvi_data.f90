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

!> kind parameters are defined in MODULE data_parameters
USE mo_kind, ONLY: wp, &
                   i4, &
                   i8

!> abort_extpar defined in MODULE utilities_extpar
USE mo_utilities_extpar, ONLY: abort_extpar

USE mo_grid_structures, ONLY: reg_lonlat_grid
                           
IMPLICIT NONE

PRIVATE

PUBLIC :: allocate_raw_ndvi_fields, &
          ndvi_raw_data_grid, &
          ndvi_field_row_mom, &
          ndvi_field_row, &
          lon_ndvi, &
          lat_ndvi, &
          ntime_ndvi,&
          deallocate_ndvi_fields

PUBLIC :: undef_ndvi, minimal_ndvi

TYPE(reg_lonlat_grid) :: ndvi_raw_data_grid
                         
!< longitide coordinates of the soil grid in the geographical (lonlat) system, dimension (nlon_reg)
REAL (KIND=wp), ALLOCATABLE  :: lon_ndvi(:)
 !< latitude coordinates of the soil grid in the geographical (lonlat) system, dimension (nlat_reg)    
REAL (KIND=wp), ALLOCATABLE  :: lat_ndvi(:)

REAL (KIND=wp), ALLOCATABLE  :: ndvi_field_row_mom(:,:)      !< field for one row of ndvi data with monthly mean values

REAL (KIND=wp), ALLOCATABLE  :: ndvi_field_row(:)      !< field for one row of ndvi data

INTEGER (KIND=i4) :: ntime_ndvi = 12 !< number of timesteps (12 for monthly mean values)

REAL (KIND=wp) :: undef_ndvi = 0.0  !< undefined value for NDVI data
REAL (KIND=wp) :: minimal_ndvi = 0.09 !< minimal NDVI value bare soil value



CONTAINS

  !> allocate raw data fields
  SUBROUTINE allocate_raw_ndvi_fields(ncolumns,nrows,nt)
  IMPLICIT NONE
  INTEGER , INTENT(IN) :: ncolumns !< number of columns
  INTEGER , INTENT(IN) :: nrows    !< number of rows
  INTEGER (KIND=i4), INTENT(IN) :: nt !< number of timesteps (12 for monthly mean values)


  INTEGER :: errorcode !< error status variable

   ALLOCATE(ndvi_field_row_mom(1:ncolumns,1:nt), STAT=errorcode) 
      IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the field ndvi_field_row_mom')
    ndvi_field_row_mom = 0. 


   ALLOCATE(ndvi_field_row(1:ncolumns), STAT=errorcode) 
      IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the field ndvi_field_row')
    ndvi_field_row = 0. 


   ALLOCATE(lat_ndvi(1:nrows), STAT=errorcode) 
      IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the field lat_ndvi')
    lat_ndvi = 0. 


   ALLOCATE(lon_ndvi(1:ncolumns), STAT=errorcode) 
      IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the field lon_ndvi')
    lon_ndvi = 0. 



  END  SUBROUTINE allocate_raw_ndvi_fields


  SUBROUTINE deallocate_ndvi_fields()

    USE mo_ndvi_tg_fields, ONLY: ndvi_field,    &
                        &        ndvi_max,      &
                        &        ndvi_field_mom,&
                        &        ndvi_ratio_mom

    IMPLICIT NONE

    INTEGER :: errorcode !< error status variable


    DEALLOCATE(lat_ndvi, STAT=errorcode) 
    IF(errorcode.NE.0) CALL abort_extpar('Cant deallocate the field lat_ndvi')
    DEALLOCATE(lon_ndvi, STAT=errorcode) 
    IF(errorcode.NE.0) CALL abort_extpar('Cant deallocate the field lon_ndvi')
    DEALLOCATE(ndvi_field_row_mom, STAT=errorcode) 
    IF(errorcode.NE.0) CALL abort_extpar('Cant deallocate the field ndvi_field_row_mom')
    DEALLOCATE(ndvi_field, STAT=errorcode) 
    IF(errorcode.NE.0) CALL abort_extpar('Cant deallocate the field ndvi_field')
    DEALLOCATE(ndvi_max, STAT=errorcode) 
    IF(errorcode.NE.0) CALL abort_extpar('Cant deallocate the field ndvi_max')
    DEALLOCATE(ndvi_field_mom, STAT=errorcode) 
    IF(errorcode.NE.0) CALL abort_extpar('Cant deallocate the field ndvi_field_mom')
    DEALLOCATE(ndvi_ratio_mom, STAT=errorcode) 
    IF(errorcode.NE.0) CALL abort_extpar('Cant deallocate the field ndvi_ratio_mom')
    DEALLOCATE(ndvi_field_row, STAT=errorcode) 
    IF(errorcode.NE.0) CALL abort_extpar('Cant deallocate the field ndvi_field_row')

    END SUBROUTINE deallocate_ndvi_fields

END MODULE mo_ndvi_data
