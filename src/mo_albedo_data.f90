!+ Fortran module with data fields for Albedo data
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_8         2013/03/12 Frank Brenner
!  introduced MODIS albedo dataset(s) as new external parameter(s)
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module with data fields for Albedo data
!> \author Frank Brenner, Hermann Asensio
MODULE mo_albedo_data

!> kind parameters are defined in MODULE data_parameters
USE mo_kind, ONLY: wp, i4

!> abort_extpar defined in MODULE utilities_extpar
USE mo_utilities_extpar, ONLY: abort_extpar

USE mo_grid_structures, ONLY: reg_lonlat_grid
                           
IMPLICIT NONE

PRIVATE

PUBLIC    allocate_raw_alb_fields, &
          deallocate_raw_alb_fields, &
          alb_raw_data_grid, &
          alb_field_row, &
          lon_alb, &
          lat_alb, &
          ntime_alb, &
          ialb_type

PUBLIC :: undef_alb, minimal_alb

PUBLIC :: undef_alb_bs, minimal_alb_dry, maximal_alb_dry, minimal_alb_sat, &
          maximal_alb_sat

PUBLIC :: zalso, wso_min, wso_max, csalb, csalbw

PUBLIC :: allocate_alb_interp_fields,   &
          alb_interp_data

TYPE(reg_lonlat_grid) :: alb_raw_data_grid
                         
REAL (KIND=wp), ALLOCATABLE  :: lon_alb(:)  !< longitide coordinates, dimension (nlon_reg)
REAL (KIND=wp), ALLOCATABLE  :: lat_alb(:)  !< latitude coordinates, dimension (nlat_reg)

REAL (KIND=wp), ALLOCATABLE  :: alb_field_row(:)      !< field for one row of albedo data

REAL (KIND=wp), ALLOCATABLE :: zalso(:,:)
REAL (KIND=wp), ALLOCATABLE :: wso_min(:), wso_max(:)
REAL (KIND=wp), ALLOCATABLE :: csalb(:), csalbw(:)


INTEGER (KIND=i4) :: ntime_alb = 12 !< number of timesteps (12 for monthly mean values)
INTEGER (KIND=i4) :: ialb_type = 1 ! two possible values 
                                         ! 1 = create background albedo
                                         ! 2 = create soil albedo set

REAL (KIND=wp) :: undef_alb = 0.0  !< undefined value for ALB data
REAL (KIND=wp) :: minimal_alb = 0.07 !< minimal ALB value open sea
REAL (KIND=wp) :: undef_alb_bs = -1.E20  !< undefined value for bare soil albedo
REAL (KIND=wp) :: minimal_alb_dry = 0.1224 !< minimal value for dry soil
REAL (KIND=wp) :: maximal_alb_dry = 0.4925 !< maximum value for dry soil
REAL (KIND=wp) :: minimal_alb_sat = 0.0612 !< minimal value for saturated soil
REAL (KIND=wp) :: maximal_alb_sat = 0.3825 !< maximum value for saturated soil


CONTAINS

  !> allocate raw data fields
  SUBROUTINE allocate_raw_alb_fields(ncolumns,nrows)
  IMPLICIT NONE
  INTEGER , INTENT(IN) :: ncolumns !< number of columns
  INTEGER , INTENT(IN) :: nrows    !< number of rows

  INTEGER :: errorcode !< error status variable


  ALLOCATE(alb_field_row(1:ncolumns), STAT=errorcode) 
  IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the field alb_field_row')
  alb_field_row = 0. 

  ALLOCATE(lat_alb(1:nrows), STAT=errorcode) 
  IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the field lat_alb')
  lat_alb = 0. 

  ALLOCATE(lon_alb(1:ncolumns), STAT=errorcode) 
  IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the field lon_alb')
  lon_alb = 0. 

  END  SUBROUTINE allocate_raw_alb_fields

  SUBROUTINE allocate_alb_interp_fields(nt)

  IMPLICIT NONE
  INTEGER (KIND=i4), INTENT(IN) :: nt !< number of timesteps (12 for monthly mean values)

  ALLOCATE(zalso(0:10,1:nt))
  ALLOCATE(wso_min(0:10))
  ALLOCATE(wso_max(0:10))
  ALLOCATE(csalb(0:10))
  ALLOCATE(csalbw(0:10))

  END SUBROUTINE allocate_alb_interp_fields

  SUBROUTINE alb_interp_data()

     wso_min(1) = 0
     wso_min(2) = 0
     wso_min(3) = 0.012
     wso_min(4) = 0.03
     wso_min(5) = 0.035
     wso_min(6) = 0.06
     wso_min(7) = 0.065
     wso_min(8) = 0.098
     wso_min(9) = 0.0

     wso_max(1) = 0
     wso_max(2) = 0
     wso_max(3) = 0.364
     wso_max(4) = 0.445
     wso_max(5) = 0.455
     wso_max(6) = 0.475
     wso_max(7) = 0.507
     wso_max(8) = 0.863
     wso_max(9) = 0

     csalb(1) = 0.7 ! - 0.25*0     !0.7
     csalb(2) = 0.3 ! - 0.25*0     !0.3
     csalb(3) = 0.3 ! - 0.25*0.44  !0.19
     csalb(4) = 0.25! - 0.25*0.27  !0.1825
     csalb(5) = 0.25! - 0.25*0.24  !0.19
     csalb(6) = 0.25! - 0.25*0.23  !0.1925
     csalb(7) = 0.25! - 0.25*0.22  !0.195
     csalb(8) = 0.2 ! - 0.25*0.1   !0.175
     csalb(9) = 0.07! - 0.25*0     !0.07

     csalbw(1) = 0
     csalbw(2) = 0
     csalbw(3) = 0.44
     csalbw(4) = 0.27
     csalbw(5) = 0.24
     csalbw(6) = 0.23
     csalbw(7) = 0.22
     csalbw(8) = 0.1
     csalbw(9) = 0
    

  END SUBROUTINE alb_interp_data

  SUBROUTINE deallocate_raw_alb_fields


  IMPLICIT NONE 
    
  INTEGER :: errorcode !< error status variable
    
  DEALLOCATE(alb_field_row, STAT=errorcode) 
  IF(errorcode.NE.0) CALL abort_extpar('Cant deallocate array alb_field_row')
    
  DEALLOCATE(lat_alb, STAT=errorcode) 
  IF(errorcode.NE.0) CALL abort_extpar('Cant deallocate array lat_alb')
  DEALLOCATE(lon_alb, STAT=errorcode) 
  IF(errorcode.NE.0) CALL abort_extpar('Cant deallocate array lon_alb')

  END SUBROUTINE deallocate_raw_alb_fields
  
END MODULE mo_albedo_data
