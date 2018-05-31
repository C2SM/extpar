 !+ Fortran module for NDVI data on target grid for external parameters
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
!> Fortran module for albedo data on target grid for external parameters 
!> \author Frank Brenner, Hermann Asensio
MODULE mo_albedo_tg_fields

  !> kind parameters are defined in MODULE data_parameters
  USE mo_kind, ONLY: wp
  USE mo_kind, ONLY: i4

  !> abort_extpar defined in MODULE utilities_extpar
  USE mo_utilities_extpar, ONLY: abort_extpar

  USE mo_grid_structures, ONLY: target_grid_def


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: alb_dry, &
            alb_sat, &
            alb_field_mom, &
            alnid_field_mom, &
            aluvd_field_mom, &
            allocate_alb_target_fields, &
            deallocate_alb_target_fields, &
            alb_interpol



  REAL(KIND=wp), ALLOCATABLE  :: alb_field_mom(:,:,:,:) !< field for monthly mean albedo data (12 months)
  REAL(KIND=wp), ALLOCATABLE  :: alnid_field_mom(:,:,:,:)
  REAL(KIND=wp), ALLOCATABLE  :: aluvd_field_mom(:,:,:,:)

  REAL(KIND=wp), ALLOCATABLE :: alb_interpol(:,:,:,:) !<  field for interpolated albedo

  REAL(KIND=wp), ALLOCATABLE :: alb_dry(:,:,:) !< field for dry soil albedo
  REAL(KIND=wp), ALLOCATABLE :: alb_sat(:,:,:) !< field for saturated soil albedo


  CONTAINS

  !> allocate fields for albedo target data 
    SUBROUTINE allocate_alb_target_fields(tg,nt,raw_id)
    IMPLICIT NONE

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
    INTEGER (KIND=i4), INTENT(IN) :: nt !< number of timesteps (12 for monthly mean values)

    INTEGER (KIND=i4), INTENT(IN) :: raw_id !< type of albedo treatment

    INTEGER :: errorcode !< error status variable

    ALLOCATE (alb_field_mom(1:tg%ie,1:tg%je,1:tg%ke,1:nt), STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant allocate array alb_field_mom')
    alb_field_mom = 0.0

  !> the following fields are always used in the interface write_netcdf_cosmo_grid_extpar
  !> and must be allocated even if not used
    IF (raw_id == 2) THEN
      ALLOCATE (alb_dry(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
      IF(errorcode.NE.0) CALL abort_extpar('Cant allocate array alb_dry')
      alb_dry = 0.0

      ALLOCATE (alb_sat(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
      IF(errorcode.NE.0) CALL abort_extpar('Cant allocate array alb_sat')
      alb_sat = 0.0

      ALLOCATE (alnid_field_mom(0,0,0,0), STAT=errorcode)
      IF(errorcode.NE.0) CALL abort_extpar('Cant allocate array alnid_field_mom')
      alnid_field_mom = 0.0

      ALLOCATE (aluvd_field_mom(0,0,0,0), STAT=errorcode)
      IF(errorcode.NE.0) CALL abort_extpar('Cant allocate array aluvd_field_mom')
      aluvd_field_mom = 0.0

      ALLOCATE (alb_interpol(0,0,0,0), STAT=errorcode)
      IF(errorcode.NE.0) CALL abort_extpar('Cant allocate array alb_interpol')
      alb_interpol = 0.0

    ELSE  

      ALLOCATE (alb_dry(0,0,0), STAT=errorcode)
      IF(errorcode.NE.0) CALL abort_extpar('Cant allocate array alb_dry')
      alb_dry = 0.0

      ALLOCATE (alb_sat(0,0,0), STAT=errorcode)
      IF(errorcode.NE.0) CALL abort_extpar('Cant allocate array alb_sat')
      alb_sat = 0.0

      ALLOCATE (alnid_field_mom(1:tg%ie,1:tg%je,1:tg%ke,1:nt), STAT=errorcode)
      IF(errorcode.NE.0) CALL abort_extpar('Cant allocate array alnid_field_mom')
      alnid_field_mom = 0.0

      ALLOCATE (aluvd_field_mom(1:tg%ie,1:tg%je,1:tg%ke,1:nt), STAT=errorcode)
      IF(errorcode.NE.0) CALL abort_extpar('Cant allocate array aluvd_field_mom')
      aluvd_field_mom = 0.0

      ALLOCATE (alb_interpol(1:tg%ie,1:tg%je,1:tg%ke,1:nt), STAT=errorcode)
      IF(errorcode.NE.0) CALL abort_extpar('Cant allocate array alb_interpol')
      alb_interpol = 0.0
    ENDIF


    END SUBROUTINE allocate_alb_target_fields

  !> deallocate fields for albedo target data 
    SUBROUTINE deallocate_alb_target_fields(raw_id)
    IMPLICIT NONE

    INTEGER (KIND=i4), INTENT(IN) :: raw_id !< type of albedo treatment

    INTEGER :: errorcode !< error status variable

    DEALLOCATE (alb_field_mom, STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant deallocate array alb_field_mom')

    DEALLOCATE (alb_dry, STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant deallocate array alb_dry')

    DEALLOCATE (alb_sat, STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant deallocate array alb_sat')

    DEALLOCATE (alnid_field_mom, STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant deallocate array alnid_field_mom')

    DEALLOCATE (aluvd_field_mom, STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant deallocate array aluvd_field_mom')

    DEALLOCATE (alb_interpol, STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant deallocate arrayalb_interpol')


    END SUBROUTINE deallocate_alb_target_fields

END Module mo_albedo_tg_fields
