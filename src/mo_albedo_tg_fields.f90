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
  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4

  USE mo_grid_structures,       ONLY: target_grid_def


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: alb_dry, &
    &       alb_sat, &
    &       alb_field_mom, &
    &       alnid_field_mom, &
    &       aluvd_field_mom, &
    &       allocate_alb_target_fields, &
    &       deallocate_alb_target_fields, &
    &       alb_interpol



  REAL(KIND=wp), ALLOCATABLE  :: alb_field_mom(:,:,:,:), & !< field for monthly mean albedo data (12 months)
    &                            alnid_field_mom(:,:,:,:), &
    &                            aluvd_field_mom(:,:,:,:), &
    &                            alb_interpol(:,:,:,:), & !<  field for interpolated albedo
    &                            alb_dry(:,:,:), & !< field for dry soil albedo
    &                            alb_sat(:,:,:) !< field for saturated soil albedo


  CONTAINS

  !> allocate fields for albedo target data
  SUBROUTINE allocate_alb_target_fields(tg,nt,raw_id, l_use_array_cache)

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
    INTEGER (KIND=i4), INTENT(IN)     :: nt, & !< number of timesteps (12 for monthly mean values)
      &                                  raw_id !< type of albedo treatment
    LOGICAL, INTENT(in)               :: l_use_array_cache
    
    INTEGER(KIND=i4)                  :: errorcode !< error status variable

    CALL logging%info('Enter routine: allocate_alb_target_fields')

if (l_use_array_cache) then
   call allocate_cached('alb_field_mom', alb_field_mom, [tg%ie,tg%je,tg%ke,nt])
else
   allocate(alb_field_mom(tg%ie,tg%je,tg%ke,nt), stat=errorcode)
endif
    IF(errorcode.NE.0) CALL logging%error('Cant allocate array alb_field_mom',__FILE__,__LINE__)
    alb_field_mom = 0.0

  !> the following fields are always used in the interface write_netcdf_cosmo_grid_extpar
  !> and must be allocated even if not used
    IF (raw_id == 2) THEN
if (l_use_array_cache) then
   call allocate_cached('alb_dry', alb_dry, [tg%ie,tg%je,tg%ke])
else
   allocate(alb_dry(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
      IF(errorcode.NE.0) CALL logging%error('Cant allocate array alb_dry',__FILE__,__LINE__)
      alb_dry = 0.0

if (l_use_array_cache) then
   call allocate_cached('alb_sat', alb_sat, [tg%ie,tg%je,tg%ke])
else
   allocate(alb_sat(tg%ie,tg%je,tg%ke), stat=errorcode)
endif
      IF(errorcode.NE.0) CALL logging%error('Cant allocate array alb_sat',__FILE__,__LINE__)
      alb_sat = 0.0

if (l_use_array_cache) then
   call allocate_cached('alnid_field_mom', alnid_field_mom, [0,0,0,0])
else
   allocate(alnid_field_mom(0,0,0,0), stat=errorcode)
endif
      IF(errorcode.NE.0) CALL logging%error('Cant allocate array alnid_field_mom',__FILE__,__LINE__)
      alnid_field_mom = 0.0

if (l_use_array_cache) then
   call allocate_cached('aluvd_field_mom', aluvd_field_mom, [0,0,0,0])
else
   allocate(aluvd_field_mom(0,0,0,0), stat=errorcode)
endif
      IF(errorcode.NE.0) CALL logging%error('Cant allocate array aluvd_field_mom',__FILE__,__LINE__)
      aluvd_field_mom = 0.0

if (l_use_array_cache) then
   call allocate_cached('alb_interpol', alb_interpol, [0,0,0,0])
else
   allocate(alb_interpol(0,0,0,0), stat=errorcode)
endif
      IF(errorcode.NE.0) CALL logging%error('Cant allocate array alb_interpol',__FILE__,__LINE__)
      alb_interpol = 0.0

    ELSE

if (l_use_array_cache) then
   call allocate_cached('alb_dry', alb_dry, [0,0,0])
else
   allocate(alb_dry(0,0,0), stat=errorcode)
endif
      IF(errorcode.NE.0) CALL logging%error('Cant allocate array alb_dry',__FILE__,__LINE__)
      alb_dry = 0.0

if (l_use_array_cache) then
   call allocate_cached('alb_sat', alb_sat, [0,0,0])
else
   allocate(alb_sat(0,0,0), stat=errorcode)
endif
      IF(errorcode.NE.0) CALL logging%error('Cant allocate array alb_sat',__FILE__,__LINE__)
      alb_sat = 0.0

if (l_use_array_cache) then
   call allocate_cached('alnid_field_mom', alnid_field_mom, [tg%ie,tg%je,tg%ke,nt])
else
   allocate(alnid_field_mom(tg%ie,tg%je,tg%ke,nt), stat=errorcode)
endif
      IF(errorcode.NE.0) CALL logging%error('Cant allocate array alnid_field_mom',__FILE__,__LINE__)
      alnid_field_mom = 0.0

if (l_use_array_cache) then
   call allocate_cached('aluvd_field_mom', aluvd_field_mom, [tg%ie,tg%je,tg%ke,nt])
else
   allocate(aluvd_field_mom(tg%ie,tg%je,tg%ke,nt), stat=errorcode)
endif
      IF(errorcode.NE.0) CALL logging%error('Cant allocate array aluvd_field_mom',__FILE__,__LINE__)
      aluvd_field_mom = 0.0

if (l_use_array_cache) then
   call allocate_cached('alb_interpol', alb_interpol, [tg%ie,tg%je,tg%ke,nt])
else
   allocate(alb_interpol(tg%ie,tg%je,tg%ke,nt), stat=errorcode)
endif
      IF(errorcode.NE.0) CALL logging%error('Cant allocate array alb_interpol',__FILE__,__LINE__)
      alb_interpol = 0.0
    ENDIF

  END SUBROUTINE allocate_alb_target_fields

  !> deallocate fields for albedo target data
  SUBROUTINE deallocate_alb_target_fields()

    INTEGER (KIND=i4) :: errorcode !< error status variable

    CALL logging%info('Enter routine: deallocate_alb_target_fields')

    DEALLOCATE (alb_field_mom, STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate array alb_field_mom',__FILE__,__LINE__)

    DEALLOCATE (alb_dry, STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate array alb_dry',__FILE__,__LINE__)

    DEALLOCATE (alb_sat, STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate array alb_sat',__FILE__,__LINE__)

    DEALLOCATE (alnid_field_mom, STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate array alnid_field_mom',__FILE__,__LINE__)

    DEALLOCATE (aluvd_field_mom, STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate array aluvd_field_mom',__FILE__,__LINE__)

    DEALLOCATE (alb_interpol, STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate arrayalb_interpol',__FILE__,__LINE__)

  END SUBROUTINE deallocate_alb_target_fields

END Module mo_albedo_tg_fields
