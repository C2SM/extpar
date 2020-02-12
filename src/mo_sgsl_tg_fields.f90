!+ Fortran module for subgrid-scale slope data on target grid
!  for external parameters 
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V4_0         2016/07/26 Daniel Luethi
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module for subgrid-scale slope data on target grid 
!> for external parameters
!> \author Daniel Luethi
MODULE mo_sgsl_tg_fields

  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4

  USE mo_grid_structures,       ONLY: target_grid_def

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: sgsl, &
    &        allocate_sgsl_target_fields
            

  PUBLIC ::   add_parameters_domain, &
    &          vertex_param, &
    &          allocate_additional_sgsl_param
       
  REAL(KIND=wp), ALLOCATABLE        :: sgsl(:,:,:)  !< field with subgrid-scale slopes
  !> data structure for parameters on vertices of Icon grid

  TYPE add_parameters_domain
     REAL(KIND=wp), ALLOCATABLE     :: sgsl_vert(:,:,:)   !< height on vertex
     INTEGER (KIND=i4), ALLOCATABLE :: npixel_vert(:,:,:) !< number of raw data pixel corresponding to vertex
  END TYPE add_parameters_domain

  TYPE(add_parameters_domain)       :: vertex_param  !< additional external parameters for ICON domain



  CONTAINS

  !> allocate fields for GLOBE target data 
  SUBROUTINE allocate_sgsl_target_fields(tg)

    IMPLICIT NONE

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
    INTEGER(KIND=i4)                  :: errorcode !< error status variable
      
    CALL logging%info('Enter routine: allocate_sgsl_target_fields')

    ALLOCATE (sgsl(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array sgsl',__FILE__,__LINE__)
    sgsl = 0.0

  END SUBROUTINE allocate_sgsl_target_fields

  !> allocate additional parameters which correspond to the vertex
  !!
  !! the target grid has the dimension nvertex
  !! for future developments (optimizations, other code structure) the target grid is 
  !! defined as a 3-dimensional matrix, but the dimension are set to (nvertex,1,1) in this case
  SUBROUTINE allocate_additional_sgsl_param(nvertex)

    IMPLICIT NONE

    INTEGER(KIND=i4), INTENT(IN) :: nvertex  !< number of vertices in target domains

    INTEGER(KIND=i4), PARAMETER  :: je = 1, & !< dummy for dimension selection
         &                          ke = 1 !< dummy for dimension selection

    INTEGER(KIND=i4)             :: errorcode !< error status variable

    CALL logging%info('Enter routine: allocate_additional_sgsl_param')

    ALLOCATE(vertex_param%sgsl_vert(1:nvertex,1:je,1:ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the vertex_param%sgsl_vert(nvertex,je,ke',__FILE__,__LINE__)
    vertex_param%sgsl_vert = 0.0

    ALLOCATE(vertex_param%npixel_vert(1:nvertex,1:je,1:ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the vertex_param%npixel_vert(nvertex,je,ke)',__FILE__,__LINE__)
    vertex_param%npixel_vert = 0

  END SUBROUTINE allocate_additional_sgsl_param

END Module mo_sgsl_tg_fields
