MODULE mo_python_data

  USE mo_kind,                  ONLY: wp,i4

  IMPLICIT NONE

  PUBLIC :: &
  ! emiss
       &    ntime_emiss,minimal_emiss, undef_emiss, &
  ! ndvi
       &    undef_ndvi, minimal_ndvi, ntime_ndvi, &                    
  !albedo
       &    ntime_alb, &
       &    ialb_type, &
       &    undef_alb, minimal_alb, &
       &    undef_alb_bs, minimal_alb_dry, maximal_alb_dry, minimal_alb_sat, &
       &    maximal_alb_sat, &
       &    zalso, wso_min, wso_max, csalb, csalbw, &
       &    allocate_alb_interp_fields,   &
       &    alb_interp_data



  INTEGER (KIND=i4)            :: &
  ! emiss
       &                          ntime_emiss = 12, & !< number of timesteps (12 for monthly mean values)
  ! ndvi
       &                          ntime_ndvi = 12, &  !< number of timesteps (12 for monthly mean values)
  ! albedo
       &                          ntime_alb = 12, &   !< number of timesteps (12 for monthly mean values)
       &                          ialb_type = 1, &    !< VIS,NIR,UV (1), soil (2), VIS (3)
  ! cru
       &                          i_t_cru_fine = 1, &
       &                          i_t_cru_coarse = 2

  REAL (KIND=wp)               :: &
  ! emiss
       &                          undef_emiss = 0.0, &   !< undefined value for EMISS data
       &                          minimal_emiss = 0.0, & !< minimal EMISS value 
  ! ndvi
       &                          undef_ndvi = 0.0, &    !< undefined value for NDVI data
       &                          minimal_ndvi = 0.09, & !< minimal NDVI value bare soil value
  ! albedo
       &                          undef_alb = 0.0, & !< undefined value for ALB data
       &                          minimal_alb = 0.07, & !< minimal ALB value open sea
       &                          undef_alb_bs = -1.E20, &  !< undefined value for bare soil albedo
       &                          minimal_alb_dry = 0.1224, & !< minimal value for dry soil
       &                          maximal_alb_dry = 0.4925, & !< maximum value for dry soil
       &                          minimal_alb_sat = 0.0612, & !< minimal value for saturated soil
       &                          maximal_alb_sat = 0.3825 !< maximum value for saturated soil

  REAL (KIND=wp), ALLOCATABLE  :: & 
  ! albedo
    &                             zalso(:,:), &
    &                             wso_min(:), wso_max(:), &
    &                             csalb(:), csalbw(:)

  CONTAINS

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

END MODULE mo_python_data
