!  Fortran module containing the subroutines to determine the SSO parameters
!
!  History:
!  Version                  Date                 Name
!  ------------------------ -------------------- -----------------
!  V2_0                     2013/04/17           Martina Messmer
!    initial release
!    calculation of SSO parameters are separated from the main program
!    the SSO parameters are only calculated if desired due to performance
!    reasons
! 
!  Code Description:
!  Language: Fortran 90
!=================================================================
MODULE mo_topo_sso

! Modules used: 

  !> kind parameters are defined in MODULE data_parameters
  USE mo_kind, ONLY: wp
  USE mo_kind, ONLY: i4
  USE mo_kind, ONLY: i8


IMPLICIT NONE

PRIVATE

PUBLIC :: auxiliary_sso_parameter, &
          calculate_sso

  CONTAINS

   SUBROUTINE auxiliary_sso_parameter_icon(d2x,d2y,j_n,j_c,j_s,hh,nc,dxrat,dhdx,dhdy,dhdxdx,dhdydy,dhdxdy)

     USE mo_topo_data, ONLY: nc_tot !< number of total GLOBE/ASTER columns un a latitude circle

     REAL(KIND=wp),   INTENT(IN) ::  d2x,dxrat ! 2 times grid distance for gradient calculation (in [m])
     REAL(KIND=wp),   INTENT(IN) ::  d2y       ! 2 times grid distance for gradient calculation (in [m])
     INTEGER(KIND=i4),INTENT(IN) :: j_n
     INTEGER(KIND=i4),INTENT(IN) :: j_c 
     INTEGER(KIND=i4),INTENT(IN) :: j_s, nc

!DR     INTEGER(KIND=i4),INTENT(INOUT) :: hh(0:nc_tot+1,1:3) !< topographic height for gradient calculations
     REAL(KIND=wp),   INTENT(INOUT) :: hh(0:nc_tot+1,1:3) !< topographic height for gradient calculations
     REAL(KIND=wp),   INTENT(OUT):: dhdxdx(1:nc_tot)  !< x-gradient square for one latitude row
     REAL(KIND=wp),   INTENT(OUT):: dhdydy(1:nc_tot)  !< y-gradient square for one latitude row
     REAL(KIND=wp),   INTENT(OUT):: dhdxdy(1:nc_tot)  !< dxdy for one latitude row

     REAL(KIND=wp),   INTENT(OUT)  :: dhdx(1:nc_tot)    !< x-gradient for one latitude row
     REAL(KIND=wp),   INTENT(OUT)  :: dhdy(1:nc_tot)    !< y-gradient for one latitude row

      
     INTEGER(KIND=i4):: i   !< counter

 
     DO i = 1,nc
       dhdx(i) = (hh(i+1,j_c) - hh(i-1,j_c))/(d2x*dxrat)  ! centered differences as gradient, except for mlat=1 and mlat= 21600
       dhdy(i) = (hh(i,j_n) - hh(i,j_s))/ABS(d2y)
!!$       dhdxdx(i) = ((hh(i+1,j_c) - hh(i-1,j_c))/d2x)**2.
!!$       dhdydy(i) = ((hh(i,j_s) - hh(i,j_n))/d2y)**2.
!!$       dhdxdy(i) = ((hh(i+1,j_c) - hh(i-1,j_c))/d2x)*((hh(i,j_s) - hh(i,j_n))/d2y)
     ENDDO

     dhdxdx(1:nc) = dhdx(1:nc) * dhdx(1:nc) ! x-gradient square
     dhdydy(1:nc) = dhdy(1:nc) * dhdy(1:nc) ! y-gradient square
     dhdxdy(1:nc) = dhdx(1:nc) * dhdy(1:nc) ! dx*dy

   END SUBROUTINE auxiliary_sso_parameter_icon

   SUBROUTINE auxiliary_sso_parameter_cosmo(d2x,d2y,j_n,j_c,j_s,hh,dhdxdx,dhdydy,dhdxdy)

     USE mo_topo_data, ONLY: nc_tot !< number of total GLOBE/ASTER columns un a latitude circle

     REAL(KIND=wp),   INTENT(IN) ::  d2x       ! 2 times grid distance for gradient calculation (in [m])
     REAL(KIND=wp),   INTENT(IN) ::  d2y       ! 2 times grid distance for gradient calculation (in [m])
     INTEGER(KIND=i4),INTENT(IN) :: j_n
     INTEGER(KIND=i4),INTENT(IN) :: j_c
     INTEGER(KIND=i4),INTENT(IN) :: j_s

     INTEGER(KIND=i4),INTENT(INOUT) :: hh(0:nc_tot+1,1:3) !< topographic height for gradient calculations
     REAL(KIND=wp),   INTENT(OUT):: dhdxdx(1:nc_tot)  !< x-gradient square for one latitude row
     REAL(KIND=wp),   INTENT(OUT):: dhdydy(1:nc_tot)  !< y-gradient square for one latitude row
     REAL(KIND=wp),   INTENT(OUT):: dhdxdy(1:nc_tot)  !< dxdy for one latitude row

     REAL(KIND=wp)   :: dhdx(1:nc_tot)    !< x-gradient for one latitude row
     REAL(KIND=wp)   :: dhdy(1:nc_tot)    !< y-gradient for one latitude row
     
     INTEGER(KIND=i4):: i   !< counter


     DO i = 1,nc_tot
       dhdx(i) = (hh(i+1,j_c) - hh(i-1,j_c))/d2x  ! centered differences as gradient, except for mlat=1 and mlat= 21600
       dhdy(i) = (hh(i,j_n) - hh(i,j_s))/d2y 
     ENDDO

     dhdxdx(1:nc_tot) = dhdx(1:nc_tot) * dhdx(1:nc_tot) ! x-gradient square
     dhdydy(1:nc_tot) = dhdy(1:nc_tot) * dhdy(1:nc_tot) ! y-gradient square
     dhdxdy(1:nc_tot) = dhdx(1:nc_tot) * dhdy(1:nc_tot) ! dx*dy


   END SUBROUTINE auxiliary_sso_parameter_cosmo

!--------------------------------------------------------------------------------------------------------------------------

   SUBROUTINE calculate_sso(tg,no_raw_data_pixel,   &
        &                   h11,h12,h22,stdh_target,&
        &                   theta_target,           &
        &                   aniso_target,           &
        &                   slope_target)

     USE mo_cosmo_grid,       ONLY: COSMO_grid !< structure which contains the definition of the COSMO grid
     USE mo_target_grid_data, ONLY: lon_geo, & !< longitude coordinates of the grid in the geographical system 
     &                              lat_geo !< latitude coordinates of the grid in the geographical system
     USE mo_utilities_extpar, ONLY: uv2uvrot
     USE mo_math_constants,   ONLY: pi
     USE mo_grid_structures,  ONLY: target_grid_def  !< Definition of data type with target grid definition
     USE mo_grid_structures,  ONLY: igrid_icon
     USE mo_grid_structures,  ONLY: igrid_cosmo
     USE mo_grid_structures,  ONLY: igrid_gme

     TYPE(target_grid_def), INTENT(IN) :: tg              !< structure with target grid description
     INTEGER(KIND=i8),      INTENT(IN) :: no_raw_data_pixel(:,:,:)
     REAL(KIND=wp),         INTENT(IN) :: h11(:,:,:)      !< help variables
     REAL(KIND=wp),         INTENT(IN) :: h12(:,:,:)      !< help variables
     REAL(KIND=wp),         INTENT(IN) :: h22(:,:,:)      !< help variables
     REAL(KIND=wp),         INTENT(IN) :: stdh_target(:,:,:)
     REAL(KIND=wp),         INTENT(OUT):: theta_target(1:tg%ie,1:tg%je,1:tg%ke) !< sso parameter, angle of principal axis
     REAL(KIND=wp),         INTENT(OUT):: aniso_target(1:tg%ie,1:tg%je,1:tg%ke) !< sso parameter, anisotropie factor
     REAL(KIND=wp),         INTENT(OUT):: slope_target(1:tg%ie,1:tg%je,1:tg%ke) !< sso parameter, mean slope
!< standard deviation of subgrid scale orographic height
     REAL(KIND=wp) :: point_lon, point_lat
     REAL(KIND=wp) :: znorm
     REAL(KIND=wp) :: zh11, zh12, zh22
     REAL(KIND=wp) :: K_lm, L_lm, M_lm                     !< Variables to determine angle of principal axis, 
                                                           !< anisotropy and slope after Lott and Miller 96
     REAL(KIND=wp) :: K_lm_prime, L_lm_prime, M_lm_prime   !< Variables to determine angle of principal axis, 
                                                           !< anisotropy and slope after Lott and Miller 96
     REAL(KIND=wp) :: theta                 ! angle of principle axis
     REAL(KIND=wp) :: theta_rot             ! angle of principle axis in the rotated system
     REAL(KIND=wp) :: theta_u, theta_v      ! help variables for the rotation of theta into the rotated system
     REAL(KIND=wp) :: theta_urot, theta_vrot! help variables for the rotation of theta into the rotated system
     REAL(KIND=wp) :: gamma_lm              ! ansisotropy factor
     REAL(KIND=wp) :: gamma_lm2             ! anisotropy factor square
     REAL(KIND=wp) :: zaehler               ! help variable
     REAL(KIND=wp) :: nenner                ! help variable
     REAL(KIND=wp) :: sigma                 ! slope parameter

     INTEGER(KIND=i4):: ke, je, ie !< counters

       print *,'SSO parameters'
      ! SSO parameters
      ! angle of principal axis
       DO ke=1, tg%ke
       DO je=1, tg%je
       DO ie=1, tg%ie
         IF (no_raw_data_pixel(ie,je,ke) > 1) THEN
           znorm = 1.0/no_raw_data_pixel(ie,je,ke) 
         ELSE
           znorm = 0.
         ENDIF
         IF (stdh_target(ie,je,ke) > 10.0) THEN ! avoid trivial case of sea point
           ! SSO parameters
           zh11 = h11(ie,je,ke) * znorm 
           zh12 = h12(ie,je,ke) * znorm 
           zh22 = h22(ie,je,ke) * znorm 
          
          ! calculation of angle of principal axis
          !----------------------------------------------------------------------------------
           ! Equation (A.1) of Lott and Miller, 1996
           K_lm = 0.5 * (zh11 + zh22)
           L_lm = 0.5 * (zh11 - zh22) 
           M_lm = zh12

          ! angle of principle axis
           theta = 0.5 * ATAN2(M_lm,L_lm)    ! Lott and Miller 1996, equation (A.2)
           SELECT CASE(tg%igrid_type)
           CASE(igrid_icon)  ! ICON GRID
             theta_target(ie,je,ke) = theta  
           CASE(igrid_cosmo)  ! COSMO GRID
              ! compute theta in the rotated grid
             theta_u = COS(theta)
             theta_v = SIN(theta)
             point_lon = lon_geo(ie,je,ke)
             point_lat = lat_geo(ie,je,ke)
             CALL uv2uvrot(theta_u, theta_v,           &
                           point_lat, point_lon,       &
                           COSMO_grid%pollat, COSMO_grid%pollon, &
                           theta_urot, theta_vrot)
             theta_rot = ATAN2(theta_vrot,theta_urot)            ! angle of principle axis in the rotated system
              ! Restrict the range of values of theta to [-pi/2,pi/2]
             IF (theta_rot < -pi/2.) theta_rot = theta_rot + pi
             IF (theta_rot >  pi/2.) theta_rot = theta_rot - pi
             theta_target(ie,je,ke) = theta_rot  ! angle of principle axis in the rotated system
           CASE(igrid_gme) ! GME GRID
             theta_target(ie,je,ke) = theta
           END SELECT
          !----------------------------------------------------------------------------------
          ! calculation of anisotropy factor
          !----------------------------------------------------------------------------------
           K_lm_prime = K_lm
           L_lm_prime = SQRT(L_lm*L_lm + M_lm*M_lm)
           M_lm_prime = 0.0
           zaehler = K_lm_prime - L_lm_prime
           nenner  = K_lm_prime + L_lm_prime
           IF (zaehler <= EPSILON(zaehler) ) zaehler = 0.0
           IF (nenner  <= EPSILON(nenner) )  nenner  = EPSILON(nenner)
           gamma_lm2 = zaehler / nenner
           gamma_lm = SQRT(gamma_lm2) ! Lott and Miller 1996, equation (A.3)
           aniso_target(ie,je,ke) = gamma_lm  ! anisotropy factor
          !----------------------------------------------------------------------------------
          ! Calculation of slope parameter
          !----------------------------------------------------------------------------------
           sigma = K_lm_prime + L_lm_prime
           sigma = SQRT(sigma) ! Lott and Miller 1996, equation (A.5)
           slope_target(ie,je,ke) = sigma
         ELSE  ! seapoints or other points with STDH <= 10 m
           theta_target(ie,je,ke) = 0.0
           aniso_target(ie,je,ke) = 0.0
           slope_target(ie,je,ke) = 0.0
         ENDIF
       ENDDO
       ENDDO
       ENDDO

     END SUBROUTINE calculate_sso

END MODULE mo_topo_sso
