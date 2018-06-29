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

  USE mo_kind,             ONLY: wp, i4, i8
  USE mo_topo_data,        ONLY: nc_tot !< number of total GLOBE/ASTER columns un a latitude circle
  USE mo_cosmo_grid,       ONLY: cosmo_grid !< structure which contains the definition of the COSMO grid
  USE mo_target_grid_data, ONLY: lon_geo, lat_geo !< longitude, latitude of the grid in the geographical system
  USE mo_utilities_extpar, ONLY: uv2uvrot
  USE mo_math_constants,   ONLY: pi
  USE mo_grid_structures,  ONLY: target_grid_def, &  !< Definition of data type with target grid definition
       &                         igrid_icon, igrid_cosmo

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: auxiliary_sso_parameter_icon, auxiliary_sso_parameter_cosmo, calculate_sso

CONTAINS

  SUBROUTINE auxiliary_sso_parameter_icon(d2x,d2y,j_n,j_c,j_s,hh,nc,dxrat,dhdx,dhdy,dhdxdx,dhdydy,dhdxdy)

    REAL(wp),    INTENT(in) :: d2x       ! 2 times grid distance for gradient calculation (in [m])
    REAL(wp),    INTENT(in) :: d2y       ! 2 times grid distance for gradient calculation (in [m])
    INTEGER(i4), INTENT(in) :: j_n
    INTEGER(i4), INTENT(in) :: j_c
    INTEGER(i4), INTENT(in) :: j_s

    REAL(wp),    INTENT(inout) :: hh(0:nc_tot+1,1:3) !< topographic height for gradient calculations

    INTEGER(i4), INTENT(in)    :: nc
    REAL(wp),    INTENT(in)    :: dxrat
    REAL(wp),    INTENT(out)   :: dhdxdx(1:nc_tot)  !< x-gradient square for one latitude row
    REAL(wp),    INTENT(out)   :: dhdydy(1:nc_tot)  !< y-gradient square for one latitude row
    REAL(wp),    INTENT(out)   :: dhdxdy(1:nc_tot)  !< dxdy for one latitude row

    REAL(wp),    INTENT(out)   :: dhdx(1:nc_tot)    !< x-gradient for one latitude row
    REAL(wp),    INTENT(out)   :: dhdy(1:nc_tot)    !< y-gradient for one latitude row

    INTEGER(i4):: i

    DO i = 1, nc
      dhdx(i) = (hh(i+1,j_c) - hh(i-1,j_c))/(d2x*dxrat)  ! centered differences as gradient, except for mlat=1 and mlat= 21600
      dhdy(i) = (hh(i,j_n) - hh(i,j_s))/ABS(d2y)
    ENDDO

    dhdxdx(1:nc) = dhdx(1:nc) * dhdx(1:nc) ! x-gradient square
    dhdydy(1:nc) = dhdy(1:nc) * dhdy(1:nc) ! y-gradient square
    dhdxdy(1:nc) = dhdx(1:nc) * dhdy(1:nc) ! dx*dy

  END SUBROUTINE auxiliary_sso_parameter_icon

  SUBROUTINE auxiliary_sso_parameter_cosmo(d2x,d2y,j_n,j_c,j_s,hh,dhdxdx,dhdydy,dhdxdy)

    REAL(wp),   INTENT(in) ::  d2x       ! 2 times grid distance for gradient calculation (in [m])
    REAL(wp),   INTENT(in) ::  d2y       ! 2 times grid distance for gradient calculation (in [m])
    INTEGER(i4),INTENT(in) :: j_n
    INTEGER(i4),INTENT(in) :: j_c
    INTEGER(i4),INTENT(in) :: j_s

    INTEGER(i4),INTENT(inout) :: hh(0:nc_tot+1,1:3) !< topographic height for gradient calculations
    REAL(wp),   INTENT(out):: dhdxdx(1:nc_tot)  !< x-gradient square for one latitude row
    REAL(wp),   INTENT(out):: dhdydy(1:nc_tot)  !< y-gradient square for one latitude row
    REAL(wp),   INTENT(out):: dhdxdy(1:nc_tot)  !< dxdy for one latitude row

    REAL(wp)   :: dhdx(1:nc_tot)    !< x-gradient for one latitude row
    REAL(wp)   :: dhdy(1:nc_tot)    !< y-gradient for one latitude row

    INTEGER(i4):: i   !< counter

    DO i = 1, nc_tot
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

    TYPE(target_grid_def), INTENT(in) :: tg              !< structure with target grid description
    INTEGER,          INTENT(in) :: no_raw_data_pixel(:,:,:)
    REAL(wp),         INTENT(in) :: h11(:,:,:)      !< help variables
    REAL(wp),         INTENT(in) :: h12(:,:,:)      !< help variables
    REAL(wp),         INTENT(in) :: h22(:,:,:)      !< help variables
    REAL(wp),         INTENT(in) :: stdh_target(:,:,:)
    REAL(wp),         INTENT(out):: theta_target(1:tg%ie,1:tg%je,1:tg%ke) !< sso parameter, angle of principal axis
    REAL(wp),         INTENT(out):: aniso_target(1:tg%ie,1:tg%je,1:tg%ke) !< sso parameter, anisotropie factor
    REAL(wp),         INTENT(out):: slope_target(1:tg%ie,1:tg%je,1:tg%ke) !< sso parameter, mean slope
    !< standard deviation of subgrid scale orographic height
    REAL(wp) :: point_lon, point_lat
    REAL(wp) :: znorm
    REAL(wp) :: zh11, zh12, zh22
    REAL(wp) :: K_lm, L_lm, M_lm                     !< Variables to determine angle of principal axis,
    !< anisotropy and slope after Lott and Miller 96
    REAL(wp) :: K_lm_prime, L_lm_prime               !< Variables to determine angle of principal axis,
    !< anisotropy and slope after Lott and Miller 96
    REAL(wp) :: theta                 ! angle of principle axis
    REAL(wp) :: theta_rot             ! angle of principle axis in the rotated system
    REAL(wp) :: theta_u, theta_v      ! help variables for the rotation of theta into the rotated system
    REAL(wp) :: theta_urot, theta_vrot! help variables for the rotation of theta into the rotated system
    REAL(wp) :: gamma_lm              ! ansisotropy factor
    REAL(wp) :: gamma_lm2             ! anisotropy factor square
    REAL(wp) :: zaehler               ! help variable
    REAL(wp) :: nenner                ! help variable
    REAL(wp) :: sigma                 ! slope parameter

    INTEGER(i4):: ke, je, ie !< counters

    print *,'SSO parameters'
    ! SSO parameters
    ! angle of principal axis
    DO ke = 1, tg%ke
      DO je = 1, tg%je
        DO ie = 1, tg%ie
          IF (stdh_target(ie,je,ke) > 10.0_wp) THEN ! avoid trivial case of sea point
            IF (no_raw_data_pixel(ie,je,ke) > 1) THEN
              znorm = 1.0_wp/no_raw_data_pixel(ie,je,ke)
            ELSE
              znorm = 0.0_wp
            ENDIF
            ! SSO parameters
            if (ie > 13873626) then
              write(*,*)ie,je,ke
              write(*,*)no_raw_data_pixel(ie,je,ke)
              write(*,*)znorm
              write(*,*)h11(ie,je,ke)
              write(*,*)h12(ie,je,ke)
              write(*,*)h22(ie,je,ke)
            endif
            zh11 = h11(ie,je,ke) * znorm
            zh12 = h12(ie,je,ke) * znorm
            zh22 = h22(ie,je,ke) * znorm
            if (ie > 13873626) then
              write(*,*)zh11,zh12,zh22
            endif
            ! calculation of angle of principal axis
            !----------------------------------------------------------------------------------
            ! Equation (A.1) of Lott and Miller, 1996
            K_lm = 0.5_wp * (zh11 + zh22)
            L_lm = 0.5_wp * (zh11 - zh22)
            M_lm = zh12

            ! angle of principle axis
            theta = 0.5_wp * ATAN2(M_lm,L_lm)    ! Lott and Miller 1996, equation (A.2)
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
            END SELECT
            !----------------------------------------------------------------------------------
            ! calculation of anisotropy factor
            !----------------------------------------------------------------------------------
            K_lm_prime = K_lm
            L_lm_prime = SQRT(L_lm*L_lm + M_lm*M_lm)
            zaehler = K_lm_prime - L_lm_prime
            nenner  = K_lm_prime + L_lm_prime
            IF (zaehler <= EPSILON(zaehler) ) zaehler = 0.0_wp
            IF (nenner  <= EPSILON(nenner) )  nenner  = EPSILON(nenner)
            gamma_lm2 = zaehler / nenner
            gamma_lm = SQRT(gamma_lm2)         ! Lott and Miller 1996, equation (A.3)
            aniso_target(ie,je,ke) = gamma_lm  ! anisotropy factor
            !----------------------------------------------------------------------------------
            ! Calculation of slope parameter
            !----------------------------------------------------------------------------------
            sigma = K_lm_prime + L_lm_prime
            sigma = SQRT(sigma) ! Lott and Miller 1996, equation (A.5)
            slope_target(ie,je,ke) = sigma
          ELSE  ! seapoints or other points with STDH <= 10 m
            theta_target(ie,je,ke) = 0.0_wp
            aniso_target(ie,je,ke) = 0.0_wp
            slope_target(ie,je,ke) = 0.0_wp
          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE calculate_sso

END MODULE mo_topo_sso
