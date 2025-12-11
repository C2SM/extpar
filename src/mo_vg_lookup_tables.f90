!+ Fortran Module with lookup-tables for vanGenuchten parameters for different soil types
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2025/11/12 Linda Schlemmer
!  Initial release
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran Module with lookup-tables for van Genuchten parameters for different soil types
!> \author Linda Schlemmer
!!
!! Description:
!! van Genuchten parameters for FAO soiltypes as backup for data points, at which values provided by the HiHydroSoil
!! dataset holds unplausible values
!!

MODULE mo_vg_lookup_tables

  USE mo_kind,                  ONLY: wp, i4

  IMPLICIT NONE


  PUBLIC :: ntype_fao

  ! small value
  REAL, PARAMETER :: eps = EPSILON(1._wp)
  INTEGER (KIND=i4), PARAMETER :: ntype_fao = 10  !< FAO has 10 soil types

  !---------------------------------------------------------------------------------------------- 
  !----------------------------------------------------------------------------------------------
  REAL (KIND=wp) :: cporv_vg(ntype_fao)  = (/ &       !< lookup table pore volume (fraction of volume)
    &   1.0, &!eps,         &
    &   1.0, &!eps,         &
    &   0.43,            &
    &   0.41,            &
    &   0.43,            &
    &   0.41,            &
    &   0.507,           &
    &   0.73 ,            &
    &   1.0, &!eps,         &
    &   1.0 /)!eps /)

  REAL (KIND=wp) :: cfcap_vg(ntype_fao)  = (/ &       !< lookup table field capacita (fraction of volume)
    &   eps,          &
    &   eps,          &
    &   0.196,           &
    &   0.260,           &
    &   0.340,           &
    &   0.370,           &
    &   0.463,           &
    &   0.70,            &
    &   eps,          &
    &   eps           /)

  REAL (KIND=wp) :: cpwp_vg(ntype_fao)   = (/ &
    &   0.0,             &
    &   0.0,             &
    &   0.046,           &
    &   0.100,           &
    &   0.110,           &
    &   0.185,           &
    &   0.257,           &
    &   0.265,           &
    &   0.0,             &
    &   0.0              /)


  REAL (KIND=wp) :: cadp_vg(ntype_fao)   = (/ &
    &   0.0,             &
    &   0.0,             &
    &   0.045,           &
    &   0.065,           &
    &   0.078,           &
    &   0.095,           &
    &   0.068,           &
    &   0.098,           &
    &   0.0,             &
    &   0.0              /)

  REAL (KIND=wp) :: ckw0_vg(ntype_fao)   = (/ &
    &   0.0,             &
    &   0.0,             &
    & 712.8,             &
    & 106.1,             &
    &  25.0,             &
    &   6.2,             &
    &   4.8,             &
    &  13.44,            &
    &    0.0,            &
    &    0.0             /)

  REAL (KIND=wp) :: n_vg(ntype_fao)      = (/ &
    &    1.E-6,          &
    &    1.E-6,          &
    &    2.3,            &
    &    1.89,           &
    &    1.56,           &
    &    1.31,           &
    &    1.09,           &
    &    1.32,           &
    &    1.E-6,          &
    &    1.E-6           /)

  REAL (KIND=wp) :: alpha_vg(ntype_fao)  = (/ &
    &    eps,         &
    &    eps,         &
    &    0.145,          &
    &    0.075,          &
    &    0.036,          &
    &    0.019,          &
    &    0.008,          &
    &    0.01,           &
    &    eps,         &
    &    eps          /)


END MODULE mo_vg_lookup_tables
