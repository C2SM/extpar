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
    &   1.e-10,          &
    &   1.e-10,             &
    &   0.43,            &
    &   0.41,            &
    &   0.43,            &
    &   0.41,            &
    &   0.507,            &
    &   0.73,            &
    &   1.e-10,           &
    &   1.e-10 /)

  REAL (KIND=wp) :: cfcap_vg(ntype_fao)  = (/ &       !< lookup table field capacita (fraction of volume)
    &   1.e-10,          &
    &   1.e-10,          &
    &   0.196,           &
    &   0.260,           &
    &   0.340,           &
    &   0.370,           &
    &   0.463,           &
    &   0.70,            &
    &   1.e-10,          &
    &   1.e-10           /)

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
    &    1.e-10,         &
    &    1.e-10,         &
    &    0.145,          &
    &    0.075,          &
    &    0.036,          &
    &    0.019,          &
    &    0.008,          &
    &    0.01,           &
    &    1.e-10,         &
    &    1.e-10          /)

   REAL (KIND=wp) :: crocg_vg(ntype_fao)  = (/ & ! * 1.E6
    &    1.92 ,         &
    &    2.10 ,         &
    &    1.28 ,         &
    &    1.35 ,         &
    &    1.42 ,         &
    &    1.50 ,         &
    &    1.63 ,         &
    &    0.58 ,         &
    &    4.18 ,         &
    &    1.92          /) 
  
   REAL (KIND=wp) :: cala0_vg(ntype_fao)  = (/ &
    &   2.26  ,         &
    &   2.41  ,         &
    &   0.30  ,         &
    &   0.28  ,         &
    &   0.25  ,         &
    &   0.21  ,         &
    &   0.18  ,         &
    &   0.06  ,         &
    &   1.0   ,         &
    &   2.26           /)

   REAL (KIND=wp) :: cala1_vg(ntype_fao)  = (/ &
    &   2.26  ,         &
    &   2.41  ,         &
    &   2.40  ,         &
    &   2.40  ,         &
    &   1.58  ,         &
    &   1.55  ,         &
    &   1.50  ,         &
    &   0.50  ,         &
    &   1.0   ,         &
    &   2.26           /)   

   REAL (KIND=wp) :: csand_vg(ntype_fao)  = (/ &
    &   0.0  ,         &
    &   0.0  ,         &
    &   90.0  ,         &
    &   65.0  ,         &
    &   40.0  ,         &
    &   35.0  ,         &
    &   15.0  ,         &
    &   90.0  ,         &
    &   0.0   ,         &
    &   0.0           /) 

   REAL (KIND=wp) :: cclay_vg(ntype_fao)  = (/ &
    &   0.0  ,         &
    &   0.0  ,         &
    &   5.0  ,         &
    &   10.0  ,         &
    &   20.0  ,         &
    &   35.0  ,         &
    &   70.0  ,         &
    &   5.0  ,         &
    &   0.0   ,         &
    &   0.0           /) 

   REAL (KIND=wp) :: csilt_vg(ntype_fao)  = (/ &
    &   0.0  ,         &
    &   0.0  ,         &
    &   5.0  ,         &
    &   25.0  ,         &
    &   40.0  ,         &
    &   30.0  ,         &
    &   15.0  ,         &
    &   5.0  ,         &
    &   0.0   ,         &
    &   0.0           /) 

   
END MODULE mo_vg_lookup_tables
