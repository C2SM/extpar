!+ Fortran module with variables and routines for TERRA_URB
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2023/02/17 Jacopo Canton
!  Initial release
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module with variables and routines for TERRA_URB
!> \author Jacopo Canton
!!
!! Description:
!! Variables and routines used to generate the additional data required by
!! TERRA_URB.
!! The code is heavily based upon Matthias Demuzere's WUDAPT-TO-COSMO scripts
!! (https://github.com/matthiasdemuzere/WUDAPT-to-COSMO) which, in turn, are
!! based on SURY (https://github.com/hendrikwout/sury).
!!
MODULE mo_terra_urb

  USE mo_kind,   ONLY: wp, i4

  USE mo_logging

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: tu_test_all

  INTEGER(KIND=i4), PARAMETER :: nr_lcz = 10, & !< number of LCZ classes
                                 nr_ucp_values = 20 !< number of values defined in the tables

  TYPE :: lcz_class
    INTEGER(KIND=i4) :: class_nr
    ! values defined in the table
    REAL(KIND=wp)    :: ISA,        &
                        AHF,        &
                        FR_PAVED,   &
                        URBAN,      &
                        URB_BLDFR,  &
                        URB_BLDH,   &
                        URB_H2W,    &
                        URB_RfALB,  &
                        URB_WaALB,  &
                        URB_RdALB,  &
                        URB_RfEMI,  &
                        URB_WaEMI,  &
                        URB_RdEMI,  &
                        URB_RfHCAP, &
                        URB_WaHCAP, &
                        URB_RdHCAP, &
                        URB_RfHCON, &
                        URB_WaHCON, &
                        URB_RdHCON
    ! derived bulk values
    REAL(KIND=wp)    :: URB_SALB_BK, &
                        URB_EMIS_BK, &
                        URB_TALB_BK, &
                        URB_HCON,    &
                        URB_HCAP,    &
                        URB_EMIS_FL, &
                        URB_SALB_FL, &
                        URB_TALB_FL, &
                        URB_SALB,    &
                        URB_TALB,    &
                        URB_EMIS
  END TYPE lcz_class

  TYPE(lcz_class) :: ucp(nr_lcz) !< key-based version of the table below

  REAL(KIND=wp), DIMENSION(nr_lcz, nr_ucp_values), PARAMETER :: lcz_ucp_default = reshape( (/ & !< LCZ_UCP_default.csv from Matthias' WUDAPT-TO-COSMO
  ! class_nr, ISA,  FR_PAVED, URBAN, URB_BLDFR, URB_BLDH, URB_H2W, URB_RfALB, URB_WaALB, URB_RdALB, URB_RfEMI, URB_WaEMI, URB_RdEMI, URB_RfHCAP, URB_WaHCAP, URB_RdHCAP, URB_RfHCON, URB_WaHCON, URB_RdHCON, AHF
     1.,      0.95, 0.95,     0.95,  0.5,       25.0,     2.5,     0.13,      0.25,      0.14,      0.91,      0.9,       0.95,      1800000.,   1800000.,   1750000.,   1.25,       1.09,       0.77,       100., &
     2.,      0.9,  0.9,      0.9,   0.5,       15.0,     1.25,    0.18,      0.2,       0.14,      0.91,      0.9,       0.95,      1800000.,   2670000.,   1680000.,   1.25,       1.50,       0.73,        35., &
     3.,      0.85, 0.85,     0.85,  0.55,       5.0,     1.25,    0.15,      0.2,       0.14,      0.91,      0.9,       0.95,      1440000.,   2050000.,   1630000.,   1.00,       1.25,       0.69,        30., &
     4.,      0.65, 0.65,     0.65,  0.3,       25.0,     1.0,     0.13,      0.25,      0.14,      0.91,      0.9,       0.95,      1800000.,   2000000.,   1540000.,   1.25,       1.45,       0.64,        30., &
     5.,      0.7,  0.7,      0.7,   0.3,       15.0,     0.5,     0.13,      0.25,      0.14,      0.91,      0.9,       0.95,      1800000.,   2000000.,   1500000.,   1.25,       1.45,       0.62,        15., &
     6.,      0.6,  0.6,      0.6,   0.3,        5.0,     0.5,     0.13,      0.25,      0.14,      0.91,      0.9,       0.95,      1440000.,   2050000.,   1470000.,   1.00,       1.25,       0.60,        10., &
     7.,      0.85, 0.85,     0.85,  0.8,        3.0,     1.5,     0.15,      0.2,       0.18,      0.28,      0.9,       0.92,      2000000.,    720000.,   1670000.,   2.00,       0.50,       0.72,        30., &
     8.,      0.85, 0.85,     0.85,  0.4,        7.0,     0.2,     0.18,      0.25,      0.14,      0.91,      0.9,       0.95,      1800000.,   1800000.,   1380000.,   1.25,       1.25,       0.51,        40., &
     9.,      0.3,  0.3,      0.3,   0.15,       5.0,     0.15,    0.13,      0.25,      0.14,      0.91,      0.9,       0.95,      1440000.,   2560000.,   1370000.,   1.00,       1.00,       0.55,         5., &
    10.,      0.55, 0.55,     0.55,  0.25,       8.5,     0.35,    0.1,       0.2,       0.14,      0.91,      0.9,       0.95,      2000000.,   1690000.,   1490000.,   2.00,       1.33,       0.61,       300.  &
    /), (/nr_lcz, nr_ucp_values/), order=(/2, 1/))

  CONTAINS

    !---------------------------------------------------------------------------
    !> Main routine to prepare fields
    SUBROUTINE terra_urb_to_extpar

      INTEGER :: i

      CALL logging%info('Enter routine: terra_urb_to_extpar')

      ! Prepare the urban canopy data
      CALL prepare_ucp_lookup()

      CALL logging%info('Exit routine: terra_urb_to_extpar')

    END SUBROUTINE terra_urb_to_extpar

    !---------------------------------------------------------------------------
    !> Helper function to prepare the urban canopy data
    !>
    !> The code generally follows SURY: https://github.com/hendrikwout/sury
    !>
    !> * Wouters, H., Demuzere, M., Blahak, U., Fortuniak, K., Maiheu., B.,
    !>     Camps, J., Tielemans, and N. P. M. van Lipzig, 2016.  The efficient
    !>     urban-canopy dependency parametrization SURY (v1.0) for atmospheric
    !>     modelling: description and application with the COSMO-CLM model
    !>     (v5.0_clm6) for a Belgian Summer, Geosci. Model Dev., 2016.
    !>
    !> Define the look-up table, based on the values of:
    !> * Stewart, I. D., & Oke, T. R. (2012). Local Climate Zones for Urban
    !>     Temperature Studies. Bulletin of the American Meteorological
    !>     Society, 93(12), 1879–1900.
    !> * Stewart, I. D., Oke, T. R., & Krayenhoff, E. S. (2014). Evaluation of
    !>     the ‘local climate zone’ scheme using temperature observations and
    !>     model simulations. International Journal of Climatology, 34(4),
    !>     1062–1080. https://doi.org/10.1002/joc.3746
    !>
    !> The latter paper describes thermal admittance values of facets only.
    !> Heat conductivity and capacity values are obtained via Scott Krayenhoff
    !> (personal communication).
    SUBROUTINE prepare_ucp_lookup

      INTEGER :: i
      LOGICAL,       PARAMETER :: saiWeight = .FALSE. !< Weigh parameters according to Surface Area Index (Default = False)
      REAL(KIND=wp), PARAMETER :: snow_f = 0.0,   &   !< snow fraction (default = 0)
                                  alb_snow = 0.7, &   !< snow albedo (default = 0.7)
                                  emi_snow = 0.997    !< emissivity albedo (default = 0.997)

      REAL(KIND=wp) :: psi_canyon, psi_bulk, &
                       alb_roof_snow, alb_road_snow, alb_wall_snow, &
                       emi_roof_snow, emi_road_snow, emi_wall_snow, &
                       SAI

      CALL logging%info('Enter routine: prepare_ucp_lookup')

      ! Fill the UCP table for each LCZ class
      DO i = 1,nr_lcz
        ucp(i)%class_nr   = int(lcz_ucp_default(i, 1))
        ucp(i)%ISA        =     lcz_ucp_default(i, 2)
        ucp(i)%FR_PAVED   =     lcz_ucp_default(i, 3)
        ucp(i)%URBAN      =     lcz_ucp_default(i, 4)
        ucp(i)%URB_BLDFR  =     lcz_ucp_default(i, 5)
        ucp(i)%URB_BLDH   =     lcz_ucp_default(i, 6)
        ucp(i)%URB_H2W    =     lcz_ucp_default(i, 7)
        ucp(i)%URB_RfALB  =     lcz_ucp_default(i, 8)
        ucp(i)%URB_WaALB  =     lcz_ucp_default(i, 9)
        ucp(i)%URB_RdALB  =     lcz_ucp_default(i, 10)
        ucp(i)%URB_RfEMI  =     lcz_ucp_default(i, 11)
        ucp(i)%URB_WaEMI  =     lcz_ucp_default(i, 12)
        ucp(i)%URB_RdEMI  =     lcz_ucp_default(i, 13)
        ucp(i)%URB_RfHCAP =     lcz_ucp_default(i, 14)
        ucp(i)%URB_WaHCAP =     lcz_ucp_default(i, 15)
        ucp(i)%URB_RdHCAP =     lcz_ucp_default(i, 16)
        ucp(i)%URB_RfHCON =     lcz_ucp_default(i, 17)
        ucp(i)%URB_WaHCON =     lcz_ucp_default(i, 18)
        ucp(i)%URB_RdHCON =     lcz_ucp_default(i, 19)
        ucp(i)%AHF        =     lcz_ucp_default(i, 20)
      END DO

      ! Compute the UCP data
      DO i = 1,nr_lcz

        ! canyon albedo reduction factor, eq. 15
        psi_canyon = exp(-0.6 * ucp(i)%URB_H2W)
        ! total albedo reduction factor, eq. 14
        psi_bulk = psi_canyon * (1. - ucp(i)%URB_BLDFR) + ucp(i)%URB_BLDFR

        ! bulk shortwave albedo, using facet information, eq. 16
        alb_roof_snow = ucp(i)%URB_RfALB * (1. - snow_f) + alb_snow * snow_f
        alb_road_snow = ucp(i)%URB_RdALB * (1. - snow_f) + alb_snow * snow_f
        alb_wall_snow = ucp(i)%URB_WaALB * (1. - snow_f) + alb_snow * snow_f
        ucp(i)%URB_SALB_BK = (alb_road_snow + 2. * ucp(i)%URB_H2W * alb_wall_snow) &
                           / (1. + 2. * ucp(i)%URB_H2W) * psi_canyon * (1. - ucp(i)%URB_BLDFR) &
                           + alb_roof_snow * ucp(i)%URB_BLDFR

        ! bulk emissivity, using facet information, eq. 16
        emi_roof_snow = (1. - ucp(i)%URB_RfEMI) * (1. - snow_f) + (1. - emi_snow) * snow_f
        emi_road_snow = (1. - ucp(i)%URB_RdEMI) * (1. - snow_f) + (1. - emi_snow) * snow_f
        emi_wall_snow = (1. - ucp(i)%URB_WaEMI) * (1. - snow_f) + (1. - emi_snow) * snow_f
        ucp(i)%URB_EMIS_BK = 1. - ((emi_road_snow + 2. * ucp(i)%URB_H2W * emi_wall_snow) &
                           / (1. + 2. * ucp(i)%URB_H2W) * psi_canyon * (1. - ucp(i)%URB_BLDFR) &
                           + emi_roof_snow * ucp(i)%URB_BLDFR)

        ! bulk thermal albedo
        ucp(i)%URB_TALB_BK = 1 - ucp(i)%URB_EMIS_BK

        ! calculate surface area index from geometrical considerations, eq. 3
        SAI = (1. + 2. * ucp(i)%URB_H2W) * (1. - ucp(i)%URB_BLDFR) + ucp(i)%URB_BLDFR

        ! get mean heat capacity and conductivity, using eq. 10, 11 and 4.
        ucp(i)%URB_HCON = ((1. - ucp(i)%URB_BLDFR) / SAI) &
                        * (2. * ucp(i)%URB_H2W * ucp(i)%URB_WaHCON + ucp(i)%URB_RdHCON) &
                        + (ucp(i)%URB_BLDFR / SAI * ucp(i)%URB_RfHCON)
        ucp(i)%URB_HCAP = ((1. - ucp(i)%URB_BLDFR) / SAI) &
                        * (2. * ucp(i)%URB_H2W * ucp(i)%URB_WaHCAP + ucp(i)%URB_RdHCAP) &
                        + (ucp(i)%URB_BLDFR / SAI * ucp(i)%URB_RfHCAP)

        ! mean facet-level albedo and emissivity based on eq. 10
        ! Only added for testing and potential comparison with other models.
        ! These values are currently not used in TERRA_URB.
        ucp(i)%URB_EMIS_FL = ((1. - ucp(i)%URB_BLDFR) / SAI) &
                           * (2. * ucp(i)%URB_H2W * ucp(i)%URB_WaEMI + ucp(i)%URB_RdEMI) &
                           + (ucp(i)%URB_BLDFR / SAI * ucp(i)%URB_RfEMI)
        ucp(i)%URB_SALB_FL = ((1. - ucp(i)%URB_BLDFR) / SAI) &
                           * (2. * ucp(i)%URB_H2W * ucp(i)%URB_WaALB + ucp(i)%URB_RdALB) &
                           + (ucp(i)%URB_BLDFR / SAI * ucp(i)%URB_RfALB)
        ucp(i)%URB_TALB_FL = 1. - ucp(i)%URB_EMIS_FL

        ! for now, TERRA-URB only reads in one average facet-level albedo.
        ! The bulk calculation from eqs. 13 is done within TERRA_URB
        ! Therefore, the bulk value needs to be reversed back to a mean
        ! facet value, so that eq. 13 is solved for alb = alb_bulk / psi_bulk
        ! The same is done for the emissivity.
        ucp(i)%URB_SALB = ucp(i)%URB_SALB_BK / psi_bulk
        ucp(i)%URB_TALB = ucp(i)%URB_TALB_BK / psi_bulk
        ucp(i)%URB_EMIS = 1. - ucp(i)%URB_TALB

        ! Also add the thermal admittance (commented out in Matthias' code)
        ! ucp(i)%URB_TADM = (ucp(i)%URB_HCAP * ucp(i)%URB_HCON)**0.5

        ! is SAI weighting requested, according to Eq. 4?
        ! This is done within TERRA_URB, so no need to do for COSMO/CLM input files.
        IF (saiWeight) THEN
          ucp(i)%URB_HCON = ucp(i)%URB_HCON * SAI
          ucp(i)%URB_HCAP = ucp(i)%URB_HCAP * SAI
          ! ucp(i)%URB_TADM = ucp(i)%URB_TADM * SAI ! (commented out in Matthias' code)
        END IF

      END DO

      CALL logging%info('Exit routine: prepare_ucp_lookup')

    END SUBROUTINE prepare_ucp_lookup

!===============================================================================
    SUBROUTINE tu_test_all

      INTEGER :: i

      call terra_urb_to_extpar()

      DO i = 1,nr_lcz
        write(*,*) 'class_nr : ', ucp(i)%class_nr
        write(*,*) '    ISA        : ', ucp(i)%ISA
        write(*,*) '    FR_PAVED   : ', ucp(i)%FR_PAVED
        write(*,*) '    URBAN      : ', ucp(i)%URBAN
        write(*,*) '    URB_BLDFR  : ', ucp(i)%URB_BLDFR
        write(*,*) '    URB_BLDH   : ', ucp(i)%URB_BLDH
        write(*,*) '    URB_H2W    : ', ucp(i)%URB_H2W
        write(*,*) '    URB_RfALB  : ', ucp(i)%URB_RfALB
        write(*,*) '    URB_WaALB  : ', ucp(i)%URB_WaALB
        write(*,*) '    URB_RdALB  : ', ucp(i)%URB_RdALB
        write(*,*) '    URB_RfEMI  : ', ucp(i)%URB_RfEMI
        write(*,*) '    URB_WaEMI  : ', ucp(i)%URB_WaEMI
        write(*,*) '    URB_RdEMI  : ', ucp(i)%URB_RdEMI
        write(*,*) '    URB_RfHCAP : ', ucp(i)%URB_RfHCAP
        write(*,*) '    URB_WaHCAP : ', ucp(i)%URB_WaHCAP
        write(*,*) '    URB_RdHCAP : ', ucp(i)%URB_RdHCAP
        write(*,*) '    URB_RfHCON : ', ucp(i)%URB_RfHCON
        write(*,*) '    URB_WaHCON : ', ucp(i)%URB_WaHCON
        write(*,*) '    URB_RdHCON : ', ucp(i)%URB_RdHCON
        write(*,*) '    AHF        : ', ucp(i)%AHF
      END DO
    END SUBROUTINE tu_test_all

END MODULE mo_terra_urb
