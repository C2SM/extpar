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

  USE mo_kind,            ONLY: wp, i4
  USE mo_grid_structures, ONLY: target_grid_def
  USE mo_io_utilities,    ONLY: dim_meta_info, &
    &                           var_meta_info, &
    &                           vartype_real,  &
    &                           netcdf_put_var

  USE mo_logging

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: l_terra_urb,               &
    &       terra_urb_start,           &
    &       terra_urb_aggregate,       &
    &       terra_urb_def_fields_meta, &
    &       terra_urb_write_netcdf,    &
    &       terra_urb_end

  ! Modules' variables
  LOGICAL                     :: l_terra_urb = .FALSE.
  INTEGER(KIND=i4), PARAMETER :: nr_lcz = 10, & !< number of LCZ classes
                                 nr_ucp_values = 20 !< number of values defined in the tables

  !> List of urban parameters written to file:
  !> - ISA
  !> - AHF
  !> - FR_PAVED
  !> - URB_BLDFR
  !> - URB_BLDH
  !> - URB_H2W
  !> - URB_SALB
  !> - URB_TALB
  !> - URB_EMIS
  !> - URB_SALB_FL
  !> - URB_TALB_FL
  !> - URB_EMIS_FL
  !> - URB_SALB_BK
  !> - URB_TALB_BK
  !> - URB_EMIS_BK
  !> - URB_HCON
  !> - URB_HCAP
  REAL (KIND=wp), ALLOCATABLE :: tu_ISA        (:,:,:), &
       &                         tu_AHF        (:,:,:), &
       &                         tu_FR_PAVED   (:,:,:), &
       &                         tu_URB_BLDFR  (:,:,:), &
       &                         tu_URB_BLDH   (:,:,:), &
       &                         tu_URB_H2W    (:,:,:), &
       &                         tu_URB_SALB   (:,:,:), &
       &                         tu_URB_TALB   (:,:,:), &
       &                         tu_URB_EMIS   (:,:,:), &
       &                         tu_URB_SALB_FL(:,:,:), &
       &                         tu_URB_TALB_FL(:,:,:), &
       &                         tu_URB_EMIS_FL(:,:,:), &
       &                         tu_URB_SALB_BK(:,:,:), &
       &                         tu_URB_TALB_BK(:,:,:), &
       &                         tu_URB_EMIS_BK(:,:,:), &
       &                         tu_URB_HCON   (:,:,:), &
       &                         tu_URB_HCAP   (:,:,:)


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
                        URB_TALB_BK, &
                        URB_EMIS_BK, &
                        URB_HCON,    &
                        URB_HCAP,    &
                        URB_SALB_FL, &
                        URB_TALB_FL, &
                        URB_EMIS_FL, &
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

  TYPE(var_meta_info) :: tu_ISA_meta, &
       &                 tu_AHF_meta, &
       &                 tu_FR_PAVED_meta, &
       &                 tu_URB_BLDFR_meta, &
       &                 tu_URB_BLDH_meta, &
       &                 tu_URB_H2W_meta, &
       &                 tu_URB_SALB_meta, &
       &                 tu_URB_TALB_meta, &
       &                 tu_URB_EMIS_meta, &
       &                 tu_URB_SALB_FL_meta, &
       &                 tu_URB_TALB_FL_meta, &
       &                 tu_URB_EMIS_FL_meta, &
       &                 tu_URB_SALB_BK_meta, &
       &                 tu_URB_TALB_BK_meta, &
       &                 tu_URB_EMIS_BK_meta, &
       &                 tu_URB_HCON_meta, &
       &                 tu_URB_HCAP_meta

  CONTAINS

    !===========================================================================
    !> Main start routine
    !>
    SUBROUTINE terra_urb_start(tg)

      TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description

      CALL logging%info('Enter routine: terra_urb_start')

      ! Prepare the urban canopy data
      CALL tu_prepare_ucp_lookup()

      ! Allocate target fields
      CALL tu_allocate_target_fields(tg)

      CALL logging%info('Exit routine: terra_urb_start')

    END SUBROUTINE terra_urb_start

    !===========================================================================
    !> Main end routine
    !>
    SUBROUTINE terra_urb_end

      CALL logging%info('Enter routine: terra_urb_end')

      ! Deallocate target fields
      CALL tu_deallocate_target_fields()

      CALL logging%info('Exit routine: terra_urb_end')

    END SUBROUTINE terra_urb_end

    !===========================================================================
    !> Value assignment routine
    !>
    SUBROUTINE terra_urb_aggregate(ie,je,ke, lcz_nr, apix)

      INTEGER (KIND=i4), INTENT(IN) :: ie, je, ke, lcz_nr
      REAL (KIND=wp),    INTENT(IN) :: apix

      tu_ISA        (ie,je,ke) = tu_ISA        (ie,je,ke) + apix * ucp(lcz_nr)%ISA
      tu_AHF        (ie,je,ke) = tu_AHF        (ie,je,ke) + apix * ucp(lcz_nr)%AHF
      tu_FR_PAVED   (ie,je,ke) = tu_FR_PAVED   (ie,je,ke) + apix * ucp(lcz_nr)%FR_PAVED
      tu_URB_BLDFR  (ie,je,ke) = tu_URB_BLDFR  (ie,je,ke) + apix * ucp(lcz_nr)%URB_BLDFR
      tu_URB_BLDH   (ie,je,ke) = tu_URB_BLDH   (ie,je,ke) + apix * ucp(lcz_nr)%URB_BLDH
      tu_URB_H2W    (ie,je,ke) = tu_URB_H2W    (ie,je,ke) + apix * ucp(lcz_nr)%URB_H2W
      tu_URB_SALB   (ie,je,ke) = tu_URB_SALB   (ie,je,ke) + apix * ucp(lcz_nr)%URB_SALB
      tu_URB_TALB   (ie,je,ke) = tu_URB_TALB   (ie,je,ke) + apix * ucp(lcz_nr)%URB_TALB
      tu_URB_EMIS   (ie,je,ke) = tu_URB_EMIS   (ie,je,ke) + apix * ucp(lcz_nr)%URB_EMIS
      tu_URB_SALB_FL(ie,je,ke) = tu_URB_SALB_FL(ie,je,ke) + apix * ucp(lcz_nr)%URB_SALB_FL
      tu_URB_TALB_FL(ie,je,ke) = tu_URB_TALB_FL(ie,je,ke) + apix * ucp(lcz_nr)%URB_TALB_FL
      tu_URB_EMIS_FL(ie,je,ke) = tu_URB_EMIS_FL(ie,je,ke) + apix * ucp(lcz_nr)%URB_EMIS_FL
      tu_URB_SALB_BK(ie,je,ke) = tu_URB_SALB_BK(ie,je,ke) + apix * ucp(lcz_nr)%URB_SALB_BK
      tu_URB_TALB_BK(ie,je,ke) = tu_URB_TALB_BK(ie,je,ke) + apix * ucp(lcz_nr)%URB_TALB_BK
      tu_URB_EMIS_BK(ie,je,ke) = tu_URB_EMIS_BK(ie,je,ke) + apix * ucp(lcz_nr)%URB_EMIS_BK
      tu_URB_HCON   (ie,je,ke) = tu_URB_HCON   (ie,je,ke) + apix * ucp(lcz_nr)%URB_HCON
      tu_URB_HCAP   (ie,je,ke) = tu_URB_HCAP   (ie,je,ke) + apix * ucp(lcz_nr)%URB_HCAP

    END SUBROUTINE terra_urb_aggregate

    !===========================================================================
    !> Write the data to file
    !>
    SUBROUTINE terra_urb_write_netcdf(ncid, undefined)

      INTEGER (KIND=i4), INTENT(IN) :: ncid
      REAL(KIND=wp),     INTENT(IN) :: undefined !< value to indicate undefined grid elements

      CALL logging%info('Enter routine: terra_urb_write_netcdf')

      CALL netcdf_put_var(ncid,tu_ISA,         tu_ISA_meta,         undefined)
      CALL netcdf_put_var(ncid,tu_AHF,         tu_AHF_meta,         undefined)
      CALL netcdf_put_var(ncid,tu_FR_PAVED,    tu_FR_PAVED_meta,    undefined)
      CALL netcdf_put_var(ncid,tu_URB_BLDFR,   tu_URB_BLDFR_meta,   undefined)
      CALL netcdf_put_var(ncid,tu_URB_BLDH,    tu_URB_BLDH_meta,    undefined)
      CALL netcdf_put_var(ncid,tu_URB_H2W,     tu_URB_H2W_meta,     undefined)
      CALL netcdf_put_var(ncid,tu_URB_SALB,    tu_URB_SALB_meta,    undefined)
      CALL netcdf_put_var(ncid,tu_URB_TALB,    tu_URB_TALB_meta,    undefined)
      CALL netcdf_put_var(ncid,tu_URB_EMIS,    tu_URB_EMIS_meta,    undefined)
      CALL netcdf_put_var(ncid,tu_URB_SALB_FL, tu_URB_SALB_FL_meta, undefined)
      CALL netcdf_put_var(ncid,tu_URB_TALB_FL, tu_URB_TALB_FL_meta, undefined)
      CALL netcdf_put_var(ncid,tu_URB_EMIS_FL, tu_URB_EMIS_FL_meta, undefined)
      CALL netcdf_put_var(ncid,tu_URB_SALB_BK, tu_URB_SALB_BK_meta, undefined)
      CALL netcdf_put_var(ncid,tu_URB_TALB_BK, tu_URB_TALB_BK_meta, undefined)
      CALL netcdf_put_var(ncid,tu_URB_EMIS_BK, tu_URB_EMIS_BK_meta, undefined)
      CALL netcdf_put_var(ncid,tu_URB_HCON,    tu_URB_HCON_meta,    undefined)
      CALL netcdf_put_var(ncid,tu_URB_HCAP,    tu_URB_HCAP_meta,    undefined)

      CALL logging%info('Exit routine: terra_urb_write_netcdf')

    END SUBROUTINE terra_urb_write_netcdf

    !---------------------------------------------------------------------------
    !> Prepare the urban canopy data
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
    SUBROUTINE tu_prepare_ucp_lookup

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

    END SUBROUTINE tu_prepare_ucp_lookup

    !---------------------------------------------------------------------------
    !> Allocate fields for target grid
    !!
    SUBROUTINE tu_allocate_target_fields(tg)

      TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description

      INTEGER(KIND=i4)                  :: errorcode !< error status variable

      errorcode = 0

      CALL logging%info('Enter routine: tu_allocate_target_fields')

      ALLOCATE (tu_ISA(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant allocate the array tu_ISA',__FILE__,__LINE__)
      tu_ISA = 0.0
      ALLOCATE (tu_AHF(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant allocate the array tu_AHF',__FILE__,__LINE__)
      tu_AHF = 0.0
      ALLOCATE (tu_FR_PAVED(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant allocate the array tu_FR_PAVED',__FILE__,__LINE__)
      tu_FR_PAVED = 0.0
      ALLOCATE (tu_URB_BLDFR(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant allocate the array tu_URB_BLDFR',__FILE__,__LINE__)
      tu_URB_BLDFR = 0.0
      ALLOCATE (tu_URB_BLDH(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant allocate the array tu_URB_BLDH',__FILE__,__LINE__)
      tu_URB_BLDH = 0.0
      ALLOCATE (tu_URB_H2W(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant allocate the array tu_URB_H2W',__FILE__,__LINE__)
      tu_URB_H2W = 0.0
      ALLOCATE (tu_URB_SALB(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant allocate the array tu_URB_SALB',__FILE__,__LINE__)
      tu_URB_SALB = 0.0
      ALLOCATE (tu_URB_TALB(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant allocate the array tu_URB_TALB',__FILE__,__LINE__)
      tu_URB_TALB = 0.0
      ALLOCATE (tu_URB_EMIS(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant allocate the array tu_URB_EMIS',__FILE__,__LINE__)
      tu_URB_EMIS = 0.0
      ALLOCATE (tu_URB_SALB_FL(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant allocate the array tu_URB_SALB_FL',__FILE__,__LINE__)
      tu_URB_SALB_FL = 0.0
      ALLOCATE (tu_URB_TALB_FL(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant allocate the array tu_URB_TALB_FL',__FILE__,__LINE__)
      tu_URB_TALB_FL = 0.0
      ALLOCATE (tu_URB_EMIS_FL(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant allocate the array tu_URB_EMIS_FL',__FILE__,__LINE__)
      tu_URB_EMIS_FL = 0.0
      ALLOCATE (tu_URB_SALB_BK(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant allocate the array tu_URB_SALB_BK',__FILE__,__LINE__)
      tu_URB_SALB_BK = 0.0
      ALLOCATE (tu_URB_TALB_BK(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant allocate the array tu_URB_TALB_BK',__FILE__,__LINE__)
      tu_URB_TALB_BK = 0.0
      ALLOCATE (tu_URB_EMIS_BK(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant allocate the array tu_URB_EMIS_BK',__FILE__,__LINE__)
      tu_URB_EMIS_BK = 0.0
      ALLOCATE (tu_URB_HCON(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant allocate the array tu_URB_HCON',__FILE__,__LINE__)
      tu_URB_HCON = 0.0
      ALLOCATE (tu_URB_HCAP(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant allocate the array tu_URB_HCAP',__FILE__,__LINE__)
      tu_URB_HCAP = 0.0

      CALL logging%info('Exit routine: tu_allocate_target_fields')

    END SUBROUTINE tu_allocate_target_fields

    !---------------------------------------------------------------------------
    !> Deallocate fields for target grid
    !!
    SUBROUTINE tu_deallocate_target_fields

      INTEGER(KIND=i4) :: errorcode !< error status variable

      errorcode = 0

      CALL logging%info('Enter routine: tu_deallocate_target_fields')

      DEALLOCATE (tu_ISA, STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant deallocate the array tu_ISA',__FILE__,__LINE__)
      DEALLOCATE (tu_AHF, STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant deallocate the array tu_AHF',__FILE__,__LINE__)
      DEALLOCATE (tu_FR_PAVED, STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant deallocate the array tu_FR_PAVED',__FILE__,__LINE__)
      DEALLOCATE (tu_URB_BLDFR, STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant deallocate the array tu_URB_BLDFR',__FILE__,__LINE__)
      DEALLOCATE (tu_URB_BLDH, STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant deallocate the array tu_URB_BLDH',__FILE__,__LINE__)
      DEALLOCATE (tu_URB_H2W, STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant deallocate the array tu_URB_H2W',__FILE__,__LINE__)
      DEALLOCATE (tu_URB_SALB, STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant deallocate the array tu_URB_SALB',__FILE__,__LINE__)
      DEALLOCATE (tu_URB_TALB, STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant deallocate the array tu_URB_TALB',__FILE__,__LINE__)
      DEALLOCATE (tu_URB_EMIS, STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant deallocate the array tu_URB_EMIS',__FILE__,__LINE__)
      DEALLOCATE (tu_URB_SALB_FL, STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant deallocate the array tu_URB_SALB_FL',__FILE__,__LINE__)
      DEALLOCATE (tu_URB_TALB_FL, STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant deallocate the array tu_URB_TALB_FL',__FILE__,__LINE__)
      DEALLOCATE (tu_URB_EMIS_FL, STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant deallocate the array tu_URB_EMIS_FL',__FILE__,__LINE__)
      DEALLOCATE (tu_URB_SALB_BK, STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant deallocate the array tu_URB_SALB_BK',__FILE__,__LINE__)
      DEALLOCATE (tu_URB_TALB_BK, STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant deallocate the array tu_URB_TALB_BK',__FILE__,__LINE__)
      DEALLOCATE (tu_URB_EMIS_BK, STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant deallocate the array tu_URB_EMIS_BK',__FILE__,__LINE__)
      DEALLOCATE (tu_URB_HCON, STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant deallocate the array tu_URB_HCON',__FILE__,__LINE__)
      DEALLOCATE (tu_URB_HCAP, STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant deallocate the array tu_URB_HCAP',__FILE__,__LINE__)

      CALL logging%info('Exit routine: tu_deallocate_target_fields')

    END SUBROUTINE tu_deallocate_target_fields

    !---------------------------------------------------------------------------
    !> Write the metadata
    !>
    SUBROUTINE terra_urb_def_fields_meta(n_dim, diminfo, gridmp, coord, dataset)

      INTEGER,                     INTENT(IN) :: n_dim
      TYPE(dim_meta_info), TARGET, INTENT(IN) :: diminfo(:)
      CHARACTER (len=80),          INTENT(IN) :: gridmp, coord, dataset
      CHARACTER (len=1), PARAMETER            :: c_undef = "-" !< default character for undefined string

      tu_ISA_meta%varname = 'ISA'
      tu_ISA_meta%n_dim = n_dim
      tu_ISA_meta%diminfo => diminfo
      tu_ISA_meta%vartype = vartype_real
      tu_ISA_meta%standard_name = c_undef
      tu_ISA_meta%long_name = 'impervious surface area'
      tu_ISA_meta%shortName = 'ISA'
      tu_ISA_meta%stepType = 'instant'
      tu_ISA_meta%units =  c_undef
      tu_ISA_meta%grid_mapping = gridmp
      tu_ISA_meta%coordinates = coord
      tu_ISA_meta%data_set = dataset

      tu_AHF_meta%varname = 'AHF'
      tu_AHF_meta%n_dim = n_dim
      tu_AHF_meta%diminfo => diminfo
      tu_AHF_meta%vartype = vartype_real
      tu_AHF_meta%standard_name = c_undef
      tu_AHF_meta%long_name = 'Anthropogenic heat flux'
      tu_AHF_meta%shortName = 'GFLUX'
      tu_AHF_meta%stepType = 'instant'
      tu_AHF_meta%units =  'W m-2'
      tu_AHF_meta%grid_mapping = gridmp
      tu_AHF_meta%coordinates = coord
      tu_AHF_meta%data_set = dataset

      tu_FR_PAVED_meta%varname = 'FR_PAVED'
      tu_FR_PAVED_meta%n_dim = n_dim
      tu_FR_PAVED_meta%diminfo => diminfo
      tu_FR_PAVED_meta%vartype = vartype_real
      tu_FR_PAVED_meta%standard_name = c_undef
      tu_FR_PAVED_meta%long_name = 'Fraction of impervious surface area'
      tu_FR_PAVED_meta%shortName = 'FR_PAVED'
      tu_FR_PAVED_meta%stepType = 'instant'
      tu_FR_PAVED_meta%units =  c_undef
      tu_FR_PAVED_meta%grid_mapping = gridmp
      tu_FR_PAVED_meta%coordinates = coord
      tu_FR_PAVED_meta%data_set = dataset

      tu_URB_BLDFR_meta%varname = 'URB_BLDFR'
      tu_URB_BLDFR_meta%n_dim = n_dim
      tu_URB_BLDFR_meta%diminfo => diminfo
      tu_URB_BLDFR_meta%vartype = vartype_real
      tu_URB_BLDFR_meta%standard_name = c_undef
      tu_URB_BLDFR_meta%long_name = 'Urban building fraction'
      tu_URB_BLDFR_meta%shortName = 'URB_BLDFR'
      tu_URB_BLDFR_meta%stepType = 'instant'
      tu_URB_BLDFR_meta%units =  c_undef
      tu_URB_BLDFR_meta%grid_mapping = gridmp
      tu_URB_BLDFR_meta%coordinates = coord
      tu_URB_BLDFR_meta%data_set = dataset

      tu_URB_BLDH_meta%varname = 'URB_BLDH'
      tu_URB_BLDH_meta%n_dim = n_dim
      tu_URB_BLDH_meta%diminfo => diminfo
      tu_URB_BLDH_meta%vartype = vartype_real
      tu_URB_BLDH_meta%standard_name = c_undef
      tu_URB_BLDH_meta%long_name = 'Urban building height'
      tu_URB_BLDH_meta%shortName = 'URB_BLDH'
      tu_URB_BLDH_meta%stepType = 'instant'
      tu_URB_BLDH_meta%units =  c_undef
      tu_URB_BLDH_meta%grid_mapping = gridmp
      tu_URB_BLDH_meta%coordinates = coord
      tu_URB_BLDH_meta%data_set = dataset

      tu_URB_H2W_meta%varname = 'URB_H2W'
      tu_URB_H2W_meta%n_dim = n_dim
      tu_URB_H2W_meta%diminfo => diminfo
      tu_URB_H2W_meta%vartype = vartype_real
      tu_URB_H2W_meta%standard_name = c_undef
      tu_URB_H2W_meta%long_name = 'Urban canyon height to width'
      tu_URB_H2W_meta%shortName = 'URB_H2W'
      tu_URB_H2W_meta%stepType = 'instant'
      tu_URB_H2W_meta%units =  c_undef
      tu_URB_H2W_meta%grid_mapping = gridmp
      tu_URB_H2W_meta%coordinates = coord
      tu_URB_H2W_meta%data_set = dataset

      tu_URB_SALB_meta%varname = 'URB_SALB'
      tu_URB_SALB_meta%n_dim = n_dim
      tu_URB_SALB_meta%diminfo => diminfo
      tu_URB_SALB_meta%vartype = vartype_real
      tu_URB_SALB_meta%standard_name = c_undef
      tu_URB_SALB_meta%long_name = 'Urban shortwave albedo'
      tu_URB_SALB_meta%shortName = 'URB_SALB'
      tu_URB_SALB_meta%stepType = 'instant'
      tu_URB_SALB_meta%units =  c_undef
      tu_URB_SALB_meta%grid_mapping = gridmp
      tu_URB_SALB_meta%coordinates = coord
      tu_URB_SALB_meta%data_set = dataset

      tu_URB_TALB_meta%varname = 'URB_TALB'
      tu_URB_TALB_meta%n_dim = n_dim
      tu_URB_TALB_meta%diminfo => diminfo
      tu_URB_TALB_meta%vartype = vartype_real
      tu_URB_TALB_meta%standard_name = c_undef
      tu_URB_TALB_meta%long_name = 'Urban thermal albedo'
      tu_URB_TALB_meta%shortName = 'URB_TALB'
      tu_URB_TALB_meta%stepType = 'instant'
      tu_URB_TALB_meta%units =  c_undef
      tu_URB_TALB_meta%grid_mapping = gridmp
      tu_URB_TALB_meta%coordinates = coord
      tu_URB_TALB_meta%data_set = dataset

      tu_URB_EMIS_meta%varname = 'URB_EMIS'
      tu_URB_EMIS_meta%n_dim = n_dim
      tu_URB_EMIS_meta%diminfo => diminfo
      tu_URB_EMIS_meta%vartype = vartype_real
      tu_URB_EMIS_meta%standard_name = c_undef
      tu_URB_EMIS_meta%long_name = 'Urban emissivity'
      tu_URB_EMIS_meta%shortName = 'URB_EMIS'
      tu_URB_EMIS_meta%stepType = 'instant'
      tu_URB_EMIS_meta%units =  c_undef
      tu_URB_EMIS_meta%grid_mapping = gridmp
      tu_URB_EMIS_meta%coordinates = coord
      tu_URB_EMIS_meta%data_set = dataset

      tu_URB_SALB_FL_meta%varname = 'URB_SALB_FL'
      tu_URB_SALB_FL_meta%n_dim = n_dim
      tu_URB_SALB_FL_meta%diminfo => diminfo
      tu_URB_SALB_FL_meta%vartype = vartype_real
      tu_URB_SALB_FL_meta%standard_name = c_undef
      tu_URB_SALB_FL_meta%long_name = 'Urban shortwave albedo (mean facet-level)'
      tu_URB_SALB_FL_meta%shortName = 'URB_SALB_FL'
      tu_URB_SALB_FL_meta%stepType = 'instant'
      tu_URB_SALB_FL_meta%units =  c_undef
      tu_URB_SALB_FL_meta%grid_mapping = gridmp
      tu_URB_SALB_FL_meta%coordinates = coord
      tu_URB_SALB_FL_meta%data_set = dataset

      tu_URB_TALB_FL_meta%varname = 'URB_TALB_FL'
      tu_URB_TALB_FL_meta%n_dim = n_dim
      tu_URB_TALB_FL_meta%diminfo => diminfo
      tu_URB_TALB_FL_meta%vartype = vartype_real
      tu_URB_TALB_FL_meta%standard_name = c_undef
      tu_URB_TALB_FL_meta%long_name = 'Urban thermal albedo (mean facet-level)'
      tu_URB_TALB_FL_meta%shortName = 'URB_TALB_FL'
      tu_URB_TALB_FL_meta%stepType = 'instant'
      tu_URB_TALB_FL_meta%units =  c_undef
      tu_URB_TALB_FL_meta%grid_mapping = gridmp
      tu_URB_TALB_FL_meta%coordinates = coord
      tu_URB_TALB_FL_meta%data_set = dataset

      tu_URB_EMIS_FL_meta%varname = 'URB_EMIS_FL'
      tu_URB_EMIS_FL_meta%n_dim = n_dim
      tu_URB_EMIS_FL_meta%diminfo => diminfo
      tu_URB_EMIS_FL_meta%vartype = vartype_real
      tu_URB_EMIS_FL_meta%standard_name = c_undef
      tu_URB_EMIS_FL_meta%long_name = 'Urban emissivity (mean facet-level)'
      tu_URB_EMIS_FL_meta%shortName = 'URB_EMIS_FL'
      tu_URB_EMIS_FL_meta%stepType = 'instant'
      tu_URB_EMIS_FL_meta%units =  c_undef
      tu_URB_EMIS_FL_meta%grid_mapping = gridmp
      tu_URB_EMIS_FL_meta%coordinates = coord
      tu_URB_EMIS_FL_meta%data_set = dataset

      tu_URB_SALB_BK_meta%varname = 'URB_SALB_BK'
      tu_URB_SALB_BK_meta%n_dim = n_dim
      tu_URB_SALB_BK_meta%diminfo => diminfo
      tu_URB_SALB_BK_meta%vartype = vartype_real
      tu_URB_SALB_BK_meta%standard_name = c_undef
      tu_URB_SALB_BK_meta%long_name = 'Urban shortwave albedo (bulk)'
      tu_URB_SALB_BK_meta%shortName = 'URB_SALB_BK'
      tu_URB_SALB_BK_meta%stepType = 'instant'
      tu_URB_SALB_BK_meta%units =  c_undef
      tu_URB_SALB_BK_meta%grid_mapping = gridmp
      tu_URB_SALB_BK_meta%coordinates = coord
      tu_URB_SALB_BK_meta%data_set = dataset

      tu_URB_TALB_BK_meta%varname = 'URB_TALB_BK'
      tu_URB_TALB_BK_meta%n_dim = n_dim
      tu_URB_TALB_BK_meta%diminfo => diminfo
      tu_URB_TALB_BK_meta%vartype = vartype_real
      tu_URB_TALB_BK_meta%standard_name = c_undef
      tu_URB_TALB_BK_meta%long_name = 'Urban thermal albedo (bulk)'
      tu_URB_TALB_BK_meta%shortName = 'URB_TALB_BK'
      tu_URB_TALB_BK_meta%stepType = 'instant'
      tu_URB_TALB_BK_meta%units =  c_undef
      tu_URB_TALB_BK_meta%grid_mapping = gridmp
      tu_URB_TALB_BK_meta%coordinates = coord
      tu_URB_TALB_BK_meta%data_set = dataset

      tu_URB_EMIS_BK_meta%varname = 'URB_EMIS_BK'
      tu_URB_EMIS_BK_meta%n_dim = n_dim
      tu_URB_EMIS_BK_meta%diminfo => diminfo
      tu_URB_EMIS_BK_meta%vartype = vartype_real
      tu_URB_EMIS_BK_meta%standard_name = c_undef
      tu_URB_EMIS_BK_meta%long_name = 'Urban emissivity (bulk)'
      tu_URB_EMIS_BK_meta%shortName = 'URB_EMIS_BK'
      tu_URB_EMIS_BK_meta%stepType = 'instant'
      tu_URB_EMIS_BK_meta%units =  c_undef
      tu_URB_EMIS_BK_meta%grid_mapping = gridmp
      tu_URB_EMIS_BK_meta%coordinates = coord
      tu_URB_EMIS_BK_meta%data_set = dataset

      tu_URB_HCON_meta%varname = 'URB_HCON'
      tu_URB_HCON_meta%n_dim = n_dim
      tu_URB_HCON_meta%diminfo => diminfo
      tu_URB_HCON_meta%vartype = vartype_real
      tu_URB_HCON_meta%standard_name = c_undef
      tu_URB_HCON_meta%long_name = 'Urban mean heat conductivity'
      tu_URB_HCON_meta%shortName = 'URB_HCON'
      tu_URB_HCON_meta%stepType = 'instant'
      tu_URB_HCON_meta%units =  c_undef
      tu_URB_HCON_meta%grid_mapping = gridmp
      tu_URB_HCON_meta%coordinates = coord
      tu_URB_HCON_meta%data_set = dataset

      tu_URB_HCAP_meta%varname = 'URB_HCAP'
      tu_URB_HCAP_meta%n_dim = n_dim
      tu_URB_HCAP_meta%diminfo => diminfo
      tu_URB_HCAP_meta%vartype = vartype_real
      tu_URB_HCAP_meta%standard_name = c_undef
      tu_URB_HCAP_meta%long_name = 'Urban mean heat capacity'
      tu_URB_HCAP_meta%shortName = 'URB_HCAP'
      tu_URB_HCAP_meta%stepType = 'instant'
      tu_URB_HCAP_meta%units =  c_undef
      tu_URB_HCAP_meta%grid_mapping = gridmp
      tu_URB_HCAP_meta%coordinates = coord
      tu_URB_HCAP_meta%data_set = dataset

    END SUBROUTINE terra_urb_def_fields_meta


END MODULE mo_terra_urb
