!+ Fortran Module with lookup-tables for the GLCC data
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran Module with lookup-tables for the GLCC data
!> \author Hermann Asensio
!!
!! Description:
!! The GLCC dataset contains the following land use classification scheme
!!
!! - operational settings of GME and COSMO (Heise, 2005)
!! - experimantal settings, lookup-table analog to Ecoclimap (Masson 2003), use only in combination with NDVI data!

MODULE mo_glcc_lookup_tables

 !> kind parameters are defined in MODULE data_parameters
  USE mo_kind, ONLY: wp
  USE mo_kind, ONLY: i4
  USE mo_kind, ONLY: i4

  !> abort_extpar defined in MODULE utilities_extpar
  USE mo_utilities_extpar, ONLY: abort_extpar
  USE mo_io_units,         ONLY: filename_max



                           
IMPLICIT NONE


PUBLIC :: init_glcc_lookup_tables
PUBLIC :: get_name_glcc_lookup_tables
PUBLIC :: glcc_look_up
PUBLIC :: glcc_legend

PUBLIC :: nclass_glcc
PUBLIC :: ilookup_table_glcc
PUBLIC :: name_lookup_table_glcc

PUBLIC :: i_gme_lookup_table_glcc, i_cosmo_lookup_table_glcc, i_experimental_lookup_table_glcc
PUBLIC :: z0_lt_glcc, lnz0_lt_glcc, plc_mn_lt_glcc, plc_mx_lt_glcc 
PUBLIC :: lai_mn_lt_glcc, lai_mx_lt_glcc, rd_lt_glcc, emiss_lt_glcc, rs_min_lt_glcc         



INTEGER (KIND=i4), PARAMETER :: nclass_glcc = 24 !< GLCC has 24 classes for the land use description

INTEGER (KIND=i4), PARAMETER :: i_gme_lookup_table_glcc = 1 
                                !< lookup_table for GLCC land use classes from the operational settings of GME (Ritter, 2007)
INTEGER (KIND=i4), PARAMETER :: i_cosmo_lookup_table_glcc = 2 
                                !< lookup_table for GLCC land use classes from the operational settings of COSMO (Heise, 2005)
INTEGER (KIND=i4), PARAMETER :: i_experimental_lookup_table_glcc = 3 
!< lookup_table for GLCC land use classes for experimental setting, analog to look-up tables of ECOCLIMAP (Masson 2003)

INTEGER (KIND=i4) :: ilookup_table_glcc !< integer switch to choose a lookup table
CHARACTER (LEN=filename_max) :: name_lookup_table_glcc !< name of lookup table




REAL (KIND=wp) :: z0_lt_glcc(nclass_glcc)      !< lookup table landuse class to roughness length [m]
REAL (KIND=wp) :: lnz0_lt_glcc(nclass_glcc)    !< corresponding natural logarithm of z0c_gme_o
REAL (KIND=wp) :: plc_mn_lt_glcc(nclass_glcc)  !< lookup table landuse class to minimal plant cover
REAL (KIND=wp) :: plc_mx_lt_glcc(nclass_glcc)  !< lookup table landuse class to maximal plant cover
REAL (KIND=wp) :: lai_mn_lt_glcc(nclass_glcc)  !< lookup table landuse class to minimal leaf area index
REAL (KIND=wp) :: lai_mx_lt_glcc(nclass_glcc)  !< lookup table landuse class to maximal leaf area index
REAL (KIND=wp) :: rd_lt_glcc(nclass_glcc)      !< lookup table landuse class to root depth [m]
REAL (KIND=wp) :: emiss_lt_glcc(nclass_glcc)   !< lookup table landuse class to surface thermal emissivity
REAL (KIND=wp) :: rs_min_lt_glcc(nclass_glcc)  !< lookup table landuse class to minimal stomata resistance

!---------------------------------------------------------------------------------------------- 
!----------------------------------------------------------------------------------------------
REAL (KIND=wp) :: z0c_gme_o(nclass_glcc)  = (/ &       !< lookup table landuse class to roughness length [m]
  &         1.00,       &       ! class urban and built-up land
  &          0.10,       &       ! dryland cropland and pasture
  &          0.10,       &       ! irrigated cropland and pasture
  &          0.10,       &       ! mixed dryland/irrigated .....
  &          0.07,       &       ! cropland/grassland mosaic
  &          0.25,       &       ! cropland/woodland  mosaic
  &          0.03,       &       ! grassland
  &          0.20,       &       ! shrubland
  &          0.15,       &       ! mixed shrubland/grassland
  &          0.15,       &       ! savanna
  &          1.00,       &       ! decidous broadleaf forest
  &          1.00,       &       ! decidous needleleaf forest
  &          1.00,       &       ! evergreen broadleaf forest
  &          1.00,       &       ! evergreen needleleaf forest
  &          1.00,       &       ! mixed forest
  &          0.0002,     &       ! water bodies
  &          0.05,       &       ! herbaceous wetland
  &          0.20,       &       ! wooded     wetland
  &          0.05,       &       ! barren or sparsely vegetated
  &          0.05,       &       ! herbaceous  tundra
  &          0.20,       &       ! wooded      tundra
  &          0.10,       &       ! mixed       tundra
  &          0.03,       &       ! bare ground tundra
  &          0.01   /)           ! snow or ice

REAL (KIND=wp) :: lnz0c_gme_o(nclass_glcc)    !< corresponding natural logarithm of z0c_gme_o


REAL (KIND=wp) :: zplcmnc_gme_o(nclass_glcc) = (/ &      !< lookup table landuse class to minimal plant cover
  &          0.05,      &           ! class urban and built-up land
  &          0.45,      &           ! dryland cropland and pasture
  &          0.50,      &           ! irrigated cropland and pasture
  &          0.45,      &           ! mixed dryland/irrigated .....
  &          0.45,      &           ! cropland/grassland mosaic
  &          0.45,      &           ! cropland/woodland  mosaic
  &          1.00,      &           ! grassland
  &          0.10,      &           ! shrubland
  &          0.10,      &           ! mixed shrubland/grassland
  &          0.20,      &           ! savanna
  &          0.00,      &           ! decidous broadleaf forest
  &          0.00,      &           ! decidous needleleaf forest
  &          1.00,      &           ! evergreen broadleaf forest
  &          1.00,      &           ! evergreen needleleaf forest
  &          0.50,      &           ! mixed forest
  &          0.00,      &           ! water bodies
  &          0.40,      &           ! herbaceous wetland
  &          0.10,      &           ! wooded     wetland
  &          0.02,      &           ! barren or sparsely vegetated
  &          0.00,      &           ! herbaceous  tundra
  &          0.20,      &           ! wooded      tundra
  &          0.10,      &           ! mixed       tundra
  &          0.00,      &           ! bare ground tundra
  &          0.00   /)              ! snow or ice



REAL (KIND=wp) :: zplcmxc_gme_o(nclass_glcc) = (/ &     !< lookup table landuse class to maximal plant cover
  &           0.20, &       ! class urban and built-up land
  &           0.70, &       ! dryland cropland and pasture
  &           0.90, &       ! irrigated cropland and pasture
  &           0.80, &       ! mixed dryland/irrigated .....
  &           0.90, &       ! cropland/grassland mosaic
  &           0.80, &       ! cropland/woodland  mosaic
  &           0.90, &       ! grassland
  &           0.80, &       ! shrubland
  &           0.80, &       ! mixed shrubland/grassland
  &           0.80, &       ! savanna
  &           0.90, &       ! decidous broadleaf forest
  &           0.90, &       ! decidous needleleaf forest
  &           0.80, &       ! evergreen broadleaf forest
  &           0.80, &       ! evergreen needleleaf forest
  &           0.90, &       ! mixed forest
  &           0.00, &       ! water bodies
  &           0.80, &       ! herbaceous wetland
  &           0.80, &       ! wooded     wetland
  &           0.05, &       ! barren or sparsely vegetated
  &           0.50, &       ! herbaceous  tundra
  &           0.50, &       ! wooded      tundra
  &           0.50, &       ! mixed       tundra
  &           0.00, &       ! bare ground tundra
  &           0.00  /)      ! snow or ice



REAL (KIND=wp) :: zlaimnc_gme_o(nclass_glcc) = (/ &      !< lookup table landuse class to minimal leaf area index
  &           0.10,  &               ! class urban and built-up land
  &           0.20,  &               ! dryland cropland and pasture
  &           0.20,  &               ! irrigated cropland and pasture
  &           0.20,  &               ! mixed dryland/irrigated .....
  &           0.35,  &               ! cropland/grassland mosaic
  &           1.20,  &               ! cropland/woodland  mosaic
  &           0.50,  &               ! grassland
  &           0.10,  &               ! shrubland
  &           0.10,  &               ! mixed shrubland/grassland
  &           1.00,  &               ! savanna
  &           0.00,  &               ! decidous broadleaf forest
  &           0.00,  &               ! decidous needleleaf forest
  &           9.00,  &               ! evergreen broadleaf forest
  &           8.00,  &               ! evergreen needleleaf forest
  &           2.25,  &               ! mixed forest
  &           0.00,  &               ! water bodies
  &           1.00,  &               ! herbaceous wetland
  &           1.00,  &               ! wooded     wetland
  &           0.50,  &               ! barren or sparsely vegetated
  &           0.00,  &               ! herbaceous  tundra
  &           0.50,  &               ! wooded      tundra
  &           0.25,  &               ! mixed       tundra
  &           0.00,  &               ! bare ground tundra
  &           0.00   /)              ! snow or ice


REAL (KIND=wp) :: zlaimxc_gme_o(nclass_glcc) = (/ &      !< lookup table landuse class to maximal leaf area index
  &           1.00,  &      ! class urban and built-up land
  &           3.30,  &      ! dryland cropland and pasture
  &           3.00,  &      ! irrigated cropland and pasture
  &           3.10,  &      ! mixed dryland/irrigated .....
  &           3.00,  &      ! cropland/grassland mosaic
  &           3.50,  &      ! cropland/woodland  mosaic
  &           2.00,  &      ! grassland
  &           2.50,  &      ! shrubland
  &           2.25,  &      ! mixed shrubland/grassland
  &           2.00,  &      ! savanna
  &           6.00,  &      ! decidous broadleaf forest
  &           5.00,  &      ! decidous needleleaf forest
  &           5.00,  &      ! evergreen broadleaf forest
  &           5.00,  &      ! evergreen needleleaf forest
  &           5.00,  &      ! mixed forest
  &           0.00,  &      ! water bodies
  &           2.00,  &      ! herbaceous wetland
  &           2.50,  &      ! wooded     wetland
  &           0.60,  &      ! barren or sparsely vegetated
  &           1.00,  &      ! herbaceous  tundra
  &           1.60,  &      ! wooded      tundra
  &           1.30,  &      ! mixed       tundra
  &           0.00,  &      ! bare ground tundra
  &           0.00  /)      ! snow or ice

REAL (KIND=wp) :: zrd_gme_o(nclass_glcc)  = (/ &         !< lookup table landuse class to root depth [m]
  &          0.60,  &            ! class urban and built-up land
  &          1.00,  &            ! dryland cropland and pasture
  &          0.60,  &            ! irrigated cropland and pasture
  &          0.80,  &            ! mixed dryland/irrigated .....
  &          1.00,  &            ! cropland/grassland mosaic
  &          1.00,  &            ! cropland/woodland  mosaic
  &          0.60,  &            ! grassland
  &          1.00,  &            ! shrubland
  &          1.00,  &            ! mixed shrubland/grassland
  &          2.00,  &            ! savanna
  &          1.00,  &            ! decidous broadleaf forest
  &          0.60,  &            ! decidous needleleaf forest
  &          1.00,  &            ! evergreen broadleaf forest
  &          0.60,  &            ! evergreen needleleaf forest
  &          0.80,  &            ! mixed forest
  &          0.00,  &            ! water bodies
  &          0.40,  &            ! herbaceous wetland
  &          0.40,  &            ! wooded     wetland
  &          0.30,  &            ! barren or sparsely vegetated
  &          0.10,  &            ! herbaceous  tundra
  &          0.10,  &            ! wooded      tundra
  &          0.10,  &            ! mixed       tundra
  &          0.00,  &            ! bare ground tundra
  &          0.00  /)            ! snow or ice




REAL (KIND=wp) :: zemiss_gme_o(nclass_glcc) = (/ &       !< lookup table landuse class to surface thermal emissivity
  &          0.960, &            ! class urban and built-up land
  &          0.985, &            ! dryland cropland and pasture
  &          0.990, &            ! irrigated cropland and pasture
  &          0.990, &            ! mixed dryland/irrigated .....
  &          0.990, &            ! cropland/grassland mosaic
  &          0.990, &            ! cropland/woodland  mosaic
  &          0.993, &            ! grassland
  &          0.985, &            ! shrubland
  &          0.985, &            ! mixed shrubland/grassland
  &          0.993, &            ! savanna
  &          0.990, &            ! decidous broadleaf forest
  &          0.990, &            ! decidous needleleaf forest
  &          0.996, &            ! evergreen broadleaf forest
  &          0.996, &            ! evergreen needleleaf forest
  &          0.993, &            ! mixed forest
  &          0.991, &            ! water bodies
  &          0.992, &            ! herbaceous wetland
  &          0.992, &            ! wooded     wetland
  &          0.950, &            ! barren or sparsely vegetated
  &          0.990, &            ! herbaceous  tundra
  &          0.990, &            ! wooded      tundra
  &          0.990, &            ! mixed       tundra
  &          0.990, &            ! bare ground tundra
  &          0.9999 /)            ! snow or ice


REAL (KIND=wp) :: zrs_min_gme_o(nclass_glcc) = (/ &      !< lookup table landuse class to minimal stomata resistance
  &          150.0, &       ! class urban and built-up land
  &          180.0, &       ! dryland cropland and pasture
  &          180.0, &       ! irrigated cropland and pasture
  &          180.0, &       ! mixed dryland/irrigated .....
  &          180.0, &       ! cropland/grassland mosaic
  &          200.0, &       ! cropland/woodland  mosaic
  &          110.0, &       ! grassland
  &          225.0, &       ! shrubland
  &          150.0, &       ! mixed shrubland/grassland
  &          100.0, &       ! savanna
  &          240.0, &       ! decidous broadleaf forest
  &          500.0, &       ! decidous needleleaf forest
  &          240.0, &       ! evergreen broadleaf forest
  &          500.0, &       ! evergreen needleleaf forest
  &          350.0, &       ! mixed forest
  &           0.00, &       ! water bodies
  &          110.0, &       ! herbaceous wetland
  &          110.0, &       ! wooded     wetland
  &          350.0, &       ! barren or sparsely vegetated
  &           80.0, &       ! herbaceous  tundra
  &           80.0, &       ! wooded      tundra
  &           80.0, &       ! mixed       tundra
  &           0.00,  &      ! bare ground tundra
  &           0.00  /)      ! snow or ice


!---------------------------------------------------------------------------------------------- 
!----------------------------------------------------------------------------------------------
REAL (KIND=wp) :: z0c_cosmo_o(nclass_glcc)  = (/ &       !< lookup table landuse class to roughness length [m]
  &         1.00,       &       ! class urban and built-up land
  &          0.10,       &       ! dryland cropland and pasture
  &          0.10,       &       ! irrigated cropland and pasture
  &          0.10,       &       ! mixed dryland/irrigated .....
  &          0.07,       &       ! cropland/grassland mosaic
  &          0.25,       &       ! cropland/woodland  mosaic
  &          0.03,       &       ! grassland
  &          0.20,       &       ! shrubland
  &          0.15,       &       ! mixed shrubland/grassland
  &          0.15,       &       ! savanna
  &          1.00,       &       ! decidous broadleaf forest
  &          1.00,       &       ! decidous needleleaf forest
  &          1.00,       &       ! evergreen broadleaf forest
  &          1.00,       &       ! evergreen needleleaf forest
  &          1.00,       &       ! mixed forest
  &          0.0002,     &       ! water bodies
  &          0.05,       &       ! herbaceous wetland
  &          0.20,       &       ! wooded     wetland
  &          0.05,       &       ! barren or sparsely vegetated
  &          0.05,       &       ! herbaceous  tundra
  &          0.20,       &       ! wooded      tundra
  &          0.10,       &       ! mixed       tundra
  &          0.03,       &       ! bare ground tundra
  &          0.01   /)           ! snow or ice

REAL (KIND=wp) :: lnz0c_cosmo_o(nclass_glcc)    !< corresponding natural logarithm of z0c_cosmo_o


REAL (KIND=wp) :: zplcmnc_cosmo_o(nclass_glcc) = (/ &      !< lookup table landuse class to minimal plant cover
  &          0.05,      &           ! class urban and built-up land
  &          0.45,      &           ! dryland cropland and pasture
  &          0.50,      &           ! irrigated cropland and pasture
  &          0.45,      &           ! mixed dryland/irrigated .....
  &          0.45,      &           ! cropland/grassland mosaic
  &          0.45,      &           ! cropland/woodland  mosaic
  &          1.00,      &           ! grassland
  &          0.10,      &           ! shrubland
  &          0.10,      &           ! mixed shrubland/grassland
  &          0.20,      &           ! savanna
  &          0.00,      &           ! decidous broadleaf forest
  &          0.00,      &           ! decidous needleleaf forest
  &          1.00,      &           ! evergreen broadleaf forest
  &          1.00,      &           ! evergreen needleleaf forest
  &          0.50,      &           ! mixed forest
  &          0.00,      &           ! water bodies
  &          0.40,      &           ! herbaceous wetland
  &          0.10,      &           ! wooded     wetland
  &          0.02,      &           ! barren or sparsely vegetated
  &          0.00,      &           ! herbaceous  tundra
  &          0.20,      &           ! wooded      tundra
  &          0.10,      &           ! mixed       tundra
  &          0.00,      &           ! bare ground tundra
  &          0.00   /)              ! snow or ice



REAL (KIND=wp) :: zplcmxc_cosmo_o(nclass_glcc) = (/ &     !< lookup table landuse class to maximal plant cover
  &           0.20, &       ! class urban and built-up land
  &           0.70, &       ! dryland cropland and pasture
  &           0.90, &       ! irrigated cropland and pasture
  &           0.80, &       ! mixed dryland/irrigated .....
  &           0.90, &       ! cropland/grassland mosaic
  &           0.80, &       ! cropland/woodland  mosaic
  &           0.90, &       ! grassland
  &           0.80, &       ! shrubland
  &           0.80, &       ! mixed shrubland/grassland
  &           0.80, &       ! savanna
  &           0.90, &       ! decidous broadleaf forest
  &           0.90, &       ! decidous needleleaf forest
  &           0.80, &       ! evergreen broadleaf forest
  &           0.80, &       ! evergreen needleleaf forest
  &           0.90, &       ! mixed forest
  &           0.00, &       ! water bodies
  &           0.80, &       ! herbaceous wetland
  &           0.80, &       ! wooded     wetland
  &           0.05, &       ! barren or sparsely vegetated
  &           0.50, &       ! herbaceous  tundra
  &           0.50, &       ! wooded      tundra
  &           0.50, &       ! mixed       tundra
  &           0.00, &       ! bare ground tundra
  &           0.00  /)      ! snow or ice



REAL (KIND=wp) :: zlaimnc_cosmo_o(nclass_glcc) = (/ &      !< lookup table landuse class to minimal leaf area index
  &           0.10,  &               ! class urban and built-up land
  &           0.20,  &               ! dryland cropland and pasture
  &           0.20,  &               ! irrigated cropland and pasture
  &           0.20,  &               ! mixed dryland/irrigated .....
  &           0.35,  &               ! cropland/grassland mosaic
  &           1.20,  &               ! cropland/woodland  mosaic
  &           0.50,  &               ! grassland
  &           0.10,  &               ! shrubland
  &           0.10,  &               ! mixed shrubland/grassland
  &           1.00,  &               ! savanna
  &           0.00,  &               ! decidous broadleaf forest
  &           0.00,  &               ! decidous needleleaf forest
  &           9.00,  &               ! evergreen broadleaf forest
  &           8.00,  &               ! evergreen needleleaf forest
  &           2.25,  &               ! mixed forest
  &           0.00,  &               ! water bodies
  &           1.00,  &               ! herbaceous wetland
  &           1.00,  &               ! wooded     wetland
  &           0.50,  &               ! barren or sparsely vegetated
  &           0.00,  &               ! herbaceous  tundra
  &           0.50,  &               ! wooded      tundra
  &           0.25,  &               ! mixed       tundra
  &           0.00,  &               ! bare ground tundra
  &           0.00   /)              ! snow or ice


REAL (KIND=wp) :: zlaimxc_cosmo_o(nclass_glcc) = (/ &      !< lookup table landuse class to maximal leaf area index
  &           1.00,  &      ! class urban and built-up land
  &           3.30,  &      ! dryland cropland and pasture
  &           3.00,  &      ! irrigated cropland and pasture
  &           3.10,  &      ! mixed dryland/irrigated .....
  &           3.00,  &      ! cropland/grassland mosaic
  &           3.50,  &      ! cropland/woodland  mosaic
  &           2.00,  &      ! grassland
  &           2.50,  &      ! shrubland
  &           2.25,  &      ! mixed shrubland/grassland
  &           2.00,  &      ! savanna
  &           6.00,  &      ! decidous broadleaf forest
  &           5.00,  &      ! decidous needleleaf forest
  &           5.00,  &      ! evergreen broadleaf forest
  &           5.00,  &      ! evergreen needleleaf forest
  &           5.00,  &      ! mixed forest
  &           0.00,  &      ! water bodies
  &           2.00,  &      ! herbaceous wetland
  &           2.50,  &      ! wooded     wetland
  &           0.60,  &      ! barren or sparsely vegetated
  &           1.00,  &      ! herbaceous  tundra
  &           1.60,  &      ! wooded      tundra
  &           1.30,  &      ! mixed       tundra
  &           0.00,  &      ! bare ground tundra
  &           0.00  /)      ! snow or ice

REAL (KIND=wp) :: zrd_cosmo_o(nclass_glcc)  = (/ &         !< lookup table landuse class to root depth [m]
  &          0.60,  &            ! class urban and built-up land
  &          1.00,  &            ! dryland cropland and pasture
  &          0.60,  &            ! irrigated cropland and pasture
  &          0.80,  &            ! mixed dryland/irrigated .....
  &          1.00,  &            ! cropland/grassland mosaic
  &          1.00,  &            ! cropland/woodland  mosaic
  &          0.60,  &            ! grassland
  &          1.00,  &            ! shrubland
  &          1.00,  &            ! mixed shrubland/grassland
  &          2.00,  &            ! savanna
  &          1.00,  &            ! decidous broadleaf forest
  &          0.60,  &            ! decidous needleleaf forest
  &          1.00,  &            ! evergreen broadleaf forest
  &          0.60,  &            ! evergreen needleleaf forest
  &          0.80,  &            ! mixed forest
  &          0.00,  &            ! water bodies
  &          0.40,  &            ! herbaceous wetland
  &          0.40,  &            ! wooded     wetland
  &          0.30,  &            ! barren or sparsely vegetated
  &          0.10,  &            ! herbaceous  tundra
  &          0.10,  &            ! wooded      tundra
  &          0.10,  &            ! mixed       tundra
  &          0.00,  &            ! bare ground tundra
  &          0.00  /)            ! snow or ice




REAL (KIND=wp) :: zemiss_cosmo_o(nclass_glcc) = (/ &       !< lookup table landuse class to surface thermal emissivity
  &          0.960, &            ! class urban and built-up land
  &          0.985, &            ! dryland cropland and pasture
  &          0.990, &            ! irrigated cropland and pasture
  &          0.990, &            ! mixed dryland/irrigated .....
  &          0.990, &            ! cropland/grassland mosaic
  &          0.990, &            ! cropland/woodland  mosaic
  &          0.993, &            ! grassland
  &          0.985, &            ! shrubland
  &          0.985, &            ! mixed shrubland/grassland
  &          0.993, &            ! savanna
  &          0.990, &            ! decidous broadleaf forest
  &          0.990, &            ! decidous needleleaf forest
  &          0.996, &            ! evergreen broadleaf forest
  &          0.996, &            ! evergreen needleleaf forest
  &          0.993, &            ! mixed forest
  &          0.991, &            ! water bodies
  &          0.992, &            ! herbaceous wetland
  &          0.992, &            ! wooded     wetland
  &          0.950, &            ! barren or sparsely vegetated
  &          0.990, &            ! herbaceous  tundra
  &          0.990, &            ! wooded      tundra
  &          0.990, &            ! mixed       tundra
  &          0.990, &            ! bare ground tundra
  &          0.9999 /)            ! snow or ice


REAL (KIND=wp) :: zrs_min_cosmo_o(nclass_glcc) = (/ &      !< lookup table landuse class to minimal stomata resistance
  &          150.0, &       ! class urban and built-up land
  &          180.0, &       ! dryland cropland and pasture
  &          180.0, &       ! irrigated cropland and pasture
  &          180.0, &       ! mixed dryland/irrigated .....
  &          180.0, &       ! cropland/grassland mosaic
  &          200.0, &       ! cropland/woodland  mosaic
  &          110.0, &       ! grassland
  &          225.0, &       ! shrubland
  &          150.0, &       ! mixed shrubland/grassland
  &          100.0, &       ! savanna
  &          240.0, &       ! decidous broadleaf forest
  &          500.0, &       ! decidous needleleaf forest
  &          240.0, &       ! evergreen broadleaf forest
  &          500.0, &       ! evergreen needleleaf forest
  &          350.0, &       ! mixed forest
  &           0.00, &       ! water bodies
  &          110.0, &       ! herbaceous wetland
  &          110.0, &       ! wooded     wetland
  &          350.0, &       ! barren or sparsely vegetated
  &           80.0, &       ! herbaceous  tundra
  &           80.0, &       ! wooded      tundra
  &           80.0, &       ! mixed       tundra
  &           0.00,  &      ! bare ground tundra
  &           0.00  /)      ! snow or ice



!---------------------------------------------------------------------------------------------- 
!----------------------------------------------------------------------------------------------
REAL (KIND=wp) :: z0c_experimental(nclass_glcc)  = (/ &       !< lookup table landuse class to roughness length [m]
  &         1.00,       &       ! class urban and built-up land
  &          0.10,       &       ! dryland cropland and pasture
  &          0.10,       &       ! irrigated cropland and pasture
  &          0.10,       &       ! mixed dryland/irrigated .....
  &          0.07,       &       ! cropland/grassland mosaic
  &          0.25,       &       ! cropland/woodland  mosaic
  &          0.03,       &       ! grassland
  &          0.20,       &       ! shrubland
  &          0.15,       &       ! mixed shrubland/grassland
  &          0.15,       &       ! savanna
  &          1.00,       &       ! decidous broadleaf forest
  &          1.00,       &       ! decidous needleleaf forest
  &          1.00,       &       ! evergreen broadleaf forest
  &          1.00,       &       ! evergreen needleleaf forest
  &          1.00,       &       ! mixed forest
  &          0.0002,     &       ! water bodies
  &          0.05,       &       ! herbaceous wetland
  &          0.20,       &       ! wooded     wetland
  &          0.05,       &       ! barren or sparsely vegetated
  &          0.05,       &       ! herbaceous  tundra
  &          0.20,       &       ! wooded      tundra
  &          0.10,       &       ! mixed       tundra
  &          0.03,       &       ! bare ground tundra
  &          0.01   /)           ! snow or ice

REAL (KIND=wp) :: lnz0c_experimental(nclass_glcc)    !< corresponding natural logarithm of z0c_experimental_o


REAL (KIND=wp) :: zplcmnc_experimental(nclass_glcc) = (/ &      !< lookup table landuse class to minimal plant cover
  &          0.05,      &           ! class urban and built-up land
  &          0.45,      &           ! dryland cropland and pasture
  &          0.50,      &           ! irrigated cropland and pasture
  &          0.45,      &           ! mixed dryland/irrigated .....
  &          0.45,      &           ! cropland/grassland mosaic
  &          0.45,      &           ! cropland/woodland  mosaic
  &          1.00,      &           ! grassland
  &          0.10,      &           ! shrubland
  &          0.10,      &           ! mixed shrubland/grassland
  &          0.20,      &           ! savanna
  &          0.00,      &           ! decidous broadleaf forest
  &          0.00,      &           ! decidous needleleaf forest
  &          1.00,      &           ! evergreen broadleaf forest
  &          1.00,      &           ! evergreen needleleaf forest
  &          0.50,      &           ! mixed forest
  &          0.00,      &           ! water bodies
  &          0.40,      &           ! herbaceous wetland
  &          0.10,      &           ! wooded     wetland
  &          0.02,      &           ! barren or sparsely vegetated
  &          0.00,      &           ! herbaceous  tundra
  &          0.20,      &           ! wooded      tundra
  &          0.10,      &           ! mixed       tundra
  &          0.00,      &           ! bare ground tundra
  &          0.00   /)              ! snow or ice



REAL (KIND=wp) :: zplcmxc_experimental(nclass_glcc) = (/ &     !< lookup table landuse class to maximal plant cover
  &           0.20, &       ! class urban and built-up land
  &           0.70, &       ! dryland cropland and pasture
  &           0.90, &       ! irrigated cropland and pasture
  &           0.80, &       ! mixed dryland/irrigated .....
  &           0.90, &       ! cropland/grassland mosaic
  &           0.80, &       ! cropland/woodland  mosaic
  &           0.90, &       ! grassland
  &           0.80, &       ! shrubland
  &           0.80, &       ! mixed shrubland/grassland
  &           0.80, &       ! savanna
  &           0.90, &       ! decidous broadleaf forest
  &           0.90, &       ! decidous needleleaf forest
  &           0.80, &       ! evergreen broadleaf forest
  &           0.80, &       ! evergreen needleleaf forest
  &           0.90, &       ! mixed forest
  &           0.00, &       ! water bodies
  &           0.80, &       ! herbaceous wetland
  &           0.80, &       ! wooded     wetland
  &           0.05, &       ! barren or sparsely vegetated
  &           0.50, &       ! herbaceous  tundra
  &           0.50, &       ! wooded      tundra
  &           0.50, &       ! mixed       tundra
  &           0.00, &       ! bare ground tundra
  &           0.00  /)      ! snow or ice



REAL (KIND=wp) :: zlaimnc_experimental(nclass_glcc) = (/ &      !< lookup table landuse class to minimal leaf area index
  &           0.10,  &               ! class urban and built-up land
  &           0.20,  &               ! dryland cropland and pasture
  &           0.20,  &               ! irrigated cropland and pasture
  &           0.20,  &               ! mixed dryland/irrigated .....
  &           0.35,  &               ! cropland/grassland mosaic
  &           1.20,  &               ! cropland/woodland  mosaic
  &           0.50,  &               ! grassland
  &           0.10,  &               ! shrubland
  &           0.10,  &               ! mixed shrubland/grassland
  &           1.00,  &               ! savanna
  &           0.00,  &               ! decidous broadleaf forest
  &           0.00,  &               ! decidous needleleaf forest
  &           9.00,  &               ! evergreen broadleaf forest
  &           8.00,  &               ! evergreen needleleaf forest
  &           2.25,  &               ! mixed forest
  &           0.00,  &               ! water bodies
  &           1.00,  &               ! herbaceous wetland
  &           1.00,  &               ! wooded     wetland
  &           0.50,  &               ! barren or sparsely vegetated
  &           0.00,  &               ! herbaceous  tundra
  &           0.50,  &               ! wooded      tundra
  &           0.25,  &               ! mixed       tundra
  &           0.00,  &               ! bare ground tundra
  &           0.00   /)              ! snow or ice


REAL (KIND=wp) :: zlaimxc_experimental(nclass_glcc) = (/ &      !< lookup table landuse class to maximal leaf area index
  &           1.00,  &      ! class urban and built-up land
  &           3.30,  &      ! dryland cropland and pasture
  &           3.00,  &      ! irrigated cropland and pasture
  &           3.10,  &      ! mixed dryland/irrigated .....
  &           3.00,  &      ! cropland/grassland mosaic
  &           3.50,  &      ! cropland/woodland  mosaic
  &           2.00,  &      ! grassland
  &           2.50,  &      ! shrubland
  &           2.25,  &      ! mixed shrubland/grassland
  &           2.00,  &      ! savanna
  &           6.00,  &      ! decidous broadleaf forest
  &           5.00,  &      ! decidous needleleaf forest
  &           5.00,  &      ! evergreen broadleaf forest
  &           5.00,  &      ! evergreen needleleaf forest
  &           5.00,  &      ! mixed forest
  &           0.00,  &      ! water bodies
  &           2.00,  &      ! herbaceous wetland
  &           2.50,  &      ! wooded     wetland
  &           0.60,  &      ! barren or sparsely vegetated
  &           1.00,  &      ! herbaceous  tundra
  &           1.60,  &      ! wooded      tundra
  &           1.30,  &      ! mixed       tundra
  &           0.00,  &      ! bare ground tundra
  &           0.00  /)      ! snow or ice

REAL (KIND=wp) :: zrd_experimental(nclass_glcc)  = (/ &         !< lookup table landuse class to root depth [m]
  &          0.60,  &            ! class urban and built-up land
  &          1.00,  &            ! dryland cropland and pasture
  &          0.60,  &            ! irrigated cropland and pasture
  &          0.80,  &            ! mixed dryland/irrigated .....
  &          1.00,  &            ! cropland/grassland mosaic
  &          1.00,  &            ! cropland/woodland  mosaic
  &          0.60,  &            ! grassland
  &          1.00,  &            ! shrubland
  &          1.00,  &            ! mixed shrubland/grassland
  &          2.00,  &            ! savanna
  &          1.00,  &            ! decidous broadleaf forest
  &          0.60,  &            ! decidous needleleaf forest
  &          1.00,  &            ! evergreen broadleaf forest
  &          0.60,  &            ! evergreen needleleaf forest
  &          0.80,  &            ! mixed forest
  &          0.00,  &            ! water bodies
  &          0.40,  &            ! herbaceous wetland
  &          0.40,  &            ! wooded     wetland
  &          0.30,  &            ! barren or sparsely vegetated
  &          0.10,  &            ! herbaceous  tundra
  &          0.10,  &            ! wooded      tundra
  &          0.10,  &            ! mixed       tundra
  &          0.00,  &            ! bare ground tundra
  &          0.00  /)            ! snow or ice




REAL (KIND=wp) :: zemiss_experimental(nclass_glcc) = (/ &       !< lookup table landuse class to surface thermal emissivity
  &          0.960, &            ! class urban and built-up land
  &          0.985, &            ! dryland cropland and pasture
  &          0.990, &            ! irrigated cropland and pasture
  &          0.990, &            ! mixed dryland/irrigated .....
  &          0.990, &            ! cropland/grassland mosaic
  &          0.990, &            ! cropland/woodland  mosaic
  &          0.993, &            ! grassland
  &          0.985, &            ! shrubland
  &          0.985, &            ! mixed shrubland/grassland
  &          0.993, &            ! savanna
  &          0.990, &            ! decidous broadleaf forest
  &          0.990, &            ! decidous needleleaf forest
  &          0.996, &            ! evergreen broadleaf forest
  &          0.996, &            ! evergreen needleleaf forest
  &          0.993, &            ! mixed forest
  &          0.991, &            ! water bodies
  &          0.992, &            ! herbaceous wetland
  &          0.992, &            ! wooded     wetland
  &          0.950, &            ! barren or sparsely vegetated
  &          0.990, &            ! herbaceous  tundra
  &          0.990, &            ! wooded      tundra
  &          0.990, &            ! mixed       tundra
  &          0.990, &            ! bare ground tundra
  &          0.9999 /)            ! snow or ice


REAL (KIND=wp) :: zrs_min_experimental(nclass_glcc) = (/ &      !< lookup table landuse class to minimal stomata resistance
  &          120.0, &       ! class urban and built-up land
  &          120.0, &       ! dryland cropland and pasture
  &          120.0, &       ! irrigated cropland and pasture
  &          120.0, &       ! mixed dryland/irrigated .....
  &          100.0, &       ! cropland/grassland mosaic
  &          120.0, &       ! cropland/woodland  mosaic
  &          40.0, &       ! grassland
  &          40.0, &       ! shrubland
  &          40.0, &       ! mixed shrubland/grassland
  &          40.0, &       ! savanna
  &          150.0, &       ! decidous broadleaf forest
  &          150.0, &       ! decidous needleleaf forest
  &          250.0, &       ! evergreen broadleaf forest
  &          150.0, &       ! evergreen needleleaf forest
  &          150.0, &       ! mixed forest
  &          120.0, &       ! water bodies
  &          40.0, &       ! herbaceous wetland
  &          120.0, &       ! wooded     wetland
  &          120.0, &       ! barren or sparsely vegetated
  &           40.0, &       ! herbaceous  tundra
  &           40.0, &       ! wooded      tundra
  &           40.0, &       ! mixed       tundra
  &          120.0,  &      ! bare ground tundra
  &          120.0  /)      ! snow or ice

!---------------------------------------------------------------------------------------------- 
!----------------------------------------------------------------------------------------------





!> legend of the GLCC vegetation classes
CHARACTER(len=29) :: glcc_legend(nclass_glcc) = (/&   
 &  'class urban and built-up land' , &  
 &  'dryland cropland and pasture ' , &  
 &  'irrigated cropland and pastur' , & 
 &  'mixed dryland/irrigated      ' , &   
 &  'cropland/grassland mosaic    ' , &  
 &  'cropland/woodland  mosaic    ' , & 
 &  'grassland                    ' , &
 &  'shrubland                    ' , &
 &  'mixed shrubland/grassland    ' , & 
 &  'savanna                      ' , &  
 &  'decidous broadleaf forest    ' , &  
 &  'decidous needleleaf forest   ' , &   
 &  'evergreen broadleaf forest   ' , &  
 &  'evergreen needleleaf forest  ' , &  
 &  'mixed forest                 ' , &  
 &  'water bodies                 ' , &  
 &  'herbaceous wetland           ' , &  
 &  'wooded     wetland           ' , &   
 &  'barren or sparsely vegetated ' , &   
 &  'herbaceous  tundra           ' , &  
 &  'wooded      tundra           ' , &   
 &  'mixed       tundra           ' , &  
 &  'bare ground tundra           ' , &
 &  'snow or ice                  ' /)  




CONTAINS

  !> define lookup table for GLCC landuse classes
  SUBROUTINE init_glcc_lookup_tables(nclass_glcc, &
    &      ilookup_table_glcc, &
    &      z0_lt_glcc,           &
    &      lnz0_lt_glcc,       &
    &      plc_mn_lt_glcc,        &
    &      plc_mx_lt_glcc,        &
    &      lai_mn_lt_glcc,        &
    &      lai_mx_lt_glcc,        &
    &      rd_lt_glcc,          &
    &      emiss_lt_glcc,       &
    &      rs_min_lt_glcc)
    INTEGER, INTENT(IN) :: nclass_glcc !< GLCC has 24 classes for the land use description
    INTEGER, INTENT(IN) :: ilookup_table_glcc  !< integer switch to choose a lookup table
    REAL (KIND=wp), INTENT(OUT) :: z0_lt_glcc(nclass_glcc)      !< lookup table landuse class to roughness length [m]
    REAL (KIND=wp), INTENT(OUT) :: lnz0_lt_glcc(nclass_glcc)    !< corresponding natural logarithm of z0c_gme_o
    REAL (KIND=wp), INTENT(OUT) :: plc_mn_lt_glcc(nclass_glcc)  !< lookup table landuse class to minimal plant cover
    REAL (KIND=wp), INTENT(OUT) :: plc_mx_lt_glcc(nclass_glcc)  !< lookup table landuse class to maximal plant cover
    REAL (KIND=wp), INTENT(OUT) :: lai_mn_lt_glcc(nclass_glcc)  !< lookup table landuse class to minimal leaf area index
    REAL (KIND=wp), INTENT(OUT) :: lai_mx_lt_glcc(nclass_glcc)  !< lookup table landuse class to maximal leaf area index
    REAL (KIND=wp), INTENT(OUT) :: rd_lt_glcc(nclass_glcc)      !< lookup table landuse class to root depth [m]
    REAL (KIND=wp), INTENT(OUT) :: emiss_lt_glcc(nclass_glcc)   !< lookup table landuse class to surface thermal emissivity
    REAL (KIND=wp), INTENT(OUT) :: rs_min_lt_glcc(nclass_glcc)  !< lookup table landuse class to minimal stomata resistance

    ! local variable
    INTEGER :: i !< counter
    REAL :: arg

      SELECT CASE (ilookup_table_glcc)
        CASE(i_gme_lookup_table_glcc)
           z0_lt_glcc = z0c_gme_o
           plc_mn_lt_glcc = zplcmnc_gme_o
           plc_mx_lt_glcc = zplcmxc_gme_o
           lai_mn_lt_glcc = zlaimnc_gme_o
           lai_mx_lt_glcc = zlaimxc_gme_o
           rd_lt_glcc = zrd_gme_o
           emiss_lt_glcc = zemiss_gme_o
           rs_min_lt_glcc = zrs_min_gme_o
        CASE(i_cosmo_lookup_table_glcc)
           z0_lt_glcc = z0c_cosmo_o
           plc_mn_lt_glcc = zplcmnc_cosmo_o
           plc_mx_lt_glcc = zplcmxc_cosmo_o
           lai_mn_lt_glcc = zlaimnc_cosmo_o
           lai_mx_lt_glcc = zlaimxc_cosmo_o
           rd_lt_glcc = zrd_cosmo_o
           emiss_lt_glcc = zemiss_cosmo_o
           rs_min_lt_glcc = zrs_min_cosmo_o
        CASE(i_experimental_lookup_table_glcc)
           z0_lt_glcc = z0c_experimental
           plc_mn_lt_glcc = zplcmnc_experimental
           plc_mx_lt_glcc = zplcmxc_experimental
           lai_mn_lt_glcc = zlaimnc_experimental
           lai_mx_lt_glcc = zlaimxc_experimental
           rd_lt_glcc = zrd_experimental
           emiss_lt_glcc = zemiss_experimental
           rs_min_lt_glcc = zrs_min_experimental
        CASE DEFAULT
           z0_lt_glcc = z0c_gme_o
           plc_mn_lt_glcc = zplcmnc_gme_o
           plc_mx_lt_glcc = zplcmxc_gme_o
           lai_mn_lt_glcc = zlaimnc_gme_o
           lai_mx_lt_glcc = zlaimxc_gme_o
           rd_lt_glcc = zrd_gme_o
           emiss_lt_glcc = zemiss_gme_o
           rs_min_lt_glcc = zrs_min_gme_o
      END SELECT

      DO i=1,nclass_glcc
        IF (z0_lt_glcc(i) > 0.) THEN
          arg = z0_lt_glcc(i)
          lnz0_lt_glcc(i) = ALOG(arg)
        ENDIF
      ENDDO

  END  SUBROUTINE init_glcc_lookup_tables

    !> define  name of lookup table for GLCc 
  SUBROUTINE get_name_glcc_lookup_tables(ilookup_table_glcc, name_lookup_table_glcc)
    INTEGER, INTENT(IN) :: ilookup_table_glcc  !< integer switch to choose a lookup table
    CHARACTER (LEN=filename_max), INTENT(OUT) :: name_lookup_table_glcc !< name of lookup table
    ! local variable
      SELECT CASE (ilookup_table_glcc)
        CASE(i_gme_lookup_table_glcc)
           name_lookup_table_glcc='Ritter_2007'
        CASE(i_cosmo_lookup_table_glcc)
           name_lookup_table_glcc='Heise_2005'
        CASE(i_experimental_lookup_table_glcc)
           name_lookup_table_glcc='Asensio_2010'
        CASE DEFAULT
           name_lookup_table_glcc='Ritter_2007'
      END SELECT

  END  SUBROUTINE get_name_glcc_lookup_tables



   !> assign the GLCC land use classes to some characteristic (more or less) physical parameters
  SUBROUTINE glcc_look_up(lu, &
    &      nclass_glcc, &
    &      lnz0_lt_glcc,          &
    &      plc_mn_lt_glcc,        &
    &      plc_mx_lt_glcc,        &
    &      lai_mn_lt_glcc,        &
    &      lai_mx_lt_glcc,        &
    &      rd_lt_glcc,            &
    &      emiss_lt_glcc,         &
    &      rs_min_lt_glcc,        &
    &      pland,          &
    &      pice,           &
    &      plnz0,          &
    &      proot,          &
    &      pmn,            &
    &      pmx,            &
    &      plaimn,         &
    &      plaimx,         &
    &      purb,           &
    &      pfor_d,         &
    &      pfor_e,         &
    &      pemissivity,    &
    &      prs_min,        &
    &      k_error)
  
  INTEGER, INTENT(IN) :: lu             !< land use class
  INTEGER, INTENT(IN) :: nclass_glcc !< GLCC has 24 classes for the land use description
  REAL (KIND=wp), INTENT(IN) :: lnz0_lt_glcc(nclass_glcc)    !< corresponding natural logarithm of z0c_gme_o
  REAL (KIND=wp), INTENT(IN) :: plc_mn_lt_glcc(nclass_glcc)  !< lookup table landuse class to minimal plant cover
  REAL (KIND=wp), INTENT(IN) :: plc_mx_lt_glcc(nclass_glcc)  !< lookup table landuse class to maximal plant cover
  REAL (KIND=wp), INTENT(IN) :: lai_mn_lt_glcc(nclass_glcc)  !< lookup table landuse class to minimal leaf area index
  REAL (KIND=wp), INTENT(IN) :: lai_mx_lt_glcc(nclass_glcc)  !< lookup table landuse class to maximal leaf area index
  REAL (KIND=wp), INTENT(IN) :: rd_lt_glcc(nclass_glcc)      !< lookup table landuse class to root depth [m]
  REAL (KIND=wp), INTENT(IN) :: emiss_lt_glcc(nclass_glcc)   !< lookup table landuse class to surface thermal emissivity
  REAL (KIND=wp), INTENT(IN) :: rs_min_lt_glcc(nclass_glcc)  !< lookup table landuse class to minimal stomata resistance

  REAL (KIND=wp), INTENT(OUT) :: pland          !< land cover                      (-)
  REAL (KIND=wp), INTENT(OUT) :: pice           !< ice fraction                    (-)
  REAL (KIND=wp), INTENT(OUT) :: plnz0          !< logarithm of roughness length   (m)
  REAL (KIND=wp), INTENT(OUT) :: proot          !< root depth                      (m)
  REAL (KIND=wp), INTENT(OUT) :: pmn            !< minimal plant cover             (-)
  REAL (KIND=wp), INTENT(OUT) :: pmx            !< maximum plant cover             (-)
  REAL (KIND=wp), INTENT(OUT) :: plaimn         !< minimal leaf area index         (m**2/m**2)
  REAL (KIND=wp), INTENT(OUT) :: plaimx         !< maximum leaf area index         (m**2/m**2)
  REAL (KIND=wp), INTENT(OUT) :: purb           !< urbanisation                    (-)
  REAL (KIND=wp), INTENT(OUT) :: pfor_d         !< deciduous forest                (-)
  REAL (KIND=wp), INTENT(OUT) :: pfor_e         !< evergreen forest                (-)
  REAL (KIND=wp), INTENT(OUT) :: pemissivity    !< surface thermal emissivity      (-)
  REAL (KIND=wp), INTENT(OUT) :: prs_min        !< minimum stomata resistance      (s/m)

  INTEGER, INTENT(OUT)        :: k_error     !< error return code
  
  ! local variables

       ! Test for true land points                     
          IF (lu>=1 .AND. lu<=24 .AND.lu/=16) THEN
            k_error     = 0
            pland       = 1.0
            plnz0       = lnz0_lt_glcc(lu)
            pmn         = plc_mn_lt_glcc(lu)
            pmx         = plc_mx_lt_glcc(lu)
            plaimn      = lai_mn_lt_glcc(lu)
            plaimx      = lai_mx_lt_glcc(lu)
            proot       = rd_lt_glcc(lu)
            prs_min     = rs_min_lt_glcc(lu)     
            pemissivity = emiss_lt_glcc(lu)
            purb    = 0.0
            pfor_d  = 0.0
            pfor_e  = 0.0
            pice    = 0.0

            IF (lu== 1            )   purb   = 1.0  ! urban pixel
            IF (lu== 11 .OR. lu== 12) pfor_d = 1.0  ! dec.forest 
            IF (lu== 13 .OR. lu== 14) pfor_e = 1.0  ! eve.forest 
            IF (lu== 15) THEN                       ! mix.forest
              pfor_d = 0.5  
              pfor_e = 0.5  
            END IF
            IF (lu==24            ) pice   = 1.0  ! ice or snow pixel
          ELSE IF (lu==16) THEN ! water
            k_error     = 0
            pland       = 0.0
            pemissivity = emiss_lt_glcc(lu)             ! emissivity is required everywhere
          ELSE
            k_error     = 1  ! not a valid land use class
            pland       = 0.0 
          END IF

  END  SUBROUTINE glcc_look_up






END MODULE mo_glcc_lookup_tables
