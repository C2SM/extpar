!+ Fortran Module with lookup-tables for the GLC2000 data
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V1_6         2011/11/29 Jan-Peter Schulz
!  Correction of values for plant cover (min and max) for land use
!  class no. 19 (bare areas) in the Heise (2005) look-up table.
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran Module with lookup-tables for the GLC2000 data
!> \author Hermann Asensio
!!
!! Description:
!! The GLC2000 dataset contains the following land use classification scheme
!!
!! class no.          description               !        IGBP-correspondence
!!
!!  01   'evergreen broadleaf tree     ' , &    !        evergreen broadleaf forest
!!  02   'deciduous broadleaf tree clos' , &    !        deciduous broadleaf forest
!!  03   'deciduous broadleaf tree open' , &    !        woody savannas & savanna
!!  04   'evergreen needleleaf tree    ' , &    !        evergreen needleleaf forest
!!  05   'deciduous needleleaf tree    ' , &    !        deciduous needleleaf forest
!!  06   'mixed leaf tree              ' , &    !        mixed forest
!!  07   'fresh water flooded tree     ' , &    !        evergreen broadleaf
!!  08   'saline water flooded tree    ' , &    !        evergreen broadleaf
!!  09   'mosaic tree / other nat.veg. ' , &    !        shrubland ?
!!  10   'burnt tree cover             ' , &    ! 
!!  11   'evergreen shrubs closed-open ' , &    !        shrubland open-closed
!!  12   'deciduous shrubs closed-open ' , &    !        savannas
!!  13   'herbaceous cover closed-open ' , &    !        grassland
!!  14   'sparse herbaceous or grass   ' , &    !
!!  15   'flooded shrub or herbaceous  ' , &    !        wetlands
!!  16   'cultivated & managed areas   ' , &    !        croplands
!!  17   'mosaic crop/tree/natural veg.' , &    !        cropland/other vegetation mosaic
!!  18   'mosaic crop/shrub or grass   ' , &    !        cropland/other vegetation mosaic
!!  19   'bare areas                   ' , &    !        barren or sparsely vegetated
!!  20   'water bodies                 ' , &    !        water
!!  21   'snow & ice                   ' , &    !        snow & ice
!!  22   'artificial surfaces          ' , &    !        urban and built-up areas
!!  23   'undefined                    ' /      ! 
!!
!! Different look-up tables are given in this module
!! - operational settings of GME (Ritter, 2007)
!! - operational settings of COSMO (Heise, 2005)
!! - experimantal settings, lookup-table analog to Ecoclimap (Masson 2003), use only in combination with NDVI data!
MODULE mo_glc2000_lookup_tables

 !> kind parameters are defined in MODULE data_parameters
  USE mo_kind, ONLY: wp
  USE mo_kind, ONLY: i8
  USE mo_kind, ONLY: i4

  !> abort_extpar defined in MODULE utilities_extpar
  USE mo_utilities_extpar, ONLY: abort_extpar
  USE mo_io_units,          ONLY: filename_max



                           
IMPLICIT NONE

PRIVATE

PUBLIC :: init_glc2000_lookup_tables
PUBLIC :: get_name_glc2000_lookup_tables
PUBLIC :: glc2000_look_up
PUBLIC :: glc2000_legend

PUBLIC :: nclass_glc2000
PUBLIC :: ilookup_table_glc2000
PUBLIC :: name_lookup_table_glc2000
PUBLIC :: i_gme_lookup_table, i_cosmo_lookup_table, i_experimental_lookup_table
PUBLIC :: z0_lt_glc2000, lnz0_lt_glc2000, plc_mn_lt_glc2000, plc_mx_lt_glc2000 
PUBLIC :: lai_mn_lt_glc2000, lai_mx_lt_glc2000, rd_lt_glc2000, emiss_lt_glc2000, rs_min_lt_glc2000         



INTEGER (KIND=i4), PARAMETER :: nclass_glc2000 = 23 !< GLC2000 has 23 classes for the land use description

INTEGER (KIND=i4), PARAMETER :: i_gme_lookup_table = 1 !< lookup_table for GLC2000 land use classes (IGBP correspondence) from
                                                       !< the operational settings of GME (Ritter, 2007)
INTEGER (KIND=i4), PARAMETER :: i_cosmo_lookup_table = 2 !< lookup_table for GLC2000 land use classes (IGBP correspondence) from
                                                         !< the operational settings of COSMO (Heise, 2005)
INTEGER (KIND=i4), PARAMETER :: i_experimental_lookup_table = 3 !< lookup_table for GLC2000 land use classes
                                                                !< (IGBP correspondence) for experimental setting, analogue to
                                                                !< look-up tables of ECOCLIMAP (Masson 2003)

INTEGER (KIND=i4) :: ilookup_table_glc2000 !< integer switch to choose a lookup table
CHARACTER (LEN=filename_max) :: name_lookup_table_glc2000 !< name of lookup table




REAL (KIND=wp) :: z0_lt_glc2000(nclass_glc2000)      !< lookup table landuse class to roughness length [m]
REAL (KIND=wp) :: lnz0_lt_glc2000(nclass_glc2000)    !< corresponding natural logarithm of z0c_gme_o
REAL (KIND=wp) :: plc_mn_lt_glc2000(nclass_glc2000)  !< lookup table landuse class to minimal plant cover
REAL (KIND=wp) :: plc_mx_lt_glc2000(nclass_glc2000)  !< lookup table landuse class to maximal plant cover
REAL (KIND=wp) :: lai_mn_lt_glc2000(nclass_glc2000)  !< lookup table landuse class to minimal leaf area index
REAL (KIND=wp) :: lai_mx_lt_glc2000(nclass_glc2000)  !< lookup table landuse class to maximal leaf area index
REAL (KIND=wp) :: rd_lt_glc2000(nclass_glc2000)      !< lookup table landuse class to root depth [m]
REAL (KIND=wp) :: emiss_lt_glc2000(nclass_glc2000)   !< lookup table landuse class to surface thermal emissivity
REAL (KIND=wp) :: rs_min_lt_glc2000(nclass_glc2000)  !< lookup table landuse class to minimal stomata resistance

!---------------------------------------------------------------------------------------------- 
!----------------------------------------------------------------------------------------------
REAL (KIND=wp) :: z0c_gme_o(nclass_glc2000)  = (/ &       !< lookup table landuse class to roughness length [m]
  &          1.00,       &       ! evergreen broadleaf forest   
  &          1.00,       &       ! deciduous broadleaf closed forest
  &          0.15,       &       ! deciduous broadleaf open   forest
  &          1.00,       &       ! evergreen needleleaf forest   
  &          1.00,       &       ! deciduous needleleaf forest
  &          1.00,       &       ! mixed leaf trees            
  &          1.00,       &       ! fresh water flooded trees
  &          1.00,       &       ! saline water flooded trees
  &          0.20,       &       ! mosaic tree / natural vegetation
  &          0.05,       &       ! burnt tree cover
  &          0.20,       &       ! evergreen shrubs closed-open
  &          0.15,       &       ! decidous shrubs closed-open
  &          0.03,       &       ! herbaceous vegetation closed-open
  &          0.05,       &       ! sparse herbaceous or grass 
  &          0.05,       &       ! flooded shrubs or herbaceous
  &          0.07,       &       ! cultivated & managed areas
  &          0.25,       &       ! mosaic crop / tree / natural vegetation
  &          0.07,       &       ! mosaic crop / shrub / grass
  &          0.05,       &       ! bare areas                       
  &          0.0002,     &       ! water
  &          0.01,       &       ! snow & ice 
  &          1.00,       &       ! artificial surface  
  &          0.         /)        ! undefined

REAL (KIND=wp) :: zplcmnc_gme_o(nclass_glc2000) = (/ &      !< lookup table landuse class to minimal plant cover
  &         0.80,  &       ! evergreen broadleaf forest   
  &         0.75,  &       ! deciduous broadleaf closed forest
  &         0.70,  &       ! deciduous broadleaf open   forest
  &         0.80,  &       ! evergreen needleleaf forest   
  &         0.75,  &       ! deciduous needleleaf forest
  &         0.75,  &       ! mixed leaf trees            
  &         0.80,  &       ! fresh water flooded trees
  &         0.80,  &       ! saline water flooded trees
  &         0.70,  &       ! mosaic tree / natural vegetation
  &         0.50,  &       ! burnt tree cover
  &         0.70,  &       ! evergreen shrubs closed-open
  &         0.70,  &       ! decidous shrubs closed-open
  &         0.75,  &       ! herbaceous vegetation closed-open
  &         0.50,  &       ! sparse herbaceous or grass 
  &         0.70,  &       ! flooded shrubs or herbaceous
  &         0.50,  &       ! cultivated & managed areas
  &         0.65,  &       ! mosaic crop / tree / natural vegetation
  &         0.50,  &       ! mosaic crop / shrub / grass
  &         0.02,  &       ! bare areas                         
  &         0.00,  &       ! water
  &         0.00,  &       ! snow & ice 
  &         0.10,  &       ! artificial surface  
  &         0.       /)    ! undefined


REAL (KIND=wp) :: zplcmxc_gme_o(nclass_glc2000) = (/ &     !< lookup table landuse class to maximal plant cover
 &            0.80, &       ! evergreen broadleaf forest   
 &            0.90, &       ! deciduous broadleaf closed forest
 &            0.80, &       ! deciduous broadleaf open   forest
 &            0.80, &       ! evergreen needleleaf forest   
 &            0.90, &       ! deciduous needleleaf forest
 &            0.90, &       ! mixed leaf trees            
 &            0.80, &       ! fresh water flooded trees
 &            0.80, &       ! saline water flooded trees
 &            0.80, &       ! mosaic tree / natural vegetation
 &            0.50, &       ! burnt tree cover
 &            0.80, &       ! evergreen shrubs closed-open
 &            0.80, &       ! decidous shrubs closed-open
 &            0.90, &       ! herbaceous vegetation closed-open
 &            0.50, &       ! sparse herbaceous or grass 
 &            0.80, &       ! flooded shrubs or herbaceous
 &            0.90, &       ! cultivated & managed areas
 &            0.80, &       ! mosaic crop / tree / natural vegetation
 &            0.90, &       ! mosaic crop / shrub / grass
 &            0.05, &       ! bare areas                       
 &            0.00, &       ! water
 &            0.00, &       ! snow & ice 
 &            0.20, &       ! artificial surface   
 &            0.       /)   ! undefined

REAL (KIND=wp) :: zlaimnc_gme_o(nclass_glc2000) = (/ &      !< lookup table landuse class to minimal leaf area index
 &          1.40, &        ! evergreen broadleaf forest   
 &          1.00, &        ! deciduous broadleaf closed forest
 &          1.00, &        ! deciduous broadleaf open   forest
 &          1.30, &        ! evergreen needleleaf forest   
 &          1.00, &        ! deciduous needleleaf forest
 &          1.10, &        ! mixed leaf trees            
 &          1.40, &        ! fresh water flooded trees
 &          1.40, &        ! saline water flooded trees
 &          0.60, &        ! mosaic tree / natural vegetation
 &          0.40, &        ! burnt tree cover
 &          0.60, &        ! evergreen shrubs closed-open
 &          1.00, &        ! decidous shrubs closed-open
 &          1.00, &        ! herbaceous vegetation closed-open
 &          0.40, &        ! sparse herbaceous or grass 
 &          1.00, &        ! flooded shrubs or herbaceous
 &          0.70, &        ! cultivated & managed areas
 &          1.00, &        ! mosaic crop / tree / natural vegetation
 &          0.70, &        ! mosaic crop / shrub / grass
 &          0.40, &        ! bare areas                       
 &          0.00, &        ! water
 &          0.00, &        ! snow & ice 
 &          0.10, &        ! artificial surface   
 &          0.       /)    ! undefined


REAL (KIND=wp) :: zlaimxc_gme_o(nclass_glc2000) = (/ &      !< lookup table landuse class to maximal leaf area index
 &            5.00, &       ! evergreen broadleaf forest   
 &            6.00, &       ! deciduous broadleaf closed forest
 &            4.00, &       ! deciduous broadleaf open   forest
 &            5.00, &       ! evergreen needleleaf forest   
 &            5.00, &       ! deciduous needleleaf forest
 &            5.00, &       ! mixed leaf trees            
 &            5.00, &       ! fresh water flooded trees
 &            5.00, &       ! saline water flooded trees
 &            2.50, &       ! mosaic tree / natural vegetation
 &            0.60, &       ! burnt tree cover
 &            3.00, &       ! evergreen shrubs closed-open
 &            1.50, &       ! decidous shrubs closed-open
 &            3.10, &       ! herbaceous vegetation closed-open
 &            0.60, &       ! sparse herbaceous or grass 
 &            2.00, &       ! flooded shrubs or herbaceous
 &            3.30, &       ! cultivated & managed areas
 &            3.00, &       ! mosaic crop / tree / natural vegetation
 &            3.50, &       ! mosaic crop / shrub / grass
 &            0.60, &       ! bare areas                       
 &            0.00, &       ! water
 &            0.00, &       ! snow & ice 
 &            1.00, &       ! artificial surface   
 &            0.      /)    ! undefined

REAL (KIND=wp) :: zrd_gme_o(nclass_glc2000)  = (/ &         !< lookup table landuse class to root depth [m]
  &          1.00,       &       ! evergreen broadleaf forest   
  &          1.00,       &       ! deciduous broadleaf closed forest
  &          2.00,       &       ! deciduous broadleaf open   forest
  &          0.60,       &       ! evergreen needleleaf forest   
  &          0.60,       &       ! deciduous needleleaf forest
  &          0.80,       &       ! mixed leaf trees            
  &          1.00,       &       ! fresh water flooded trees
  &          1.00,       &       ! saline water flooded trees
  &          1.00,       &       ! mosaic tree / natural vegetation
  &          0.30,       &       ! burnt tree cover
  &          1.00,       &       ! evergreen shrubs closed-open
  &          2.00,       &       ! decidous shrubs closed-open
  &          0.60,       &       ! herbaceous vegetation closed-open
  &          0.30,       &       ! sparse herbaceous or grass 
  &          0.40,       &       ! flooded shrubs or herbaceous
  &          1.00,       &       ! cultivated & managed areas
  &          1.00,       &       ! mosaic crop / tree / natural vegetation
  &          1.00,       &       ! mosaic crop / shrub / grass
  &          0.30,       &       ! bare areas                       
  &          0.00,       &       ! water
  &          0.00,       &       ! snow & ice 
  &          0.60,       &       ! artificial surface   
  &          0.         /)       ! undefined


REAL (KIND=wp) :: zemiss_gme_o(nclass_glc2000) = (/ &       !< lookup table landuse class to surface thermal emissivity
  &          0.996,      &       ! evergreen broadleaf forest   
  &          0.990,      &       ! deciduous broadleaf closed forest
  &          0.993,      &       ! deciduous broadleaf open   forest
  &          0.996,      &       ! evergreen needleleaf forest   
  &          0.990,      &       ! deciduous needleleaf forest
  &          0.993,      &       ! mixed leaf trees            
  &          0.996,      &       ! fresh water flooded trees
  &          0.996,      &       ! saline water flooded trees
  &          0.985,      &       ! mosaic tree / natural vegetation
  &          0.950,      &       ! burnt tree cover
  &          0.985,      &       ! evergreen shrubs closed-open
  &          0.993,      &       ! decidous shrubs closed-open
  &          0.993,      &       ! herbaceous vegetation closed-open
  &          0.950,      &       ! sparse herbaceous or grass 
  &          0.992,      &       ! flooded shrubs or herbaceous
  &          0.990,      &       ! cultivated & managed areas
  &          0.990,      &       ! mosaic crop / tree / natural vegetation
  &          0.990,      &       ! mosaic crop / shrub / grass
  &          0.950,      &       ! bare areas                       
  &          0.991,      &       ! water
  &          0.9999,     &       ! snow & ice 
  &          0.960 ,     &       ! artificial surface   
  &          0.         /)       ! undefined

REAL (KIND=wp) :: zrs_min_gme_o(nclass_glc2000) = (/ &      !< lookup table landuse class to minimal stomata resistance
  &          175.0,      &       ! evergreen broadleaf forest   
  &          240.0,      &       ! deciduous broadleaf closed forest
  &          240.0,      &       ! deciduous broadleaf open   forest
  &          500.0,      &       ! evergreen needleleaf forest   
  &          500.0,      &       ! deciduous needleleaf forest
  &          350.0,      &       ! mixed leaf trees            
  &          350.0,      &       ! fresh water flooded trees
  &          350.0,      &       ! saline water flooded trees
  &          300.0,      &       ! mosaic tree / natural vegetation
  &          300.0,      &       ! burnt tree cover
  &          225.0,      &       ! evergreen shrubs closed-open
  &          225.0,      &       ! decidous shrubs closed-open
  &          150.0,      &       ! herbaceous vegetation closed-open
  &          110.0,      &       ! sparse herbaceous or grass 
  &          110.0,      &       ! flooded shrubs or herbaceous
  &          180.0,      &       ! cultivated & managed areas
  &          200.0,      &       ! mosaic crop / tree / natural vegetation
  &          150.0,      &       ! mosaic crop / shrub / grass
  &          150.0,      &       ! bare areas                       
  &          150.0,      &       ! water
  &          150.0,      &       ! snow & ice 
  &          150.0,      &       ! artificial surface   
  &          0.         /)       ! undefined

!---------------------------------------------------------------------------------------------- 
!----------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------- 
!----------------------------------------------------------------------------------------------


REAL (KIND=wp) :: z0c_cosmo_o(nclass_glc2000)  = (/ &       !< lookup table landuse class to roughness length [m]
  &          1.00,       &       ! evergreen broadleaf forest   
  &          1.00,       &       ! deciduous broadleaf closed forest
  &          0.15,       &       ! deciduous broadleaf open   forest
  &          1.00,       &       ! evergreen needleleaf forest   
  &          1.00,       &       ! deciduous needleleaf forest
  &          1.00,       &       ! mixed leaf trees            
  &          1.00,       &       ! fresh water flooded trees
  &          1.00,       &       ! saline water flooded trees
  &          0.20,       &       ! mosaic tree / natural vegetation
  &          0.05,       &       ! burnt tree cover
  &          0.20,       &       ! evergreen shrubs closed-open
  &          0.15,       &       ! decidous shrubs closed-open
  &          0.03,       &       ! herbaceous vegetation closed-open
  &          0.05,       &       ! sparse herbaceous or grass 
  &          0.05,       &       ! flooded shrubs or herbaceous
  &          0.07,       &       ! cultivated & managed areas
  &          0.25,       &       ! mosaic crop / tree / natural vegetation
  &          0.07,       &       ! mosaic crop / shrub / grass
  &          0.05,       &       ! bare areas                       
  &          0.0002,     &       ! water
  &          0.01,       &       ! snow & ice 
  &          1.00,       &       ! artificial surface  
  &          0.         /)        ! undefined

REAL (KIND=wp) :: zplcmnc_cosmo_o(nclass_glc2000)  = (/ &      !< lookup table landuse class to minimal plant cover
  &         0.80,  &       ! evergreen broadleaf forest   
  &         0.75,  &       ! deciduous broadleaf closed forest
  &         0.70,  &       ! deciduous broadleaf open   forest
  &         0.80,  &       ! evergreen needleleaf forest   
  &         0.75,  &       ! deciduous needleleaf forest
  &         0.75,  &       ! mixed leaf trees            
  &         0.80,  &       ! fresh water flooded trees
  &         0.80,  &       ! saline water flooded trees
  &         0.70,  &       ! mosaic tree / natural vegetation
  &         0.50,  &       ! burnt tree cover
  &         0.70,  &       ! evergreen shrubs closed-open
  &         0.70,  &       ! decidous shrubs closed-open
  &         0.75,  &       ! herbaceous vegetation closed-open
  &         0.50,  &       ! sparse herbaceous or grass 
  &         0.70,  &       ! flooded shrubs or herbaceous
  &         0.50,  &       ! cultivated & managed areas
  &         0.65,  &       ! mosaic crop / tree / natural vegetation
  &         0.50,  &       ! mosaic crop / shrub / grass
  &         0.02,  &       ! bare areas                         
  &         0.00,  &       ! water
  &         0.00,  &       ! snow & ice 
  &         0.10,  &       ! artificial surface  
  &         0.       /)    ! undefined

REAL (KIND=wp) :: zplcmxc_cosmo_o(nclass_glc2000)  = (/ &     !< lookup table landuse class to maximal plant cover
 &            0.80, &       ! evergreen broadleaf forest   
 &            0.90, &       ! deciduous broadleaf closed forest
 &            0.80, &       ! deciduous broadleaf open   forest
 &            0.80, &       ! evergreen needleleaf forest   
 &            0.90, &       ! deciduous needleleaf forest
 &            0.90, &       ! mixed leaf trees            
 &            0.80, &       ! fresh water flooded trees
 &            0.80, &       ! saline water flooded trees
 &            0.80, &       ! mosaic tree / natural vegetation
 &            0.50, &       ! burnt tree cover
 &            0.80, &       ! evergreen shrubs closed-open
 &            0.80, &       ! decidous shrubs closed-open
 &            0.90, &       ! herbaceous vegetation closed-open
 &            0.50, &       ! sparse herbaceous or grass 
 &            0.80, &       ! flooded shrubs or herbaceous
 &            0.90, &       ! cultivated & managed areas
 &            0.80, &       ! mosaic crop / tree / natural vegetation
 &            0.90, &       ! mosaic crop / shrub / grass
 &            0.05, &       ! bare areas                       
 &            0.00, &       ! water
 &            0.00, &       ! snow & ice 
 &            0.20, &       ! artificial surface   
 &            0.       /)   ! undefined
REAL (KIND=wp) :: zlaimnc_cosmo_o(nclass_glc2000)  = (/ &      !< lookup table landuse class to minimal leaf area index
 &          1.40, &        ! evergreen broadleaf forest   
 &          1.00, &        ! deciduous broadleaf closed forest
 &          1.00, &        ! deciduous broadleaf open   forest
 &          1.30, &        ! evergreen needleleaf forest   
 &          1.00, &        ! deciduous needleleaf forest
 &          1.10, &        ! mixed leaf trees            
 &          1.40, &        ! fresh water flooded trees
 &          1.40, &        ! saline water flooded trees
 &          0.60, &        ! mosaic tree / natural vegetation
 &          0.40, &        ! burnt tree cover
 &          0.60, &        ! evergreen shrubs closed-open
 &          1.00, &        ! decidous shrubs closed-open
 &          1.00, &        ! herbaceous vegetation closed-open
 &          0.40, &        ! sparse herbaceous or grass 
 &          1.00, &        ! flooded shrubs or herbaceous
 &          0.70, &        ! cultivated & managed areas
 &          1.00, &        ! mosaic crop / tree / natural vegetation
 &          0.70, &        ! mosaic crop / shrub / grass
 &          0.40, &        ! bare areas                       
 &          0.00, &        ! water
 &          0.00, &        ! snow & ice 
 &          0.10, &        ! artificial surface   
 &          0.       /)    ! undefined

REAL (KIND=wp) :: zlaimxc_cosmo_o(nclass_glc2000)  = (/ &      !< lookup table landuse class to maximal leaf area index
 &            2.40, &       ! evergreen broadleaf forest   
 &            3.40, &       ! deciduous broadleaf closed forest
 &            2.00, &       ! deciduous broadleaf open   forest
 &            3.80, &       ! evergreen needleleaf forest   
 &            3.80, &       ! deciduous needleleaf forest
 &            3.40, &       ! mixed leaf trees            
 &            2.40, &       ! fresh water flooded trees
 &            2.40, &       ! saline water flooded trees
 &            1.50, &       ! mosaic tree / natural vegetation
 &            0.60, &       ! burnt tree cover
 &            1.50, &       ! evergreen shrubs closed-open
 &            2.00, &       ! decidous shrubs closed-open
 &            3.10, &       ! herbaceous vegetation closed-open
 &            0.60, &       ! sparse herbaceous or grass 
 &            2.00, &       ! flooded shrubs or herbaceous
 &            3.30, &       ! cultivated & managed areas
 &            2.10, &       ! mosaic crop / tree / natural vegetation
 &            3.30, &       ! mosaic crop / shrub / grass
 &            0.60, &       ! bare areas                       
 &            0.00, &       ! water
 &            0.00, &       ! snow & ice 
 &            1.00, &       ! artificial surface   
 &            0.      /)    ! undefined
REAL (KIND=wp) :: zrd_cosmo_o(nclass_glc2000)    = (/ &         !< lookup table landuse class to root depth [m]
  &          1.00,       &       ! evergreen broadleaf forest   
  &          1.00,       &       ! deciduous broadleaf closed forest
  &          2.00,       &       ! deciduous broadleaf open   forest
  &          0.60,       &       ! evergreen needleleaf forest   
  &          0.60,       &       ! deciduous needleleaf forest
  &          0.80,       &       ! mixed leaf trees            
  &          1.00,       &       ! fresh water flooded trees
  &          1.00,       &       ! saline water flooded trees
  &          1.00,       &       ! mosaic tree / natural vegetation
  &          0.30,       &       ! burnt tree cover
  &          1.00,       &       ! evergreen shrubs closed-open
  &          2.00,       &       ! decidous shrubs closed-open
  &          0.60,       &       ! herbaceous vegetation closed-open
  &          0.30,       &       ! sparse herbaceous or grass 
  &          0.40,       &       ! flooded shrubs or herbaceous
  &          1.00,       &       ! cultivated & managed areas
  &          1.00,       &       ! mosaic crop / tree / natural vegetation
  &          1.00,       &       ! mosaic crop / shrub / grass
  &          0.30,       &       ! bare areas                       
  &          0.00,       &       ! water
  &          0.00,       &       ! snow & ice 
  &          0.60,       &       ! artificial surface   
  &          0.         /)       ! undefined
REAL (KIND=wp) :: zemiss_cosmo_o(nclass_glc2000)  = (/ &       !< lookup table landuse class to surface thermal emissivity
  &          0.996,      &       ! evergreen broadleaf forest   
  &          0.990,      &       ! deciduous broadleaf closed forest
  &          0.993,      &       ! deciduous broadleaf open   forest
  &          0.996,      &       ! evergreen needleleaf forest   
  &          0.990,      &       ! deciduous needleleaf forest
  &          0.993,      &       ! mixed leaf trees            
  &          0.996,      &       ! fresh water flooded trees
  &          0.996,      &       ! saline water flooded trees
  &          0.985,      &       ! mosaic tree / natural vegetation
  &          0.950,      &       ! burnt tree cover
  &          0.985,      &       ! evergreen shrubs closed-open
  &          0.993,      &       ! decidous shrubs closed-open
  &          0.993,      &       ! herbaceous vegetation closed-open
  &          0.950,      &       ! sparse herbaceous or grass 
  &          0.992,      &       ! flooded shrubs or herbaceous
  &          0.990,      &       ! cultivated & managed areas
  &          0.990,      &       ! mosaic crop / tree / natural vegetation
  &          0.990,      &       ! mosaic crop / shrub / grass
  &          0.950,      &       ! bare areas                       
  &          0.991,      &       ! water
  &          0.9999,     &       ! snow & ice 
  &          0.960 ,     &       ! artificial surface   
  &          0.         /)       ! undefined
REAL (KIND=wp) :: zrs_min_cosmo_o(nclass_glc2000)  = (/ &      !< lookup table landuse class to minimal stomata resistance
  &          175.0,      &       ! evergreen broadleaf forest   
  &          240.0,      &       ! deciduous broadleaf closed forest
  &          240.0,      &       ! deciduous broadleaf open   forest
  &          500.0,      &       ! evergreen needleleaf forest   
  &          500.0,      &       ! deciduous needleleaf forest
  &          350.0,      &       ! mixed leaf trees            
  &          350.0,      &       ! fresh water flooded trees
  &          350.0,      &       ! saline water flooded trees
  &          300.0,      &       ! mosaic tree / natural vegetation
  &          300.0,      &       ! burnt tree cover
  &          225.0,      &       ! evergreen shrubs closed-open
  &          225.0,      &       ! decidous shrubs closed-open
  &          150.0,      &       ! herbaceous vegetation closed-open
  &          110.0,      &       ! sparse herbaceous or grass 
  &          110.0,      &       ! flooded shrubs or herbaceous
  &          180.0,      &       ! cultivated & managed areas
  &          200.0,      &       ! mosaic crop / tree / natural vegetation
  &          150.0,      &       ! mosaic crop / shrub / grass
  &          150.0,      &       ! bare areas                       
  &          150.0,      &       ! water
  &          150.0,      &       ! snow & ice 
  &          150.0,      &       ! artificial surface   
  &          0.         /)       ! undefined

!---------------------------------------------------------------------------------------------- 
!----------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------- 
!----------------------------------------------------------------------------------------------


REAL (KIND=wp) :: z0c_experimental(nclass_glc2000)   = (/ &       !< lookup table landuse class to roughness length [m]
  &          1.00,       &       ! evergreen broadleaf forest   
  &          1.00,       &       ! deciduous broadleaf closed forest
  &          0.15,       &       ! deciduous broadleaf open   forest
  &          1.00,       &       ! evergreen needleleaf forest   
  &          1.00,       &       ! deciduous needleleaf forest
  &          1.00,       &       ! mixed leaf trees            
  &          1.00,       &       ! fresh water flooded trees
  &          1.00,       &       ! saline water flooded trees
  &          0.20,       &       ! mosaic tree / natural vegetation
  &          0.05,       &       ! burnt tree cover
  &          0.20,       &       ! evergreen shrubs closed-open
  &          0.15,       &       ! decidous shrubs closed-open
  &          0.03,       &       ! herbaceous vegetation closed-open
  &          0.05,       &       ! sparse herbaceous or grass 
  &          0.05,       &       ! flooded shrubs or herbaceous
  &          0.07,       &       ! cultivated & managed areas
  &          0.25,       &       ! mosaic crop / tree / natural vegetation
  &          0.07,       &       ! mosaic crop / shrub / grass
  &          0.05,       &       ! bare areas                       
  &          0.0002,     &       ! water
  &          0.01,       &       ! snow & ice 
  &          1.00,       &       ! artificial surface  
  &          0.         /)        ! undefined

REAL (KIND=wp) :: zplcmnc_experimental(nclass_glc2000)  = (/ &      !< lookup table landuse class to minimal plant cover
  &         0.80,  &       ! evergreen broadleaf forest   
  &         0.75,  &       ! deciduous broadleaf closed forest
  &         0.70,  &       ! deciduous broadleaf open   forest
  &         0.80,  &       ! evergreen needleleaf forest   
  &         0.75,  &       ! deciduous needleleaf forest
  &         0.75,  &       ! mixed leaf trees            
  &         0.80,  &       ! fresh water flooded trees
  &         0.80,  &       ! saline water flooded trees
  &         0.70,  &       ! mosaic tree / natural vegetation
  &         0.50,  &       ! burnt tree cover
  &         0.70,  &       ! evergreen shrubs closed-open
  &         0.70,  &       ! decidous shrubs closed-open
  &         0.75,  &       ! herbaceous vegetation closed-open
  &         0.50,  &       ! sparse herbaceous or grass 
  &         0.70,  &       ! flooded shrubs or herbaceous
  &         0.50,  &       ! cultivated & managed areas
  &         0.65,  &       ! mosaic crop / tree / natural vegetation
  &         0.50,  &       ! mosaic crop / shrub / grass
  &         0.02,  &       ! bare areas                         
  &         0.00,  &       ! water
  &         0.00,  &       ! snow & ice 
  &         0.10,  &       ! artificial surface  
  &         0.       /)    ! undefined

REAL (KIND=wp) :: zplcmxc_experimental(nclass_glc2000)  = (/ &     !< lookup table landuse class to maximal plant cover
 &            0.80, &       ! evergreen broadleaf forest   
 &            0.90, &       ! deciduous broadleaf closed forest
 &            0.80, &       ! deciduous broadleaf open   forest
 &            0.80, &       ! evergreen needleleaf forest   
 &            0.90, &       ! deciduous needleleaf forest
 &            0.90, &       ! mixed leaf trees            
 &            0.80, &       ! fresh water flooded trees
 &            0.80, &       ! saline water flooded trees
 &            0.80, &       ! mosaic tree / natural vegetation
 &            0.50, &       ! burnt tree cover
 &            0.80, &       ! evergreen shrubs closed-open
 &            0.80, &       ! decidous shrubs closed-open
 &            0.90, &       ! herbaceous vegetation closed-open
 &            0.50, &       ! sparse herbaceous or grass 
 &            0.80, &       ! flooded shrubs or herbaceous
 &            0.90, &       ! cultivated & managed areas
 &            0.80, &       ! mosaic crop / tree / natural vegetation
 &            0.90, &       ! mosaic crop / shrub / grass
 &            0.05, &       ! bare areas                       
 &            0.00, &       ! water
 &            0.00, &       ! snow & ice 
 &            0.20, &       ! artificial surface   
 &            0.       /)   ! undefined

 !< lookup table landuse class to maximal plant cover
REAL (KIND=wp) :: zlaimnc_experimental(nclass_glc2000)   = (/ &      !< lookup table landuse class to minimal leaf area index
 &          1.40, &        ! evergreen broadleaf forest   
 &          1.00, &        ! deciduous broadleaf closed forest
 &          1.00, &        ! deciduous broadleaf open   forest
 &          1.30, &        ! evergreen needleleaf forest   
 &          1.00, &        ! deciduous needleleaf forest
 &          1.10, &        ! mixed leaf trees            
 &          1.40, &        ! fresh water flooded trees
 &          1.40, &        ! saline water flooded trees
 &          0.60, &        ! mosaic tree / natural vegetation
 &          0.40, &        ! burnt tree cover
 &          0.60, &        ! evergreen shrubs closed-open
 &          1.00, &        ! decidous shrubs closed-open
 &          1.00, &        ! herbaceous vegetation closed-open
 &          0.40, &        ! sparse herbaceous or grass 
 &          1.00, &        ! flooded shrubs or herbaceous
 &          0.70, &        ! cultivated & managed areas
 &          1.00, &        ! mosaic crop / tree / natural vegetation
 &          0.70, &        ! mosaic crop / shrub / grass
 &          0.40, &        ! bare areas                       
 &          0.00, &        ! water
 &          0.00, &        ! snow & ice 
 &          0.10, &        ! artificial surface   
 &          0.       /)    ! undefined

REAL (KIND=wp) :: zlaimxc_experimental(nclass_glc2000)  = (/ &      !< lookup table landuse class to maximal leaf area index
 &            5.00, &       ! evergreen broadleaf forest   
 &            6.00, &       ! deciduous broadleaf closed forest
 &            4.00, &       ! deciduous broadleaf open   forest
 &            5.00, &       ! evergreen needleleaf forest   
 &            5.00, &       ! deciduous needleleaf forest
 &            5.00, &       ! mixed leaf trees            
 &            5.00, &       ! fresh water flooded trees
 &            5.00, &       ! saline water flooded trees
 &            2.50, &       ! mosaic tree / natural vegetation
 &            0.60, &       ! burnt tree cover
 &            3.00, &       ! evergreen shrubs closed-open
 &            1.50, &       ! decidous shrubs closed-open
 &            3.10, &       ! herbaceous vegetation closed-open
 &            0.60, &       ! sparse herbaceous or grass 
 &            2.00, &       ! flooded shrubs or herbaceous
 &            3.30, &       ! cultivated & managed areas
 &            3.00, &       ! mosaic crop / tree / natural vegetation
 &            3.50, &       ! mosaic crop / shrub / grass
 &            0.60, &       ! bare areas                       
 &            0.00, &       ! water
 &            0.00, &       ! snow & ice 
 &            1.00, &       ! artificial surface   
 &            0.      /)    ! undefined

REAL (KIND=wp) :: zrd_experimental(nclass_glc2000)   = (/ &         !< lookup table landuse class to root depth [m]
  &          1.00,       &       ! evergreen broadleaf forest   
  &          1.00,       &       ! deciduous broadleaf closed forest
  &          2.00,       &       ! deciduous broadleaf open   forest
  &          0.60,       &       ! evergreen needleleaf forest   
  &          0.60,       &       ! deciduous needleleaf forest
  &          0.80,       &       ! mixed leaf trees            
  &          1.00,       &       ! fresh water flooded trees
  &          1.00,       &       ! saline water flooded trees
  &          1.00,       &       ! mosaic tree / natural vegetation
  &          0.30,       &       ! burnt tree cover
  &          1.00,       &       ! evergreen shrubs closed-open
  &          2.00,       &       ! decidous shrubs closed-open
  &          0.60,       &       ! herbaceous vegetation closed-open
  &          0.30,       &       ! sparse herbaceous or grass 
  &          0.40,       &       ! flooded shrubs or herbaceous
  &          1.00,       &       ! cultivated & managed areas
  &          1.00,       &       ! mosaic crop / tree / natural vegetation
  &          1.00,       &       ! mosaic crop / shrub / grass
  &          0.30,       &       ! bare areas                       
  &          0.00,       &       ! water
  &          0.00,       &       ! snow & ice 
  &          0.60,       &       ! artificial surface   
  &          0.         /)       ! undefined

REAL (KIND=wp) :: zemiss_experimental(nclass_glc2000)  = (/ &       !< lookup table landuse class to surface thermal emissivity
  &          0.996,      &       ! evergreen broadleaf forest   
  &          0.990,      &       ! deciduous broadleaf closed forest
  &          0.993,      &       ! deciduous broadleaf open   forest
  &          0.996,      &       ! evergreen needleleaf forest   
  &          0.990,      &       ! deciduous needleleaf forest
  &          0.993,      &       ! mixed leaf trees            
  &          0.996,      &       ! fresh water flooded trees
  &          0.996,      &       ! saline water flooded trees
  &          0.985,      &       ! mosaic tree / natural vegetation
  &          0.950,      &       ! burnt tree cover
  &          0.985,      &       ! evergreen shrubs closed-open
  &          0.993,      &       ! decidous shrubs closed-open
  &          0.993,      &       ! herbaceous vegetation closed-open
  &          0.950,      &       ! sparse herbaceous or grass 
  &          0.992,      &       ! flooded shrubs or herbaceous
  &          0.990,      &       ! cultivated & managed areas
  &          0.990,      &       ! mosaic crop / tree / natural vegetation
  &          0.990,      &       ! mosaic crop / shrub / grass
  &          0.950,      &       ! bare areas                       
  &          0.991,      &       ! water
  &          0.9999,     &       ! snow & ice 
  &          0.960 ,     &       ! artificial surface   
  &          0.         /)       ! undefined

REAL (KIND=wp) :: zrs_min_experimental(nclass_glc2000) =(/ &
  &          250.0,      &       ! evergreen broadleaf forest   
  &          150.0,      &       ! deciduous broadleaf closed forest
  &          150.0,      &       ! deciduous broadleaf open   forest
  &          150.0,      &       ! evergreen needleleaf forest   
  &          150.0,      &       ! deciduous needleleaf forest
  &          150.0,      &       ! mixed leaf trees            
  &          150.0,      &       ! fresh water flooded trees
  &          150.0,      &       ! saline water flooded trees
  &          150.0,      &       ! mosaic tree / natural vegetation
  &          150.0,      &       ! burnt tree cover
  &          120.0,      &       ! evergreen shrubs closed-open
  &          120.0,      &       ! decidous shrubs closed-open
  &          40.0,      &       ! herbaceous vegetation closed-open
  &          40.0,      &       ! sparse herbaceous or grass 
  &          40.0,      &       ! flooded shrubs or herbaceous
  &          120.0,      &       ! cultivated & managed areas
  &          120.0,      &       ! mosaic crop / tree / natural vegetation
  &          100.0,      &       ! mosaic crop / shrub / grass
  &          120.0,      &       ! bare areas                       
  &          120.0,      &       ! water
  &          120.0,      &       ! snow & ice 
  &          120.0,     &       ! artificial surface   
  &          0.         /)       ! undefined



!> legend of the GLC2000 vegetation classes
CHARACTER(len=29) :: glc2000_legend(nclass_glc2000) = (/&    ! No.        IGBP correspondence
 &  'evergreen broadleaf tree     ' , &    ! 1.         evergreen broadleaf forest
 &  'deciduous broadleaf tree clos' , &    ! 2.         deciduous broadleaf forest
 &  'deciduous broadleaf tree open' , &    ! 3.         woody savannas & savanna 
 &  'evergreen needleleaf tree    ' , &    ! 4.         evergreen needleleaf forest
 &  'deciduous needleleaf tree    ' , &    ! 5.         deciduous needleleaf forest
 &  'mixed leaf tree              ' , &    ! 6.         mixed forest
 &  'fresh water flooded tree     ' , &    ! 7.         evergreen broadleaf
 &  'saline water flooded tree    ' , &    ! 8.         evergreen broadleaf
 &  'mosaic tree / other nat.veg. ' , &    ! 9.         shrubland ?
 &  'burnt tree cover             ' , &    ! 10.
 &  'evergreen shrubs closed-open ' , &    ! 11.        shrubland open-closed
 &  'deciduous shrubs closed-open ' , &    ! 12.        savannas
 &  'herbaceous cover closed-open ' , &    ! 13.        grassland
 &  'sparse herbaceous or grass   ' , &    ! 14.
 &  'flooded shrub or herbaceous  ' , &    ! 15.        wetlands
 &  'cultivated & managed areas   ' , &    ! 16.        croplands
 &  'mosaic crop/tree/natural veg.' , &    ! 17.        cropland/other vegetation mosaic
 &  'mosaic crop/shrub or grass   ' , &    ! 18.        cropland/other vegetation mosaic
 &  'bare areas                   ' , &    ! 19.        barren or sparsely vegetated
 &  'water bodies                 ' , &    ! 20.        water
 &  'snow & ice                   ' , &    ! 21.        snow & ice
 &  'artificial surfaces          ' , &    ! 22.        urban and built-up areas
 &  'undefined                    ' /)     ! 23.




CONTAINS

  !> define lookup table for GLC2000 landuse classes
  SUBROUTINE init_glc2000_lookup_tables(nclass_glc2000, &
    &      ilookup_table_glc2000, &
    &      z0_lt_glc2000,           &
    &      lnz0_lt_glc2000,       &
    &      plc_mn_lt_glc2000,        &
    &      plc_mx_lt_glc2000,        &
    &      lai_mn_lt_glc2000,        &
    &      lai_mx_lt_glc2000,        &
    &      rd_lt_glc2000,          &
    &      emiss_lt_glc2000,       &
    &      rs_min_lt_glc2000)
    INTEGER, INTENT(IN) :: nclass_glc2000 !< GLC2000 has 23 classes for the land use description
    INTEGER, INTENT(IN) :: ilookup_table_glc2000  !< integer switch to choose a lookup table
    REAL (KIND=wp), INTENT(OUT) :: z0_lt_glc2000(nclass_glc2000)      !< lookup table landuse class to roughness length [m]
    REAL (KIND=wp), INTENT(OUT) :: lnz0_lt_glc2000(nclass_glc2000)    !< corresponding natural logarithm of z0c_gme_o
    REAL (KIND=wp), INTENT(OUT) :: plc_mn_lt_glc2000(nclass_glc2000)  !< lookup table LU class to minimal plant cover
    REAL (KIND=wp), INTENT(OUT) :: plc_mx_lt_glc2000(nclass_glc2000)  !< lookup table LU class to maximal plant cover
    REAL (KIND=wp), INTENT(OUT) :: lai_mn_lt_glc2000(nclass_glc2000)  !< lookup table LU class to minimal leaf area index
    REAL (KIND=wp), INTENT(OUT) :: lai_mx_lt_glc2000(nclass_glc2000)  !< lookup table LU class to maximal leaf area index
    REAL (KIND=wp), INTENT(OUT) :: rd_lt_glc2000(nclass_glc2000)      !< lookup table LU class to root depth [m]
    REAL (KIND=wp), INTENT(OUT) :: emiss_lt_glc2000(nclass_glc2000)   !< lookup table LU class to surface thermal emissivity
    REAL (KIND=wp), INTENT(OUT) :: rs_min_lt_glc2000(nclass_glc2000)  !< lookup table LU class to minimal stomata resistance

    ! local variable
    INTEGER :: i !< counter
    REAL :: arg

      SELECT CASE (ilookup_table_glc2000)
        CASE(i_gme_lookup_table)
           z0_lt_glc2000 = z0c_gme_o
           plc_mn_lt_glc2000 = zplcmnc_gme_o
           plc_mx_lt_glc2000 = zplcmxc_gme_o
           lai_mn_lt_glc2000 = zlaimnc_gme_o
           lai_mx_lt_glc2000 = zlaimxc_gme_o
           rd_lt_glc2000 = zrd_gme_o
           emiss_lt_glc2000 = zemiss_gme_o
           rs_min_lt_glc2000 = zrs_min_gme_o
        CASE(i_cosmo_lookup_table)
           z0_lt_glc2000 = z0c_cosmo_o
           plc_mn_lt_glc2000 = zplcmnc_cosmo_o
           plc_mx_lt_glc2000 = zplcmxc_cosmo_o
           lai_mn_lt_glc2000 = zlaimnc_cosmo_o
           lai_mx_lt_glc2000 = zlaimxc_cosmo_o
           rd_lt_glc2000 = zrd_cosmo_o
           emiss_lt_glc2000 = zemiss_cosmo_o
           rs_min_lt_glc2000 = zrs_min_cosmo_o
        CASE(i_experimental_lookup_table)
           z0_lt_glc2000 = z0c_experimental
           plc_mn_lt_glc2000 = zplcmnc_experimental
           plc_mx_lt_glc2000 = zplcmxc_experimental
           lai_mn_lt_glc2000 = zlaimnc_experimental
           lai_mx_lt_glc2000 = zlaimxc_experimental
           rd_lt_glc2000 = zrd_experimental
           emiss_lt_glc2000 = zemiss_experimental
           rs_min_lt_glc2000 = zrs_min_experimental
        CASE DEFAULT
           z0_lt_glc2000 = z0c_gme_o
           plc_mn_lt_glc2000 = zplcmnc_gme_o
           plc_mx_lt_glc2000 = zplcmxc_gme_o
           lai_mn_lt_glc2000 = zlaimnc_gme_o
           lai_mx_lt_glc2000 = zlaimxc_gme_o
           rd_lt_glc2000 = zrd_gme_o
           emiss_lt_glc2000 = zemiss_gme_o
           rs_min_lt_glc2000 = zrs_min_gme_o
      END SELECT

      DO i=1,nclass_glc2000
        IF (z0_lt_glc2000(i) > 0.) THEN
          arg = z0_lt_glc2000(i)
          lnz0_lt_glc2000(i) = ALOG(arg)
        ENDIF
      ENDDO

  END  SUBROUTINE init_glc2000_lookup_tables

  
  !> define  name of lookup table for GLC2000 
  SUBROUTINE get_name_glc2000_lookup_tables(ilookup_table_glc2000, name_lookup_table_glc2000)
    INTEGER, INTENT(IN) :: ilookup_table_glc2000  !< integer switch to choose a lookup table
    CHARACTER (LEN=filename_max), INTENT(OUT) :: name_lookup_table_glc2000 !< name of lookup table
    ! local variable
      SELECT CASE (ilookup_table_glc2000)
        CASE(i_gme_lookup_table)
           name_lookup_table_glc2000='Ritter_2007'
        CASE(i_cosmo_lookup_table)
           name_lookup_table_glc2000='Heise_2005'
        CASE(i_experimental_lookup_table)
           name_lookup_table_glc2000='Asensio_2010'
        CASE DEFAULT
           name_lookup_table_glc2000='Ritter_2007'
      END SELECT

  END  SUBROUTINE get_name_glc2000_lookup_tables



   !> assign the GLC2000 land use classes to some characteristic (more or less) physical parameters
  SUBROUTINE glc2000_look_up(lu, &
    &      nclass_glc2000, &
    &      lnz0_lt_glc2000,          &
    &      plc_mn_lt_glc2000,        &
    &      plc_mx_lt_glc2000,        &
    &      lai_mn_lt_glc2000,        &
    &      lai_mx_lt_glc2000,        &
    &      rd_lt_glc2000,            &
    &      emiss_lt_glc2000,         &
    &      rs_min_lt_glc2000,        &
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
  INTEGER, INTENT(IN) :: nclass_glc2000 !< GLC2000 has 23 classes for the land use description
  REAL (KIND=wp), INTENT(IN) :: lnz0_lt_glc2000(nclass_glc2000)    !< corresponding natural logarithm of z0c_gme_o
  REAL (KIND=wp), INTENT(IN) :: plc_mn_lt_glc2000(nclass_glc2000)  !< lookup table landuse class to minimal plant cover
  REAL (KIND=wp), INTENT(IN) :: plc_mx_lt_glc2000(nclass_glc2000)  !< lookup table landuse class to maximal plant cover
  REAL (KIND=wp), INTENT(IN) :: lai_mn_lt_glc2000(nclass_glc2000)  !< lookup table landuse class to minimal leaf area index
  REAL (KIND=wp), INTENT(IN) :: lai_mx_lt_glc2000(nclass_glc2000)  !< lookup table landuse class to maximal leaf area index
  REAL (KIND=wp), INTENT(IN) :: rd_lt_glc2000(nclass_glc2000)      !< lookup table landuse class to root depth [m]
  REAL (KIND=wp), INTENT(IN) :: emiss_lt_glc2000(nclass_glc2000)   !< lookup table landuse class to surface thermal emissivity
  REAL (KIND=wp), INTENT(IN) :: rs_min_lt_glc2000(nclass_glc2000)  !< lookup table landuse class to minimal stomata resistance

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
          IF (lu>=1 .AND. lu<=22 .AND.lu/=20) THEN
            k_error     = 0
            pland       = 1.0
            plnz0       = lnz0_lt_glc2000(lu)
            pmn         = plc_mn_lt_glc2000(lu)
            pmx         = plc_mx_lt_glc2000(lu)
            plaimn      = lai_mn_lt_glc2000(lu)
            plaimx      = lai_mx_lt_glc2000(lu)
            proot       = rd_lt_glc2000(lu)
            prs_min     = rs_min_lt_glc2000(lu)     
            pemissivity = emiss_lt_glc2000(lu)
            purb    = 0.0
            pfor_d  = 0.0
            pfor_e  = 0.0
            pice    = 0.0

            IF (lu==22            ) purb   = 1.0  ! urban pixel
            IF (lu== 2 .OR. lu== 5) pfor_d = 1.0  ! dec.forest 
            IF (lu== 1 .OR. lu== 4) pfor_e = 1.0  ! eve.forest 
            IF (lu== 6) THEN                      ! mix.forest
              pfor_d = 0.5  
              pfor_e = 0.5  
            END IF
            IF (lu==21            ) pice   = 1.0  ! ice or snow pixel
          ELSE IF (lu==20) THEN ! water
            k_error     = 0
            pland       = 0.0
            pemissivity = emiss_lt_glc2000(lu)             ! emissivity is required everywhere
          ELSE
            k_error     = 1  ! not a valid land use class
            pland       = 0.0 
          END IF

  END  SUBROUTINE glc2000_look_up






END MODULE mo_glc2000_lookup_tables
