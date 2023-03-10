!==========================================================
!+ Fortran Module with lookup-tables for the ecosg data
!!
!! Description:
!! The ECOCLIMAP-SG dataset contains the following land use
!! classification scheme
! class=value        description
!-----------------------------------------------------
! 01 'sea and oceans                                      '
! 02 'lakes                                               '
! 03 'rivers                                              '
!-----------------------------------------------------
! 04 'bare land                                           '
! 05 'bare rock                                           '
! 06 'permanent snow                                      '
! 07 'boreal broadleaf deciduous                          '
! 08 'temperate broadleaf deciduous                       '
! 09 'tropical broadleaf deciduous                        '
! 10 'temperate broadleaf evergreen                       '
! 11 'tropical broadleaf evergreen                        '
! 12 'boreal needleleaf evergreen                         '
! 13 'temperate needleleaf evergreen                      '
! 14 'boreal needleleaf deciduous                         '
! 15 'shrubs                                              '
! 16 'boreal grassland                                    '
! 17 'temperate grassland                                 '
! 18 'tropical grassland                                  '
! 19 'winter C3 crops (lower temp. & greater water avail.)'
! 20 'summer C3 crops                                     '
! 21 'C4 crops (warmer environments)                      '
! 22 'flooded trees                                       '
! 23 'flooded grassland                                   '
!-----------------------------------------------------
! 24 'LCZ1: compact high-rise                             '
! 25 'LCZ2: compact midrise                               '
! 26 'LCZ3: compact low-rise                              '
! 27 'LCZ4: open high-rise                                '
! 28 'LCZ5: open midrise                                  '
! 29 'LCZ6: open low-rise                                 '
! 30 'LCZ7: lightweight low-rise                          '
! 31 'LCZ8: large low-rise                                '
! 32 'LCZ9: sparsely built                                '
! 33 'LCZ10: heavy industry                               '
!-----------------------------------------------------
!! lookup-table compiled by Andrzej Wyszogrodzki, Adam Jaczewski
!-----------------------------------------------------
MODULE mo_ecosg_lookup_tables

  USE mo_kind,                  ONLY: wp, i4
  USE mo_io_units,              ONLY: filename_max

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_ecosg_lookup_tables, &
       &    get_name_ecosg_lookup_tables, &
       &    ecosg_look_up, &
       &    get_ecosg_idx, &
       &    ecosg_legend, &
       &    ecosg_value, &
       &    nclass_ecosg, &
       &    i_extpar_lookup_table, &
       &    name_lookup_table_ecosg, &
       &    z0_lt_ecosg, &
       &    lnz0_lt_ecosg, &
       &    plc_mn_lt_ecosg, &
       &    plc_mx_lt_ecosg, &
       &    lai_mn_lt_ecosg, &
       &    lai_mx_lt_ecosg, &
       &    rd_lt_ecosg, &
       &    skinc_lt_ecosg, &
       &    emiss_lt_ecosg, &
       &    rs_min_lt_ecosg

  INTEGER (KIND=i4), PARAMETER :: nclass_ecosg = 33, &  !< ecosg has 23 classes for the land use description
       &                          i_extpar_lookup_table = 1 !< lookup_table for ecosg land use classes

  CHARACTER (LEN=filename_max) :: name_lookup_table_ecosg !< name of lookup table

  REAL (KIND=wp)               :: z0_lt_ecosg(nclass_ecosg), &       !< lookup table landuse class to roughness length [m]
       &                          lnz0_lt_ecosg(nclass_ecosg), &     !< corresponding natural logarithm of z0c_extpar_o
       &                          plc_mn_lt_ecosg(nclass_ecosg), &   !< lookup table landuse class to minimal plant cover
       &                          plc_mx_lt_ecosg(nclass_ecosg), &   !< lookup table landuse class to maximal plant cover
       &                          lai_mn_lt_ecosg(nclass_ecosg), &   !< lookup table landuse class to minimal leaf area index
       &                          lai_mx_lt_ecosg(nclass_ecosg), &   !< lookup table landuse class to maximal leaf area index
       &                          rd_lt_ecosg(nclass_ecosg), &       !< lookup table landuse class to root depth [m]
       &                          skinc_lt_ecosg(nclass_ecosg), &    !< lookup table landuse class to skin conductivity [W m-2 K-1]
       &                          emiss_lt_ecosg(nclass_ecosg), &    !< lookup table landuse class to surface thermal emissivity
       &                          rs_min_lt_ecosg(nclass_ecosg)  !< lookup table landuse class to minimal stomata resistance

  !----------------------------------------------------------------------------------------------
  ! Tables based on:
  ! Globcover2009
  ! Oke et al. 2017 / Demuzere et al. 2019
  ! WRF landuse tavls: surface roughbess = 0.8
  !----------------------------------------------------------------------------------------------

  REAL (KIND=wp) :: z0c_extpar_o(nclass_ecosg)  = (/ & !< lookup table landuse class to roughness length [m]
   &  0.0002, & ! 01 'sea and oceans                   '
   &  0.0002, & ! 02 'lakes                            '
   &  0.0002, & ! 03 'rivers --------------------------'
   &  0.05  , & ! 04 'bare land                        '
   &  0.05  , & ! 05 'bare rock                        '
   &  0.01  , & ! 06 'permanent snow                   '
   &  1.0   , & ! 07 'boreal broadleaf deciduous       '
   &  0.575 , & ! 08 'temperate broadleaf deciduous    '
   &  1.0   , & ! 09 'tropical broadleaf deciduous     '
   &  0.15  , & ! 10 'temperate broadleaf evergreen    '
   &  1.0   , & ! 11 'tropical broadleaf evergreen     '
   &  1.0   , & ! 12 'boreal needleleaf evergreen      '
   &  1.0   , & ! 13 'temperate needleleaf evergreen   '
   &  1.0   , & ! 14 'boreal needleleaf deciduous      '
   &  0.15  , & ! 15 'shrubs                           '
   &  0.03  , & ! 16 'boreal grassland                 '
   &  0.03  , & ! 17 'temperate grassland              '
   &  0.07  , & ! 18 'tropical grassland               '
   &  0.07  , & ! 19 'winter C3 crops                  '
   &  0.07  , & ! 20 'summer C3 crops                  '
   &  1.0   , & ! 21 'C4 crops (warmer environments)   '
   &  0.03  , & ! 22 'flooded trees                    '
   &  0.05  , & ! 23 'flooded grassland ---------------'
   &  2.5   , & ! 24 'LCZ1: compact high-rise          ' ! WRF Landuse TBL - > Oke modified
   &  1.5   , & ! 25 'LCZ2: compact midrise            ' ! WRF Landuse TBL - > Oke modified
   &  0.8   , & ! 26 'LCZ3: compact low-rise           ' ! WRF Landuse TBL - > Oke modified
   &  2.5   , & ! 27 'LCZ4: open high-rise             ' ! WRF Landuse TBL - > Oke modified
   &  1.5   , & ! 28 'LCZ5: open midrise               ' ! WRF Landuse TBL - > Oke modified
   &  0.8   , & ! 29 'LCZ6: open low-rise              ' ! WRF Landuse TBL - > Oke modified
   &  0.3   , & ! 30 'LCZ7: lightweight low-rise       ' ! WRF Landuse TBL - > Oke modified
   &  0.8   , & ! 31 'LCZ8: large low-rise             ' ! WRF Landuse TBL - > Oke modified
   &  0.8   , & ! 32 'LCZ9: sparsely built             ' ! WRF Landuse TBL - > Oke modified
   &  1.0 /)    ! 33 'LCZ10: heavy industry            ' ! WRF Landuse TBL - > Oke modified
   !------------------------------------------------------
   ! WRF Landuse TBL = 0.8 unchanged when Oke in range 3-10
   !------------------------------------------------------
   ! modification by Oke 2017: roughness = surface elements x 0.1
   ! > 25
   ! 25-10
   !  3-10
   ! > 25
   ! 10-25
   ! 3-10
   ! 2-4
   ! 3-10
   ! 3-10
   ! 5-15
   !------------------------------------------------------

  REAL (KIND=wp) :: zplcmnc_extpar_o(nclass_ecosg) = (/ &      !< lookup table landuse class to minimal plant cover
   &  0.0   , & ! 01 'sea and oceans                   '
   &  0.0   , & ! 02 'lakes                            '
   &  0.0   , & ! 03 'rivers --------------------------'
   &  0.0   , & ! 04 'bare land                        '
   &  0.0   , & ! 05 'bare rock                        '
   &  0.75  , & ! 06 'permanent snow                   '
   &  0.725 , & ! 07 'boreal broadleaf deciduous       '
   &  0.7   , & ! 08 'temperate broadleaf deciduous    '
   &  0.8   , & ! 09 'tropical broadleaf deciduous     '
   &  0.8   , & ! 10 'temperate broadleaf evergreen    '
   &  0.8   , & ! 11 'tropical broadleaf evergreen     '
   &  0.8   , & ! 12 'boreal needleleaf evergreen      '
   &  0.0   , & ! 13 'temperate needleleaf evergreen   '
   &  0.75  , & ! 14 'boreal needleleaf deciduous      '
   &  0.70  , & ! 15 'shrubs                           '
   &  0.75  , & ! 16 'boreal grassland                 '
   &  0.75  , & ! 17 'temperate grassland              '
   &  0.75  , & ! 18 'tropical grassland               '
   &  0.5   , & ! 19 'winter C3 crops                  '
   &  0.5   , & ! 20 'summer C3 crops                  '
   &  0.5   , & ! 21 'C4 crops (warmer environments)   '
   &  0.8   , & ! 22 'flooded trees                    '
   &  0.5   , & ! 23 'flooded grassland ---------------'
   &  0.40  , & ! 24 'LCZ1: compact high-rise          ' ! Oke 2017 fraction minimal [%]
   &  0.30  , & ! 25 'LCZ2: compact midrise            ' ! Oke 2017 fraction minimal [%]
   &  0.20  , & ! 26 'LCZ3: compact low-rise           ' ! Oke 2017 fraction minimal [%]
   &  0.30  , & ! 27 'LCZ4: open high-rise             ' ! Oke 2017 fraction minimal [%]
   &  0.30  , & ! 28 'LCZ5: open midrise               ' ! Oke 2017 fraction minimal [%]
   &  0.20  , & ! 29 'LCZ6: open low-rise              ' ! Oke 2017 fraction minimal [%]
   &  0.10  , & ! 30 'LCZ7: lightweight low-rise       ' ! Oke 2017 fraction minimal [%]
   &  0.40  , & ! 31 'LCZ8: large low-rise             ' ! Oke 2017 fraction minimal [%]
   &  0.10  , & ! 32 'LCZ9: sparsely built             ' ! Oke 2017 fraction minimal [%]
   &  0.20 /)   ! 33 'LCZ10: heavy industry            ' ! Oke 2017 fraction minimal [%]
   !----------------------------------
   ! Oke 2017 Imprevious Plant fraction [%]
   !----------------------------------
   ! 40-60
   ! 30-50
   ! 20-50
   ! 30-40
   ! 30-50
   ! 20-50
   ! <20
   ! 40-50
   ! <20
   ! 20-40
   !----------------------------------

  REAL (KIND=wp) :: zplcmxc_extpar_o(nclass_ecosg) = (/ &     !< lookup table landuse class to maximal plant cover
   &  0.0  , & ! 01 'sea and oceans                   '
   &  0.0  , & ! 02 'lakes                            '
   &  0.0  , & ! 03 'rivers --------------------------'
   &  0.05 , & ! 04 'bare land                        '
   &  0.05 , & ! 05 'bare rock                        '
   &  0.0  , & ! 06 'permanent snow                   '
   &  0.9  , & ! 07 'boreal broadleaf deciduous       '
   &  0.85 , & ! 08 'temperate broadleaf deciduous    '
   &  0.8  , & ! 09 'tropical broadleaf deciduous     '
   &  0.8  , & ! 10 'temperate broadleaf evergreen    '
   &  0.8  , & ! 11 'tropical broadleaf evergreen     '
   &  0.8  , & ! 12 'boreal needleleaf evergreen      '
   &  0.8  , & ! 13 'temperate needleleaf evergreen   '
   &  0.9  , & ! 14 'boreal needleleaf deciduous      '
   &  0.8  , & ! 15 'shrubs                           '
   &  0.9  , & ! 16 'boreal grassland                 '
   &  0.9  , & ! 17 'temperate grassland              '
   &  0.9  , & ! 18 'tropical grassland               '
   &  0.9  , & ! 19 'winter C3 crops                  '
   &  0.9  , & ! 20 'summer C3 crops                  '
   &  0.9  , & ! 21 'C4 crops (warmer environments)   '
   &  0.8  , & ! 22 'flooded trees                    '
   &  0.8  , & ! 23 'flooded grassland ---------------'
   &  0.6  , & ! 24 'LCZ1: compact high-rise          ' ! Oke 2017 fraction maxmal [%]
   &  0.5  , & ! 25 'LCZ2: compact midrise            ' ! Oke 2017 fraction maxmal [%]
   &  0.5  , & ! 26 'LCZ3: compact low-rise           ' ! Oke 2017 fraction maxmal [%]
   &  0.4  , & ! 27 'LCZ4: open high-rise             ' ! Oke 2017 fraction maxmal [%]
   &  0.5  , & ! 28 'LCZ5: open midrise               ' ! Oke 2017 fraction maxmal [%]
   &  0.5  , & ! 29 'LCZ6: open low-rise              ' ! Oke 2017 fraction maxmal [%]
   &  0.15 , & ! 30 'LCZ7: lightweight low-rise       ' ! Oke 2017 fraction maxmal [%]
   &  0.5  , & ! 31 'LCZ8: large low-rise             ' ! Oke 2017 fraction maxmal [%]
   &  0.15 , & ! 32 'LCZ9: sparsely built             ' ! Oke 2017 fraction maxmal [%]
   &  0.4 /)   ! 33 'LCZ10: heavy industry            ' ! Oke 2017 fraction maxmal [%]
   !----------------------------------
   ! Oke 2017 Imprevious Plant fraction [%]
   !----------------------------------
   ! 40-60
   ! 30-50
   ! 20-50
   ! 30-40
   ! 30-50
   ! 20-50
   ! <20
   ! 40-50
   ! <20
   ! 20-40
   !----------------------------------

  REAL (KIND=wp) :: zlaimnc_extpar_o(nclass_ecosg) = (/ &      !< lookup table landuse class to minimal leaf area index
   &  0.0  , & ! 01 'sea and oceans                   '
   &  0.0  , & ! 02 'lakes                            '
   &  0.0  , & ! 03 'rivers --------------------------'
   &  0.4  , & ! 04 'bare land                        '
   &  0.4  , & ! 05 'bare rock                        '
   &  0.0  , & ! 06 'permanent snow                   '
   &  1.0  , & ! 07 'boreal broadleaf deciduous       '
   &  1.0  , & ! 08 'temperate broadleaf deciduous    '
   &  1.0  , & ! 09 'tropical broadleaf deciduous     '
   &  1.4  , & ! 10 'temperate broadleaf evergreen    '
   &  1.4  , & ! 11 'tropical broadleaf evergreen     '
   &  1.3  , & ! 12 'boreal needleleaf evergreen      '
   &  1.3  , & ! 13 'temperate needleleaf evergreen   '
   &  1.0  , & ! 14 'boreal needleleaf deciduous      '
   &  0.6  , & ! 15 'shrubs                           '
   &  1.0  , & ! 16 'boreal grassland                 '
   &  1.0  , & ! 17 'temperate grassland              '
   &  1.0  , & ! 18 'tropical grassland               '
   &  0.7  , & ! 19 'winter C3 crops                  '
   &  0.7  , & ! 20 'summer C3 crops                  '
   &  0.7  , & ! 21 'C4 crops (warmer environments)   '
   &  1.4  , & ! 22 'flooded trees                    '
   &  1.0  , & ! 23 'flooded grassland ---------------'
   &  0.1  , & ! 24 'LCZ1: compact high-rise          ' ! dense mix of buildings, few or no trees
   &  0.2  , & ! 25 'LCZ2: compact midrise            ' ! dense mix of buildings, few or no trees
   &  0.4  , & ! 26 'LCZ3: compact low-rise           ' ! dense mix of buildings, few or no trees
   &  0.8  , & ! 27 'LCZ4: open high-rise             ' ! opne/tall buildings, low plants, scattered trees
   &  1.0  , & ! 28 'LCZ5: open midrise               ' ! open/mid buildings, low plants, scattered trees
   &  1.2  , & ! 29 'LCZ6: open low-rise              ' ! open/low buildings, low plants, scattered trees
   &  0.8  , & ! 30 'LCZ7: lightweight low-rise       ' ! dense/single story, few or no trees
   &  0.3  , & ! 31 'LCZ8: large low-rise             ' ! open/lowrise, few or no trees
   &  1.3  , & ! 32 'LCZ9: sparsely built             ' ! sparse/mid buildings, low plants, scattered trees
   &  0.1  /)  ! 33 'LCZ10: heavy industry            ' ! low/mid industrial, few or no trees


  REAL (KIND=wp) :: zlaimxc_extpar_o(nclass_ecosg) = (/ &      !< lookup table landuse class to maximal leaf area index
   &  0.0  , & ! 01 'sea and oceans                   '
   &  0.0  , & ! 02 'lakes                            '
   &  0.0  , & ! 03 'rivers --------------------------'
   &  0.6  , & ! 04 'bare land                        '
   &  0.6  , & ! 05 'bare rock                        '
   &  0.0  , & ! 06 'permanent snow                   '
   &  3.4  , & ! 07 'boreal broadleaf deciduous       '
   &  2.7  , & ! 08 'temperate broadleaf deciduous    '
   &  2.0  , & ! 09 'tropical broadleaf deciduous     '
   &  2.4  , & ! 10 'temperate broadleaf evergreen    '
   &  2.4  , & ! 11 'tropical broadleaf evergreen     '
   &  3.8  , & ! 12 'boreal needleleaf evergreen      '
   &  3.8  , & ! 13 'temperate needleleaf evergreen   '
   &  3.8  , & ! 14 'boreal needleleaf deciduous      '
   &  1.5  , & ! 15 'shrubs                           '
   &  3.1  , & ! 16 'boreal grassland                 '
   &  3.1  , & ! 17 'temperate grassland              '
   &  3.1  , & ! 18 'tropical grassland               '
   &  3.3  , & ! 19 'winter C3 crops                  '
   &  3.3  , & ! 20 'summer C3 crops                  '
   &  3.3  , & ! 21 'C4 crops (warmer environments)   '
   &  2.4  , & ! 22 'flooded trees                    '
   &  2.0  , & ! 23 'flooded grassland ---------------'
   &  0.2  , & ! 24 'LCZ1: compact high-rise          ' ! dense mix of buildings, few or no trees
   &  0.35 , & ! 25 'LCZ2: compact midrise            ' ! dense mix of buildings, few or no trees
   &  0.5  , & ! 26 'LCZ3: compact low-rise           ' ! dense mix of buildings, few or no trees
   &  1.2  , & ! 27 'LCZ4: open high-rise             ' ! opne/tall buildings, low plants, scattered trees
   &  1.4  , & ! 28 'LCZ5: open midrise               ' ! open/mid buildings, low plants, scattered trees
   &  1.5  , & ! 29 'LCZ6: open low-rise              ' ! open/low buildings, low plants, scattered trees
   &  1.2  , & ! 30 'LCZ7: lightweight low-rise       ' ! dense/single story, few or no trees
   &  0.6  , & ! 31 'LCZ8: large low-rise             ' ! open/lowrise, few or no trees
   &  1.8  , & ! 32 'LCZ9: sparsely built             ' ! sparse/mid buildings, low plants, scattered trees
   &  0.5  /)  ! 33 'LCZ10: heavy industry            ' ! low/mid industrial, few or no trees


  REAL (KIND=wp) :: zrd_extpar_o(nclass_ecosg)  = (/ &         !< lookup table landuse class to root depth [m]
   &  0.0  , & ! 01 'sea and oceans                   '
   &  0.0  , & ! 02 'lakes                            '
   &  0.0  , & ! 03 'rivers --------------------------'
   &  0.3  , & ! 04 'bare land                        '
   &  0.3  , & ! 05 'bare rock                        '
   &  0.0  , & ! 06 'permanent snow                   '
   &  1.0  , & ! 07 'boreal broadleaf deciduous       '
   &  1.5  , & ! 08 'temperate broadleaf deciduous    '
   &  2.0  , & ! 09 'tropical broadleaf deciduous     '
   &  1.0  , & ! 10 'temperate broadleaf evergreen    '
   &  1.0  , & ! 11 'tropical broadleaf evergreen     '
   &  0.6  , & ! 12 'boreal needleleaf evergreen      '
   &  0.6  , & ! 13 'temperate needleleaf evergreen   '
   &  0.6  , & ! 14 'boreal needleleaf deciduous      '
   &  1.5  , & ! 15 'shrubs                           '
   &  0.6  , & ! 16 'boreal grassland                 '
   &  0.6  , & ! 17 'temperate grassland              '
   &  0.6  , & ! 18 'tropical grassland               '
   &  1.0  , & ! 19 'winter C3 crops                  '
   &  1.0  , & ! 20 'summer C3 crops                  '
   &  1.0  , & ! 21 'C4 crops (warmer environments)   '
   &  1.0  , & ! 22 'flooded trees                    '
   &  1.0  , & ! 23 'flooded grassland ---------------'
   &  0.6  , & ! 24 'LCZ1: compact high-rise          ' ! dense mix of buildings, few or no trees
   &  0.7  , & ! 25 'LCZ2: compact midrise            ' ! dense mix of buildings, few or no trees
   &  0.8  , & ! 26 'LCZ3: compact low-rise           ' ! dense mix of buildings, few or no trees
   &  1.6  , & ! 27 'LCZ4: open high-rise             ' ! opne/tall buildings, low plants, scattered trees
   &  1.4  , & ! 28 'LCZ5: open midrise               ' ! open/mid buildings, low plants, scattered trees
   &  1.2  , & ! 29 'LCZ6: open low-rise              ' ! open/low buildings, low plants, scattered trees
   &  1.0  , & ! 30 'LCZ7: lightweight low-rise       ' ! dense/single story, few or no trees
   &  0.8  , & ! 31 'LCZ8: large low-rise             ' ! open/lowrise, few or no trees
   &  2.0  , & ! 32 'LCZ9: sparsely built             ' ! sparse/mid buildings, low plants, scattered trees
   &  0.6  /)  ! 33 'LCZ10: heavy industry            ' ! low/mid industrial, few or no trees
   !-------------------------------------------
   ! Oke 2017 = 0.6 - 2.0
   !-------------------------------------------

  REAL (KIND=wp) :: zskinc_extpar_o(nclass_ecosg) = (/ &       !< lookup table landuse class to skin conductivity
   & 200.0  , & ! 01 'sea and oceans                   '
   & 200.0  , & ! 02 'lakes                            '
   & 200.0  , & ! 03 'rivers --------------------------'
   & 200.0  , & ! 04 'bare land                        '
   & 200.0  , & ! 05 'bare rock                        '
   & 200.0  , & ! 06 'permanent snow                   '
   &  50.0  , & ! 07 'boreal broadleaf deciduous       '
   &  40.0  , & ! 08 'temperate broadleaf deciduous    '
   &  30.0  , & ! 09 'tropical broadleaf deciduous     '
   &  50.0  , & ! 10 'temperate broadleaf evergreen    '
   &  50.0  , & ! 11 'tropical broadleaf evergreen     '
   &  50.0  , & ! 12 'boreal needleleaf evergreen      '
   &  50.0  , & ! 13 'temperate needleleaf evergreen   '
   &  50.0  , & ! 14 'boreal needleleaf deciduous      '
   &  50.0  , & ! 15 'shrubs                           '
   &  30.0  , & ! 16 'boreal grassland                 '
   &  30.0  , & ! 17 'temperate grassland              '
   &  30.0  , & ! 18 'tropical grassland               '
   &  30.0  , & ! 19 'winter C3 crops                  '
   &  30.0  , & ! 20 'summer C3 crops                  '
   &  30.0  , & ! 21 'C4 crops (warmer environments)   '
   &  50.0  , & ! 22 'flooded trees                    '
   &  30.0  , & ! 23 'flooded grassland ---------------'
   & 200.0  , & ! 24 'LCZ1: compact high-rise          ' ! GlobCover 200 W/m2 -> antropogenic heat flux
   &  75.0  , & ! 25 'LCZ2: compact midrise            ' ! GlobCover 200 W/m2 -> antropogenic heat flux
   &  75.0  , & ! 26 'LCZ3: compact low-rise           ' ! GlobCover 200 W/m2 -> antropogenic heat flux
   &  50.0  , & ! 27 'LCZ4: open high-rise             ' ! GlobCover 200 W/m2 -> antropogenic heat flux
   &  25.0  , & ! 28 'LCZ5: open midrise               ' ! GlobCover 200 W/m2 -> antropogenic heat flux
   &  25.0  , & ! 29 'LCZ6: open low-rise              ' ! GlobCover 200 W/m2 -> antropogenic heat flux
   &  35.0  , & ! 30 'LCZ7: lightweight low-rise       ' ! GlobCover 200 W/m2 -> antropogenic heat flux
   &  50.0  , & ! 31 'LCZ8: large low-rise             ' ! GlobCover 200 W/m2 -> antropogenic heat flux
   &  10.0  , & ! 32 'LCZ9: sparsely built             ' ! GlobCover 200 W/m2 -> antropogenic heat flux
   & 300.0 /)   ! 33 'LCZ10: heavy industry            ' ! GlobCover 200 W/m2 -> antropogenic heat flux
   !-------------------------------------------
   ! Oke 2017 antropogenic heat flux [W/m2]
   !-------------------------------------------
   ! 50-300
   ! < 75
   ! <75
   ! <50
   ! <25
   ! <25
   ! <35
   ! <50
   ! <10
   ! >300
   !-------------------------------------------

  REAL (KIND=wp) :: zemiss_extpar_o(nclass_ecosg) = (/ &       !< lookup table landuse class to surface thermal emissivity
   & 0.991  , & ! 01 'sea and oceans                   '
   & 0.991  , & ! 02 'lakes                            '
   & 0.991  , & ! 03 'rivers --------------------------'
   & 0.950  , & ! 04 'bare land                        '
   & 0.950  , & ! 05 'bare rock                        '
   & 0.990  , & ! 06 'permanent snow                   '
   & 0.9915 , & ! 07 'boreal broadleaf deciduous       '
   & 0.9999 , & ! 08 'temperate broadleaf deciduous    '
   & 0.993  , & ! 09 'tropical broadleaf deciduous     '
   & 0.996  , & ! 10 'temperate broadleaf evergreen    '
   & 0.996  , & ! 11 'tropical broadleaf evergreen     '
   & 0.996  , & ! 12 'boreal needleleaf evergreen      '
   & 0.996  , & ! 13 'temperate needleleaf evergreen   '
   & 0.990  , & ! 14 'boreal needleleaf deciduous      '
   & 0.990  , & ! 15 'shrubs                           '
   & 0.993  , & ! 16 'boreal grassland                 '
   & 0.993  , & ! 17 'temperate grassland              '
   & 0.993  , & ! 18 'tropical grassland               '
   & 0.990  , & ! 19 'winter C3 crops                  '
   & 0.990  , & ! 20 'summer C3 crops                  '
   & 0.990  , & ! 21 'C4 crops (warmer environments)   '
   & 0.996  , & ! 22 'flooded trees                    '
   & 0.992  , & ! 23 'flooded grassland ---------------'
   & 0.860  , & ! 24 'LCZ1: compact high-rise          ' ! dense mix of buildings, few or no trees
   & 0.880  , & ! 25 'LCZ2: compact midrise            ' ! dense mix of buildings, few or no trees
   & 0.900  , & ! 26 'LCZ3: compact low-rise           ' ! dense mix of buildings, few or no trees
   & 0.920  , & ! 27 'LCZ4: open high-rise             ' ! opne/tall buildings, low plants, scattered trees
   & 0.940  , & ! 28 'LCZ5: open midrise               ' ! open/mid buildings, low plants, scattered trees
   & 0.960  , & ! 29 'LCZ6: open low-rise              ' ! open/low buildings, low plants, scattered trees
   & 0.900  , & ! 30 'LCZ7: lightweight low-rise       ' ! dense/single story, few or no trees
   & 0.890  , & ! 31 'LCZ8: large low-rise             ' ! open/lowrise, few or no trees
   & 0.980  , & ! 32 'LCZ9: sparsely built             ' ! sparse/mid buildings, low plants, scattered trees
   & 0.850 /)   ! 33 'LCZ10: heavy industry            ' ! low/mid industrial, few or no trees
   !---------------------------------
   ! Oke 2017:
   !---------------------------------
   ! low vege, grass   = 0.9-0.98
   ! roads asphalt     = 0.89-0.96
   ! concrete          = 0.85-0.97
   ! rofs/shingles     = 0.90-0.92
   !------------------------------------

  REAL (KIND=wp) :: zrs_min_extpar_o(nclass_ecosg) = (/ &      !< lookup table landuse class to minimal stomata resistance
   &  120.  , & ! 01 'sea and oceans                   '
   &  120.  , & ! 02 'lakes                            '
   &  120.  , & ! 03 'rivers --------------------------'
   &  120.  , & ! 04 'bare land                        '
   &  120.  , & ! 05 'bare rock                        '
   &  120.  , & ! 06 'permanent snow                   '
   &  150.  , & ! 07 'boreal broadleaf deciduous       '
   &  150.  , & ! 08 'temperate broadleaf deciduous    '
   &  250.  , & ! 09 'tropical broadleaf deciduous     '
   &  250.  , & ! 10 'temperate broadleaf evergreen    '
   &  150.  , & ! 11 'tropical broadleaf evergreen     '
   &  150.  , & ! 12 'boreal needleleaf evergreen      '
   &  150.  , & ! 13 'temperate needleleaf evergreen   '
   &  150.  , & ! 14 'boreal needleleaf deciduous      '
   &  120.  , & ! 15 'shrubs                           '
   &   40.  , & ! 16 'boreal grassland                 '
   &   40.  , & ! 17 'temperate grassland              '
   &   40.  , & ! 18 'tropical grassland               '
   &  120.  , & ! 19 'winter C3 crops                  '
   &  120.  , & ! 20 'summer C3 crops                  '
   &  120.  , & ! 21 'C4 crops (warmer environments)   '
   &   40.  , & ! 22 'flooded trees                    '
   &  150.  , & ! 23 'flooded grassland ---------------'
   & 1200.  , & ! 24 'LCZ1: compact high-rise          ' ! dense mix of buildings, few or no trees
   &  800.  , & ! 25 'LCZ2: compact midrise            ' ! dense mix of buildings, few or no trees
   &  400.  , & ! 26 'LCZ3: compact low-rise           ' ! dense mix of buildings, few or no trees
   &  250.  , & ! 27 'LCZ4: open high-rise             ' ! opne/tall buildings, low plants, scattered trees
   &  100.  , & ! 28 'LCZ5: open midrise               ' ! open/mid buildings, low plants, scattered trees
   &   50.  , & ! 29 'LCZ6: open low-rise              ' ! open/low buildings, low plants, scattered trees
   &  300.  , & ! 30 'LCZ7: lightweight low-rise       ' ! dense/single story, few or no trees
   &  600.  , & ! 31 'LCZ8: large low-rise             ' ! open/lowrise, few or no trees
   &   30.  , & ! 32 'LCZ9: sparsely built             ' ! sparse/mid buildings, low plants, scattered trees
   & 1200. /)   ! 33 'LCZ10: heavy industry            ' ! low/mid industrial, few or no trees
   !---------------------------------
   ! Oke 2017:
   !---------------------------------
   ! open warer        = 0-5
   ! grasslan          = 30
   ! plants, lil water = 50-100
   ! open-lowrise, wet = 250 (watered gardens)
   ! open-lowrise, dry = 1200
   !------------------------------------


  !> legend of the ecosg vegetation classes
  CHARACTER(len=33) :: ecosg_legend(nclass_ecosg) = (/&    ! No.
   &  'sea and oceans                   ' , &    ! 1.
   &  'lakes                            ' , &    ! 2.
   &  'rivers --------------------------' , &    ! 3.
   &  'bare land                        ' , &    ! 4.
   &  'bare rock                        ' , &    ! 5.
   &  'permanent snow                   ' , &    ! 6.
   &  'boreal broadleaf deciduous       ' , &    ! 7.
   &  'temperate broadleaf deciduous    ' , &    ! 8.
   &  'tropical broadleaf deciduous     ' , &    ! 9.
   &  'temperate broadleaf evergreen    ' , &    ! 10.
   &  'tropical broadleaf evergreen     ' , &    ! 11.
   &  'boreal needleleaf evergreen      ' , &    ! 12.
   &  'temperate needleleaf evergreen   ' , &    ! 13.
   &  'boreal needleleaf deciduous      ' , &    ! 14.
   &  'shrubs                           ' , &    ! 15.
   &  'boreal grassland                 ' , &    ! 16.
   &  'temperate grassland              ' , &    ! 17.
   &  'tropical grassland               ' , &    ! 18.
   &  'winter C3 crops                  ' , &    ! 19.
   &  'summer C3 crops                  ' , &    ! 20.
   &  'C4 crops (warmer environments)   ' , &    ! 21.
   &  'flooded trees                    ' , &    ! 22.
   &  'flooded grassland ---------------' , &    ! 23.
   &  'LCZ1: compact high-rise          ' , &    ! 24.
   &  'LCZ2: compact midrise            ' , &    ! 25.
   &  'LCZ3: compact low-rise           ' , &    ! 26.
   &  'LCZ4: open high-rise             ' , &    ! 27.
   &  'LCZ5: open midrise               ' , &    ! 28.
   &  'LCZ6: open low-rise              ' , &    ! 29.
   &  'LCZ7: lightweight low-rise       ' , &    ! 30.
   &  'LCZ8: large low-rise             ' , &    ! 31.
   &  'LCZ9: sparsely built             ' , &    ! 32.
   &  'LCZ10: heavy industry            ' /)     ! 33.

  !> values of the ecosg landuse classes
  INTEGER :: ecosg_value(nclass_ecosg) =          (/&    ! No.
   &   10  , & ! 01 'sea and oceans                   '
   &   20  , & ! 02 'lakes                            '
   &   30  , & ! 03 'rivers --------------------------'
   &   40  , & ! 04 'bare land                        '
   &   50  , & ! 05 'bare rock                        '
   &   60  , & ! 06 'permanent snow                   '
   &   70  , & ! 07 'boreal broadleaf deciduous       '
   &   80  , & ! 08 'temperate broadleaf deciduous    '
   &   90  , & ! 09 'tropical broadleaf deciduous     '
   &  100  , & ! 10 'temperate broadleaf evergreen    '
   &  110  , & ! 11 'tropical broadleaf evergreen     '
   &  120  , & ! 12 'boreal needleleaf evergreen      '
   &  130  , & ! 13 'temperate needleleaf evergreen   '
   &  140  , & ! 14 'boreal needleleaf deciduous      '
   &  150  , & ! 15 'shrubs                           '
   &  160  , & ! 16 'boreal grassland                 '
   &  170  , & ! 17 'temperate grassland              '
   &  180  , & ! 18 'tropical grassland               '
   &  190  , & ! 19 'winter C3 crops                  '
   &  200  , & ! 20 'summer C3 crops                  '
   &  210  , & ! 21 'C4 crops (warmer environments)   '
   &  220  , & ! 22 'flooded trees                    '
   &  230  , & ! 23 'flooded grassland ---------------'
   &  240  , & ! 24 'LCZ1: compact high-rise          '
   &  250  , & ! 25 'LCZ2: compact midrise            '
   &  260  , & ! 26 'LCZ3: compact low-rise           '
   &  270  , & ! 27 'LCZ4: open high-rise             '
   &  280  , & ! 28 'LCZ5: open midrise               '
   &  290  , & ! 29 'LCZ6: open low-rise              '
   &  300  , & ! 30 'LCZ7: lightweight low-rise       '
   &  310  , & ! 31 'LCZ8: large low-rise             '
   &  320  , & ! 32 'LCZ9: sparsely built             '
   &  330 /)   ! 33 'LCZ10: heavy industry            '

  CONTAINS

  !> define lookup table for ecosg landuse classes
  SUBROUTINE init_ecosg_lookup_tables(nclass_ecosg,      &
       &                            ilookup_table_ecosg, &
       &                            z0_lt_ecosg,         &
       &                            lnz0_lt_ecosg,       &
       &                            plc_mn_lt_ecosg,     &
       &                            plc_mx_lt_ecosg,     &
       &                            lai_mn_lt_ecosg,     &
       &                            lai_mx_lt_ecosg,     &
       &                            rd_lt_ecosg,         &
       &                            skinc_lt_ecosg,      &
       &                            emiss_lt_ecosg,      &
       &                            rs_min_lt_ecosg)

    INTEGER(KIND=i4), INTENT(IN) :: nclass_ecosg, &  !< ecosg has 23 classes for the land use description
         &                          ilookup_table_ecosg  !< integer switch to choose a lookup table

    REAL (KIND=wp), INTENT(OUT)  :: z0_lt_ecosg(nclass_ecosg), &       !< lookup table LU class to roughness length [m]
         &                          lnz0_lt_ecosg(nclass_ecosg), &     !< corresponding natural logarithm of z0c_extpar_o
         &                          plc_mn_lt_ecosg(nclass_ecosg), &   !< lookup table landuse class to minimal plant cover
         &                          plc_mx_lt_ecosg(nclass_ecosg), &   !< lookup table landuse class to maximal plant cover
         &                          lai_mn_lt_ecosg(nclass_ecosg), &   !< lookup table landuse class to minimal LAI
         &                          lai_mx_lt_ecosg(nclass_ecosg), &   !< lookup table landuse class to maximal LAI
         &                          rd_lt_ecosg(nclass_ecosg), &       !< lookup table LU class to root depth [m]
         &                          skinc_lt_ecosg(nclass_ecosg), &    !< lookup table landuse class to skin conductivity
         &                          emiss_lt_ecosg(nclass_ecosg), &    !< lookup table LU class to surface thermal emiss.
         &                          rs_min_lt_ecosg(nclass_ecosg)      !< lookup table LU class to minimal stomata resis.
    ! local variables
    INTEGER(KIND=i4)              :: i !< counter
    REAL(KIND=wp)                 :: arg

    SELECT CASE (ilookup_table_ecosg)
      CASE(i_extpar_lookup_table)
             z0_lt_ecosg =     z0c_extpar_o
         plc_mn_lt_ecosg = zplcmnc_extpar_o
         plc_mx_lt_ecosg = zplcmxc_extpar_o
         lai_mn_lt_ecosg = zlaimnc_extpar_o
         lai_mx_lt_ecosg = zlaimxc_extpar_o
             rd_lt_ecosg =     zrd_extpar_o
          skinc_lt_ecosg =  zskinc_extpar_o
          emiss_lt_ecosg =  zemiss_extpar_o
         rs_min_lt_ecosg = zrs_min_extpar_o
      CASE DEFAULT
             z0_lt_ecosg =     z0c_extpar_o
         plc_mn_lt_ecosg = zplcmnc_extpar_o
         plc_mx_lt_ecosg = zplcmxc_extpar_o
         lai_mn_lt_ecosg = zlaimnc_extpar_o
         lai_mx_lt_ecosg = zlaimxc_extpar_o
             rd_lt_ecosg =     zrd_extpar_o
          skinc_lt_ecosg =  zskinc_extpar_o
          emiss_lt_ecosg =  zemiss_extpar_o
         rs_min_lt_ecosg = zrs_min_extpar_o
    END SELECT

    lnz0_lt_ecosg = 0.
    DO i=1,nclass_ecosg
      IF (z0_lt_ecosg(i) > 0.) THEN
        arg = z0_lt_ecosg(i)
        lnz0_lt_ecosg(i) = LOG(arg)
      ENDIF
    ENDDO

  END  SUBROUTINE init_ecosg_lookup_tables

  !> define  name of lookup table for ecosg
  SUBROUTINE get_name_ecosg_lookup_tables(ilookup_table_ecosg, name_lookup_table_ecosg)

    INTEGER(KIND=i4), INTENT(IN)              :: ilookup_table_ecosg  !< integer switch to choose a lookup table

    CHARACTER (LEN=filename_max), INTENT(OUT) :: name_lookup_table_ecosg !< name of lookup table

      SELECT CASE (ilookup_table_ecosg)
        CASE(i_extpar_lookup_table)
         name_lookup_table_ecosg='Globcover_SG Modified'
        CASE DEFAULT
         name_lookup_table_ecosg='Globcover_SG Modified'
      END SELECT

  END  SUBROUTINE get_name_ecosg_lookup_tables

  !> assign the ecosg land use classes to some characteristic (more or less) physical parameters
  SUBROUTINE ecosg_look_up(lu, &
       &                   nclass_ecosg, &
       &                   lnz0_lt_ecosg,          &
       &                   plc_mn_lt_ecosg,        &
       &                   plc_mx_lt_ecosg,        &
       &                   lai_mn_lt_ecosg,        &
       &                   lai_mx_lt_ecosg,        &
       &                   rd_lt_ecosg,            &
       &                   skinc_lt_ecosg,         &
       &                   emiss_lt_ecosg,         &
       &                   rs_min_lt_ecosg,        &
       &                   pland,          &
       &                   pice,           &
       &                   plnz0,          &
       &                   proot,          &
       &                   pmn,            &
       &                   pmx,            &
       &                   plaimn,         &
       &                   plaimx,         &
       &                   purb,           &
       &                   pfor_d,         &
       &                   pfor_e,         &
       &                   pskinc,         &
       &                   pemissivity,    &
       &                   prs_min,        &
       &                   k_error)

    INTEGER(KIND=i4), INTENT(IN)  :: lu, &              !< land use class
         &                           nclass_ecosg !< ecosg has 23 classes for the land use description + 10 urban classes

    REAL (KIND=wp), INTENT(IN)    :: lnz0_lt_ecosg(nclass_ecosg), &     !< corresponding natural logarithm of z0c_extpar_o
         &                           plc_mn_lt_ecosg(nclass_ecosg), &   !< lookup table landuse class to minimal plant cover
         &                           plc_mx_lt_ecosg(nclass_ecosg), &   !< lookup table landuse class to maximal plant cover
         &                           lai_mn_lt_ecosg(nclass_ecosg), &   !< lookup table landuse class to minimal LAI
         &                           lai_mx_lt_ecosg(nclass_ecosg), &   !< lookup table landuse class to maximal LAI
         &                           rd_lt_ecosg(nclass_ecosg), &       !< lookup table landuse class to root depth [m]
         &                           skinc_lt_ecosg(nclass_ecosg), &    !< lookup table landuse class to skin conductivity [W m-2 K-1]
         &                           emiss_lt_ecosg(nclass_ecosg), &    !< lookup table landuse class to surface thermal emiss.
         &                           rs_min_lt_ecosg(nclass_ecosg)  !< lookup table landuse class to minimal stomata resis.

    REAL (KIND=wp), INTENT(OUT)   :: pland, &           !< land cover                      (-)
         &                           pice, &            !< ice fraction                    (-)
         &                           plnz0, &           !< logarithm of roughness length   (m)
         &                           proot, &           !< root depth                      (m)
         &                           pmn, &             !< minimal plant cover             (-)
         &                           pmx, &             !< maximum plant cover             (-)
         &                           plaimn, &          !< minimal leaf area index         (m**2/m**2)
         &                           plaimx, &          !< maximum leaf area index         (m**2/m**2)
         &                           purb, &            !< urbanisation                    (-)
         &                           pfor_d, &          !< deciduous forest                (-)
         &                           pfor_e, &          !< evergreen forest                (-)
         &                           pskinc, &          !< skin conductivity               (W m-2 K-1)
         &                           pemissivity, &     !< surface thermal emissivity      (-)
         &                           prs_min        !< minimum stomata resistance      (s/m)

    INTEGER(KIND=i4), INTENT(OUT) :: k_error     !< error return code

    ! local variables
    INTEGER(KIND=i4)              :: nclass !< position of landuse class in arrays

!   CALL get_ecosg_idx(lu,nclass)
    nclass=lu

    ! Test for true land points
!   IF (nclass>=1 .AND. nclass<=22 .AND.nclass/=21) THEN
    IF (nclass>=1 .AND. nclass<=33) THEN
      k_error     = 0
      plnz0       = lnz0_lt_ecosg(nclass)
      pmn         = plc_mn_lt_ecosg(nclass)
      pmx         = plc_mx_lt_ecosg(nclass)
      plaimn      = lai_mn_lt_ecosg(nclass)
      plaimx      = lai_mx_lt_ecosg(nclass)
      proot       = rd_lt_ecosg(nclass)
      pskinc      = skinc_lt_ecosg(nclass)
      prs_min     = rs_min_lt_ecosg(nclass)
      pemissivity = emiss_lt_ecosg(nclass)

      pland   = 1.0
      purb    = 0.0
      pfor_d  = 0.0
      pfor_e  = 0.0
      pice    = 0.0

      IF (lu>=24) purb  = 1.0  ! urban surfaces
      IF (lu<=3) pland  = 0.0  ! water body (sea, lake, river)
      IF (lu==6)  pice  = 1.0  ! permanent snow

!----------------------------------------------------------
! 07 'boreal broadleaf deciduous                          '
! 08 'temperate broadleaf deciduous                       '
! 09 'tropical broadleaf deciduous                        '
! 14 'boreal needleleaf deciduous                         '
!----------------------------------------------------------
      IF (lu== 7  .OR. lu== 8  .OR. lu== 9  .OR. lu== 14) pfor_d = 1.0  ! deciduous forest

!----------------------------------------------------------
! 10 'temperate broadleaf evergreen                       '
! 11 'tropical broadleaf evergreen                        '
! 12 'boreal needleleaf evergreen                         '
! 13 'temperate needleleaf evergreen                      '
!----------------------------------------------------------
      IF (lu== 10 .OR. lu== 11 .OR. lu== 12 .OR. lu== 13) pfor_e = 1.0  ! evergreen forest

!     IF (lu== 10) THEN                      ! mixed forest
!       pfor_d = 0.5
!       pfor_e = 0.5
!     END IF

    ELSE
      k_error     = 1  ! not a valid land use class
      pland       = 0.0
      pskinc      = 200.0
    END IF

  END  SUBROUTINE ecosg_look_up


  SUBROUTINE get_ecosg_idx(lu,nclass)

    INTEGER(KIND=i4), INTENT(IN)  :: lu     !< land use class
    INTEGER(KIND=i4), INTENT(OUT) :: nclass !< position of landuse class in arrays

    SELECT CASE(lu)
      CASE (1) ! ecosg_value(1)
        nclass = 1
      CASE (2) !ecosg_value(2)
        nclass = 2
      CASE (3) ! ecosg_value(3)
        nclass = 3
      CASE (4) ! ecosg_value(4)
        nclass = 4
      CASE (5) ! ecosg_value(5)
        nclass = 5
      CASE (6) ! ecosg_value(6)
        nclass = 6
      CASE (7) ! ecosg_value(7)
        nclass = 7
      CASE (8) ! ecosg_value(8)
        nclass = 8
      CASE (9) ! ecosg_value(9)
        nclass = 9
      CASE (10) ! ecosg_value(10)
        nclass = 10
      CASE (11) ! ecosg_value(11)
        nclass = 11
      CASE (12) ! ecosg_value(12)
        nclass = 12
      CASE (13) ! ecosg_value(13)
        nclass = 13
      CASE (14) ! ecosg_value(14)
        nclass = 14
      CASE (15) ! ecosg_value(15)
        nclass = 15
      CASE (16) ! ecosg_value(16)
        nclass = 16
      CASE (17) ! ecosg_value(17)
        nclass = 17
      CASE (18) ! ecosg_value(18)
        nclass = 18
      CASE (19) ! ecosg_value(19)
        nclass = 19
      CASE (20) ! ecosg_value(20)
        nclass = 20
      CASE (21) ! ecosg_value(21)
        nclass = 21
      CASE (22) ! ecosg_value(22)
        nclass = 22
      CASE (23) ! ecosg_value(23)
        nclass = 23

      CASE (24) ! ecosg_value(24) ! LCZ1
        nclass = 24
      CASE (25) ! ecosg_value(25) ! LCZ2
        nclass = 25
      CASE (26) ! ecosg_value(26) ! LCZ3
        nclass = 26
      CASE (27) ! ecosg_value(27) ! LCZ4
        nclass = 27
      CASE (28) ! ecosg_value(28) ! LCZ5
        nclass = 28
      CASE (29) ! ecosg_value(29) ! LCZ6
        nclass = 29
      CASE (30) ! ecosg_value(30) ! LCZ7
        nclass = 30
      CASE (31) ! ecosg_value(31) ! LCZ8
        nclass = 31
      CASE (32) ! ecosg_value(32) ! LCZ9
        nclass = 32
      CASE (33) ! ecosg_value(33) ! LCZ10
        nclass = 33

      CASE DEFAULT
        nclass = 33 !
    END SELECT

  END SUBROUTINE get_ecosg_idx

END MODULE mo_ecosg_lookup_tables
