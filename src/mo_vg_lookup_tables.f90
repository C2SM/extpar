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

  USE mo_io_units,              ONLY: filename_max

  IMPLICIT NONE


  PUBLIC :: init_vg_lookup_tables, & 
       &    get_name_vg_lookup_tables, & 
       &    vg_look_up, & 
       &    vg_legend, & 
       &    ntype_fao, & 
       &    name_lookup_table_vg, & 
       &    z0_lt_glcc, lnz0_lt_glcc, plc_mn_lt_glcc, plc_mx_lt_glcc , & 
       &    lai_mn_lt_glcc, lai_mx_lt_glcc, rd_lt_glcc, emiss_lt_glcc, rs_min_lt_glcc         



  INTEGER (KIND=i4), PARAMETER :: ntype_fao = 10  !< FAO has 10 soil types

  REAL (KIND=wp)               :: cporv_vg(ntype_fao),       &     !< pore volume (fraction of volume)
       &                          cfcap_vg(ntype_fao),       &     !< field capacity (fraction of volume)
       &                          cpwp_vg(ntype_fao),        &     !< plant wilting point (fraction of volume)
       &                          cadp_vg(ntype_fao),        &     !< air dryness point (fraction of volume)
       &                          ckw0_vg(ntype_fao),        &     !< parameter for hydrological conductivity
       &                          n_vg(ntype_fao),           &     !< n parameter in van Genuchten equation
       &                          alpha_vg(ntype_fao)              !< alpha parameter in van Genuchten equation

  CHARACTER (LEN=filename_max) :: name_lookup_table_vg !< name of lookup table


  !---------------------------------------------------------------------------------------------- 
  !----------------------------------------------------------------------------------------------
  REAL (KIND=wp) :: cporv_vg(ntype_fao)  = (/ &       !< lookup table pore volume (fraction of volume)
    &   1.E-10,          &
    &   1.E-10,          &
    &   0.43,            &
    &   0.41,            &
    &   0.43,            &
    &   0.41,            &
    &   0.507,           &
    &   0.73,            &
    &   1.E-10,          &
    &   1.E-10_wp /)

  REAL (KIND=wp) :: cfcap_vg(ntype_fao)  = (/ &       !< lookup table field capacita (fraction of volume)
    &   1.E-10,          &
    &   1.E-10,          &
    &   0.196,           &
    &   0.260,           &
    &   0.340,           &
    &   0.370,           &
    &   0.463,           &
    &   0.70,            &
    &   1.E-10,          &
    &   1.E-10           /)

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
    &    1.E-10,         &
    &    1.E-10,         &
    &    0.145,          &
    &    0.075,          &
    &    0.036,          &
    &    0.019,          &
    &    0.008,          &
    &    0.01,           &
    &    1.E-10,         &
    &    1.E-10          /)

  CONTAINS

  !> define lookup table for GLCC landuse classes
  SUBROUTINE init_glcc_lookup_tables(nclass_glcc, &
       &                             ilookup_table_glcc, &
       &                             z0_lt_glcc,           &
       &                             lnz0_lt_glcc,       &
       &                             plc_mn_lt_glcc,        &
       &                             plc_mx_lt_glcc,        &
       &                             lai_mn_lt_glcc,        &
       &                             lai_mx_lt_glcc,        &
       &                             rd_lt_glcc,          &
       &                             emiss_lt_glcc,       &
       &                             rs_min_lt_glcc)

    INTEGER(KIND=i4), INTENT(IN) :: nclass_glcc, &  !< GLCC has 24 classes for the land use description
         &                          ilookup_table_glcc  !< integer switch to choose a lookup table

    REAL (KIND=wp), INTENT(OUT)  :: z0_lt_glcc(nclass_glcc), &       !< lookup table landuse class to roughness length [m]
         &                          lnz0_lt_glcc(nclass_glcc), &     !< corresponding natural logarithm of z0c_gme_o
         &                          plc_mn_lt_glcc(nclass_glcc), &   !< lookup table landuse class to minimal plant cover
         &                          plc_mx_lt_glcc(nclass_glcc), &   !< lookup table landuse class to maximal plant cover
         &                          lai_mn_lt_glcc(nclass_glcc), &   !< lookup table landuse class to minimal leaf area index
         &                          lai_mx_lt_glcc(nclass_glcc), &   !< lookup table landuse class to maximal leaf area index
         &                          rd_lt_glcc(nclass_glcc), &       !< lookup table landuse class to root depth [m]
         &                          emiss_lt_glcc(nclass_glcc), &    !< lookup table landuse class to surface thermal emissivity
         &                          rs_min_lt_glcc(nclass_glcc)  !< lookup table landuse class to minimal stomata resistance

    ! local variables
    INTEGER(KIND=i4)             :: i !< counter
    REAL(KIND=wp)                :: arg

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
        lnz0_lt_glcc(i) = LOG(arg)
      ENDIF
    ENDDO

  END  SUBROUTINE init_glcc_lookup_tables

    !> define  name of lookup table for GLCc 
  SUBROUTINE get_name_glcc_lookup_tables(ilookup_table_glcc, name_lookup_table_glcc)

    INTEGER(KIND=i4), INTENT(IN)              :: ilookup_table_glcc  !< integer switch to choose a lookup table
    
    CHARACTER (LEN=filename_max), INTENT(OUT) :: name_lookup_table_glcc !< name of lookup table

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
                   
    INTEGER, INTENT(IN)         :: lu, &              !< land use class
         &                         nclass_glcc !< GLCC has 24 classes for the land use description
                                
    REAL (KIND=wp), INTENT(IN)  :: lnz0_lt_glcc(nclass_glcc), &     !< corresponding natural logarithm of z0c_gme_o
         &                         plc_mn_lt_glcc(nclass_glcc), &   !< lookup table landuse class to minimal plant cover
         &                         plc_mx_lt_glcc(nclass_glcc), &   !< lookup table landuse class to maximal plant cover
         &                         lai_mn_lt_glcc(nclass_glcc), &   !< lookup table landuse class to minimal leaf area index
         &                         lai_mx_lt_glcc(nclass_glcc), &   !< lookup table landuse class to maximal leaf area index
         &                         rd_lt_glcc(nclass_glcc), &       !< lookup table landuse class to root depth [m]
         &                         emiss_lt_glcc(nclass_glcc), &    !< lookup table landuse class to surface thermal emissivity
         &                         rs_min_lt_glcc(nclass_glcc)  !< lookup table landuse class to minimal stomata resistance

    REAL (KIND=wp), INTENT(OUT) :: pland, &           !< land cover                      (-)
         &                         pice, &            !< ice fraction                    (-)
         &                         plnz0, &           !< logarithm of roughness length   (m)
         &                         proot, &           !< root depth                      (m)
         &                         pmn, &             !< minimal plant cover             (-)
         &                         pmx, &             !< maximum plant cover             (-)
         &                         plaimn, &          !< minimal leaf area index         (m**2/m**2)
         &                         plaimx, &          !< maximum leaf area index         (m**2/m**2)
         &                         purb, &            !< urbanisation                    (-)
         &                         pfor_d, &          !< deciduous forest                (-)
         &                         pfor_e, &          !< evergreen forest                (-)
         &                         pemissivity, &     !< surface thermal emissivity      (-)
         &                         prs_min        !< minimum stomata resistance      (s/m)

    INTEGER(KIND=i4), INTENT(OUT):: k_error     !< error return code
  
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

END MODULE mo_vg_lookup_tables
