!+ Fortran Module with lookup-tables for the ecoclimap data
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_3         2011/04/19 Hermann Asensio
!  Initial release
! V2_0_3       2014/09/17 Burkhardt Rockel
!  Added use of directory information to access raw data files
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran Module with lookup-tables for the ecoclimap data
!> \author Hermann Asensio
!! ECOCLIMAP option gs_08.03.12  programm not operational yet
!! Description:
!! lookup-table from ECOCLIMAP 
MODULE mo_ecoclimap_lookup_tables

  USE mo_logging

  USE mo_kind,                  ONLY: wp, i4

  USE mo_io_units,              ONLY: filename_max

  USE mo_io_utilities,          ONLY: join_path 

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_ecoclimap_lookup_tables, & 
       &    get_name_ecoclimap_lookup_tables, & 
       &    ecoclimap_look_up, & 
       &    ecoclimap_legend, & 
       &    ecoclimap_value, & 
       &    nclass_ecoclimap, & 
       &    ilookup_table_ecoclimap, & 
       &    i_extpar_lookup_table , & 
       &    name_lookup_table_ecoclimap, & 
       &    z012_lt_ecoclimap, lnz012_lt_ecoclimap, plc12_lt_ecoclimap, & 
       &    lai12_lt_ecoclimap, rd_lt_ecoclimap, emiss12_lt_ecoclimap, rs_min_lt_ecoclimap, & 
       &    forest_type_ecoclimap

  INTEGER (KIND=i4), PARAMETER :: nclass_ecoclimap = 243, &  !< ecoclimap has 243 classes for the land use description
       &                          ntime_ecoclimap = 12, & 
       &                          i_extpar_lookup_table = 1 !< lookup_table for ecoclimap land use classes


  INTEGER (KIND=i4)            :: ilookup_table_ecoclimap !< integer switch to choose a lookup table
  CHARACTER (LEN=filename_max) :: name_lookup_table_ecoclimap !< name of lookup table

  !< lookup table landuse class to roughness length [m]
  REAL (KIND=wp)               :: z012_lt_ecoclimap(ntime_ecoclimap, nclass_ecoclimap), & 
                     &            lnz012_lt_ecoclimap(ntime_ecoclimap, nclass_ecoclimap), &
                !< lookup table landuse class to minimal plant cover
                     &            plc12_lt_ecoclimap(ntime_ecoclimap, nclass_ecoclimap), &
                !< lookup table landuse class to minimal leaf area index
                     &            lai12_lt_ecoclimap(ntime_ecoclimap, nclass_ecoclimap), &
                !< lookup table landuse class to root depth [m]
                     &            rd_lt_ecoclimap(nclass_ecoclimap) , &
                !< lookup table landuse class to surface thermal emissivity
                     &            emiss12_lt_ecoclimap(ntime_ecoclimap, nclass_ecoclimap), &
                !< lookup table landuse class to minimal stomata resistance
                     &            rs_min_lt_ecoclimap(nclass_ecoclimap), &
                !< lookup table landuse class to forest type
                     &            forest_type_ecoclimap(3,nclass_ecoclimap)  

  !> legend of the ecoclimap vegetation classes
  CHARACTER(len=45)            :: ecoclimap_legend(nclass_ecoclimap)  ! No. 

  !> values of the ecoclimap landuse classes
  INTEGER(KIND=i4)             :: ecoclimap_value(nclass_ecoclimap)  ! No.

  CONTAINS

  !> define lookup table for ecoclimap landuse classes
  SUBROUTINE init_ecoclimap_lookup_tables(raw_data_lu_path, & !_br 17.09.14
       &                                  nclass_ecoclimap, &      !_br 17.09.14
       &                                  z012_lt_ecoclimap,           &
       &                                  lnz012_lt_ecoclimap,       &
       &                                  plc12_lt_ecoclimap,        &
       &                                  lai12_lt_ecoclimap,        &
       &                                  rd_lt_ecoclimap,          &
       &                                  emiss12_lt_ecoclimap,       &
       &                                  rs_min_lt_ecoclimap,       &
       &                                  forest_type_ecoclimap)

    CHARACTER (LEN=filename_max):: raw_data_lu_path        !< path to raw data !_br 17.09.14
    INTEGER, INTENT(IN)         :: nclass_ecoclimap !< ecoclimap has 243 classes for the land use description
    !< lookup table landuse class to roughness length [m]
    REAL (KIND=wp), INTENT(OUT) :: z012_lt_ecoclimap(ntime_ecoclimap, nclass_ecoclimap), &
    !< corresponding natural logarithm of z0c_extpar_o
         &                         lnz012_lt_ecoclimap(ntime_ecoclimap,nclass_ecoclimap), &
    !< lookup table landuse class to minimal plant cover
         &                         plc12_lt_ecoclimap(ntime_ecoclimap, nclass_ecoclimap), &
    !< lookup table landuse class to minimal leaf area index
         &                         lai12_lt_ecoclimap(ntime_ecoclimap, nclass_ecoclimap), &
    !< lookup table landuse class to root depth [m]
         &                         rd_lt_ecoclimap(nclass_ecoclimap), &
    !< lookup table landuse class to surface thermal emissivity
         &                         emiss12_lt_ecoclimap(ntime_ecoclimap,nclass_ecoclimap), &    
    !< lookup table landuse class to minimal stomata resistance
         &                         rs_min_lt_ecoclimap(nclass_ecoclimap), &   
    !< lookup table landuse class to forest type
         &                         forest_type_ecoclimap(3, nclass_ecoclimap)  
    
    ! local variable
    INTEGER(KIND=i4)            :: i,j,k,io_error !< counter
    REAL(KIND=wp)               :: arg
    CHARACTER (LEN=10)          :: dum,dum1

    !READ LOOK UP TABLES
    CALL logging%info('READ LOOK UP TABLES ECOCLIMAP')
    OPEN (UNIT=10, FILE=join_path(raw_data_lu_path,'ecoclimap_lookup.TAB'),ACTION='read', STATUS='old', iostat=io_error) !_br 17.09.14
      IF ( io_error == 0 ) THEN    
        READ (10, *) dum
!LAI    
        DO j = 1, nclass_ecoclimap
          READ (10, *) dum, (lai12_lt_ecoclimap(i,j), i=1,ntime_ecoclimap),rd_lt_ecoclimap(j),dum1,dum1  
        END DO
!PLCR
        READ (10, *) dum
        READ (10, *) dum
          DO j = 1, nclass_ecoclimap
            READ (10, *) dum, (plc12_lt_ecoclimap(i,j), i=1,12)
          END DO
!Z0
        READ (10, *) dum
        READ (10, *) dum
          DO j = 1, nclass_ecoclimap
             READ (10, *) dum, (z012_lt_ecoclimap(i,j), i=1,12)
          END DO
!VEG
        READ (10, *) dum
        READ (10, *) dum
          DO j = 1, nclass_ecoclimap
             READ (10, *) dum, (forest_type_ecoclimap(i,j), i=1,3)
          END DO
!RSM
        READ (10, *) dum
        READ (10, *) dum
        DO j = 1, nclass_ecoclimap
           READ (10, *) dum,dum,dum, rs_min_lt_ecoclimap(j)
        END DO
!EMI
        READ (10, *) dum
        READ (10, *) dum
          DO j = 1, nclass_ecoclimap
            READ (10, *) dum, (emiss12_lt_ecoclimap(i,j), i=1,12)
          END DO

      ELSE
        WRITE(message_text,*)'Could not open ecoclimap_lookup.TAB ', io_error
        CALL  logging%error(message_text,__FILE__,__LINE__)
      END IF
      CLOSE (unit=10)

      nz012_lt_ecoclimap = 0.
      DO k = 1, 12
        DO i=1,nclass_ecoclimap
          IF (z012_lt_ecoclimap(k,i) > 0.) THEN
            arg = z012_lt_ecoclimap(k,i)
            lnz012_lt_ecoclimap(k,i) = LOG(arg)
          ENDIF
        ENDDO
      END DO
    END SUBROUTINE init_ecoclimap_lookup_tables

  !> define  name of lookup table for ecoclimap
  SUBROUTINE get_name_ecoclimap_lookup_tables(ilookup_table_ecoclimap, name_lookup_table_ecoclimap)
    INTEGER, INTENT(IN) :: ilookup_table_ecoclimap  !< integer switch to choose a lookup table
    CHARACTER (LEN=filename_max), INTENT(OUT) :: name_lookup_table_ecoclimap !< name of lookup table
    ! local variable
      SELECT CASE (ilookup_table_ecoclimap)
        CASE(i_extpar_lookup_table)
           name_lookup_table_ecoclimap='ECOCLIMAP '
        CASE DEFAULT
           name_lookup_table_ecoclimap='ECOCLIMAP'
      END SELECT

  END  SUBROUTINE get_name_ecoclimap_lookup_tables



   !> assign the ecoclimap land use classes to some characteristic (more or less) physical parameters
  SUBROUTINE ecoclimap_look_up(lu, &
    &      nclass_ecoclimap, &
    &      plc12_lt_ecoclimap,        &
    &      lai12_lt_ecoclimap,        &
    &      rd_lt_ecoclimap,            &
    &      emiss12_lt_ecoclimap,       &
    &      rs_min_lt_ecoclimap,        &
    &      forest_type_ecoclimap,      &
    &      pland,          &
    &      pice,           &
    &      plnz0,          &
    &      proot,          &
    &      p12,            &
    &      plai12,         &
    &      purb,           &
    &      pfor_d,         &
    &      pfor_e,         &
    &      pemissivity12,  &
    &      prs_min,        &
    &      k_error)

  INTEGER, INTENT(IN) :: lu             !< land use class
  INTEGER, INTENT(IN) :: nclass_ecoclimap !< ecoclimap has 23 classes for the land use description
  !< lookup table landuse class to minimal plant cover
  REAL (KIND=wp), INTENT(IN) :: plc12_lt_ecoclimap(ntime_ecoclimap, nclass_ecoclimap)  
  !< lookup table landuse class to minimal leaf area index
  REAL (KIND=wp), INTENT(IN) :: lai12_lt_ecoclimap(ntime_ecoclimap, nclass_ecoclimap)  
  !< lookup table landuse class to root depth [m]
  REAL (KIND=wp), INTENT(IN) :: rd_lt_ecoclimap(nclass_ecoclimap)      
  !< lookup table landuse class to surface thermal emissivity
  REAL (KIND=wp), INTENT(IN) :: emiss12_lt_ecoclimap(ntime_ecoclimap,nclass_ecoclimap)   
  !< lookup table landuse class to minimal stomata resistance
  REAL (KIND=wp), INTENT(IN) :: rs_min_lt_ecoclimap(nclass_ecoclimap)  
  !< lookup table landuse class to forest type
  REAL (KIND=wp), INTENT(IN) :: forest_type_ecoclimap(3, nclass_ecoclimap)  


  REAL (KIND=wp), INTENT(OUT) :: pland          !< land cover                      (-)
  REAL (KIND=wp), INTENT(OUT) :: pice           !< ice fraction                    (-)
  REAL (KIND=wp), INTENT(OUT) :: plnz0(ntime_ecoclimap)          !< logarithm of roughness length   (m)
  REAL (KIND=wp), INTENT(OUT) :: proot          !< root depth                      (m)
  REAL (KIND=wp), INTENT(OUT) :: p12(ntime_ecoclimap)            !<  plant cover             (-)
  REAL (KIND=wp), INTENT(OUT) :: plai12(ntime_ecoclimap)   !<  leaf area index         (m**2/m**2)
  REAL (KIND=wp), INTENT(OUT) :: purb           !< urbanisation                    (-)
  REAL (KIND=wp), INTENT(OUT) :: pfor_d         !< deciduous forest                (-)
  REAL (KIND=wp), INTENT(OUT) :: pfor_e         !< evergreen forest                (-)
  REAL (KIND=wp), INTENT(OUT) :: pemissivity12(ntime_ecoclimap)   !< surface thermal emissivity      (-)
  REAL (KIND=wp), INTENT(OUT) :: prs_min        !< minimum stomata resistance      (s/m)

  INTEGER, INTENT(OUT)        :: k_error     !< error return code

  ! local variables
  INTEGER :: k !<  time stesp

       ! Test for true land points

           pice    = 0.0
           purb    = 0.0

          IF (lu>=4 .AND. lu<=244) THEN
            k_error     = 0
            pland       = 1.0
            DO k = 1, ntime_ecoclimap
              plnz0(k)         = lnz012_lt_ecoclimap(k,lu)
              p12(k)           = plc12_lt_ecoclimap(k,lu)
              plai12(k)        = lai12_lt_ecoclimap(k, lu)
              pemissivity12(k) = emiss12_lt_ecoclimap(k,lu)
            END DO
            proot       = rd_lt_ecoclimap(lu)

            prs_min     = rs_min_lt_ecoclimap(lu)
            pfor_d  = forest_type_ecoclimap(1,lu) 
            pfor_e  = forest_type_ecoclimap(2,lu) + forest_type_ecoclimap(3,lu)

            IF (lu==7              ) purb   = 0.60  ! artificial surface
            IF (lu==151            ) purb   = 0.90  ! artificial surface
            IF (lu==152            ) purb   = 0.60  ! artificial surface
            IF (lu==153            ) purb   = 0.60  ! artificial surface
            IF (lu==154            ) purb   = 0.60  ! artificial surface
            IF (lu==155            ) purb   = 0.90  ! artificial surface
            IF (lu==156            ) purb   = 0.90  ! artificial surface
            IF (lu==157            ) purb   = 0.90  ! artificial surface
            IF (lu==158            ) purb   = 0.30  ! artificial surface
            IF (lu==159            ) purb   = 0.10  ! artificial surface
            IF (lu==160            ) purb   = 0.10  ! artificial surface
            IF (lu==161            ) purb   = 0.20  ! artificial surface

            IF (lu==6              ) pice   = 1.0  ! ice or snow pixel

         ELSE IF (lu > 0 .AND. lu< 4) THEN ! water
            k_error     = 0
            pland       = 0.0
         ELSE
            k_error     = 1  ! not a valid land use class
            pland       = 0.0
         END IF
  END  SUBROUTINE ecoclimap_look_up

END MODULE mo_ecoclimap_lookup_tables
