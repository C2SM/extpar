!+ Fortran module with data fields for soil data
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V2_0         24/4/2013 Juergen Helmert / Martina Messmer
!  Initial release
!  The soiltype is converted from the HWSD world code to the TERRA code.
! V1_14        2014-07-18 Juergen Helmert
!  Combined COSMO Release
! V2_1         2015-01-12 Juergen Helmert
!  Bugfix correction covers CSCS SVN r5907-r6359
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module that performs calculation of the soil type and bulk density
!> \author Juergen Helmert / Martina Messmer
MODULE mo_soil_consistency

  USE mo_logging
  USE mo_kind,                    ONLY: wp, i4

  USE mo_utilities_extpar,        ONLY: free_un
                                  
  USE mo_io_units,                ONLY: filename_max
                                  
  USE mo_grid_structures,         ONLY: target_grid_def
  USE mo_soil_tg_fields,          ONLY: soiltype_fao, soiltype_hwsd!, soiltype_deep

  IMPLICIT NONE

  PUBLIC :: calculate_soiltype, &
       &    read_namelists_extpar_HWSD_index


  CONTAINS

  SUBROUTINE calculate_soiltype(tg,            &
        &                       ldeep_soil,    &
        &                       soiltype_fao,  &
        &                       soiltype_hwsd,  &
        &                       fr_sand,       &
        &                       fr_silt,       &
        &                       fr_clay,       &
        &                       fr_oc,         &
        &                       fr_bd,         &
        &                       fr_sand_deep,  &
        &                       fr_silt_deep,  &
        &                       fr_clay_deep,  &
        &                       fr_oc_deep,    &
        &                       fr_bd_deep     )

    TYPE(target_grid_def), INTENT(IN)    :: tg
    LOGICAL, INTENT(IN)                  :: ldeep_soil
    INTEGER (KIND=i4), INTENT(INOUT)     :: soiltype_hwsd(:,:,:), &!(1:tg%ie,1:tg%je,1:tg%ke)
         &                                  soiltype_fao(:,:,:)!(1:tg%ie,1:tg%je,1:tg%ke)
   
    REAL(KIND=wp), INTENT(OUT)           :: fr_sand(1:tg%ie,1:tg%je,1:tg%ke), & !< fraction sand due to HWSD
         &                                  fr_silt(1:tg%ie,1:tg%je,1:tg%ke), & !< fraction silt due to HWSD
         &                                  fr_clay(1:tg%ie,1:tg%je,1:tg%ke), & !< fraction clay due to HWSD
         &                                  fr_oc(1:tg%ie,1:tg%je,1:tg%ke), & !< fraction oc due to HWSD
         &                                  fr_bd(1:tg%ie,1:tg%je,1:tg%ke) !< fraction bd due to HWSD

    REAL(KIND=wp), INTENT(OUT), OPTIONAL :: fr_sand_deep(1:tg%ie,1:tg%je,1:tg%ke), & !< fraction sand due to HWSD
         &                                  fr_silt_deep(1:tg%ie,1:tg%je,1:tg%ke), & !< fraction silt due to HWSD
         &                                  fr_clay_deep(1:tg%ie,1:tg%je,1:tg%ke), & !< fraction clay due to HWSD
         &                                  fr_oc_deep(1:tg%ie,1:tg%je,1:tg%ke), & !< fraction oc due to HWSD
         &                                  fr_bd_deep(1:tg%ie,1:tg%je,1:tg%ke) !< fraction bd due to HWSD

    CHARACTER (len=filename_max)         :: namelist_file, & !< filename with namelists for for EXTPAR settings
    ! HWSD idex files                    
         &                                  path_HWSD_index_files, &
         &                                  lookup_table_HWSD, &   
         &                                  HWSD_data, &  
         &                                  HWSD_data_deep, &     
         &                                  HWSD_data_extpar, & 
         &                                  path_lookup_table_HWSD, &
         &                                  path_HWSD_data, &
         &                                  path_HWSD_data_deep, &
         &                                  path_HWSD_data_extpar

    INTEGER, PARAMETER                  :: n_soil=16102, &
         &                                 n_soil_db=48148

    INTEGER (KIND=i4)                    :: ic, &
         &                                  HWSD_ID(n_soil), &
         &                                  i_soil_db, &
         &                                  HWSD_SU_DB(n_soil_db), &
         &                                  HWSD_SU_DB_S(n_soil_db), &
         &                                  T_SAND(n_soil_db),T_SILT(n_soil_db),T_CLAY(n_soil_db), &
         &                                  S_SAND(n_soil_db),S_SILT(n_soil_db),S_CLAY(n_soil_db), &
         &                                  nuin, &
         &                                  i,j,k !<counters
    
    REAL(KIND=wp)                       :: T_OC(n_soil_db),T_BD(n_soil_db), &
         &                                 HWSD_SU(n_soil), HWSD_TERRA(n_soil), &
         &                                 S_OC(n_soil_db),S_BD(n_soil_db)

    REAL(KIND=wp), PARAMETER            :: horizon_mid = 65., &                    ! horizon mid-point of the subsoil in cm 
         &                                 minimum     = 0.1, &
         &                                 psi_fcap    = 0.1, &                    !bar
         &                                 psi_pwp     = 15., &                    !bar
         &                                 zadp        = 0.01, &
         &                                 topsoil     = 1., &
         &                                 subsoil     = 0.

    namelist_file = 'INPUT_SOIL'

    CALL read_namelists_extpar_HWSD_index(namelist_file,        &
                                          path_HWSD_index_files,&
                                          lookup_table_HWSD,   &
                                          HWSD_data,            &
                                          HWSD_data_deep,       &
                                          HWSD_data_extpar)

    nuin = free_un()
    path_lookup_table_HWSD = TRIM(path_HWSD_index_files)//TRIM(lookup_table_HWSD)
    OPEN(nuin,file=TRIM(path_lookup_table_HWSD), status='old')
    READ(nuin,*) !header

    DO i=1,n_soil
      READ(nuin,*) HWSD_ID(i),HWSD_SU(i),HWSD_TERRA(i)
    END DO
    CLOSE(nuin)
    
    nuin = free_un()
    path_HWSD_data = TRIM(path_HWSD_index_files)//TRIM(HWSD_data)
    OPEN(nuin,file=TRIM(path_HWSD_data), status='old')
    READ(nuin,*) !header
    
    DO i=1,n_soil_db
      READ(nuin,*) HWSD_SU_DB(i),T_SAND(i),T_SILT(i),T_CLAY(i),T_OC(i),T_BD(i)
    END DO
    CLOSE(nuin)

    IF (ldeep_soil) THEN
      nuin = free_un()
      path_HWSD_data_deep = TRIM(path_HWSD_index_files)//TRIM(HWSD_data_deep)
      OPEN(nuin,file=TRIM(path_HWSD_data_deep), status='old')
      READ(nuin,*) !header
      DO i=1,n_soil_db
        READ(nuin,*) HWSD_SU_DB_S(i),S_SAND(i),S_SILT(i),S_CLAY(i),S_OC(i),S_BD(i)
      END DO
      CLOSE(nuin)
    ENDIF

    nuin = free_un()
    path_HWSD_data_extpar = TRIM(path_HWSD_index_files) // TRIM(HWSD_data_extpar)
    OPEN(nuin,file=TRIM(path_HWSD_data_extpar), status='unknown')

    DO k=1,tg%ke
      DO j=1,tg%je
        DO i=1,tg%ie
          ic=soiltype_hwsd(i,j,k)
          soiltype_hwsd(i,j,k)=INT(HWSD_TERRA(ic))
          DO i_soil_db=1,n_soil_db
            IF(INT(HWSD_SU(ic))==HWSD_SU_DB(i_soil_db)) THEN
              fr_sand(i,j,k)=T_SAND(i_soil_db)
              fr_silt(i,j,k)=T_SILT(i_soil_db)
              fr_clay(i,j,k)=T_CLAY(i_soil_db)
              fr_oc(i,j,k)=T_OC(i_soil_db)
              fr_bd(i,j,k)=T_BD(i_soil_db)
            END IF
            IF (ldeep_soil) THEN
              IF(INT(HWSD_SU(ic))==HWSD_SU_DB_S(i_soil_db)) THEN
              
                fr_sand_deep(i,j,k)=S_SAND(i_soil_db)
                fr_silt_deep(i,j,k)=S_SILT(i_soil_db)
                fr_clay_deep(i,j,k)=S_CLAY(i_soil_db)
                fr_oc_deep(i,j,k)=S_OC(i_soil_db)
                fr_bd_deep(i,j,k)=S_BD(i_soil_db)
                              
              END IF
            ENDIF
            
            IF(HWSD_TERRA(ic)>0._wp.AND.HWSD_TERRA(ic).le.9._wp) soiltype_fao(i,j,k)=INT(HWSD_TERRA(ic))
            IF(HWSD_TERRA(ic)>9._wp) soiltype_fao(i,j,k)=5 ! for undef soiltypes (<9 and 255) use loam

            IF(INT(HWSD_SU(ic))==HWSD_SU_DB(i_soil_db)) EXIT ! leave loop - only first entry will be considered
          END DO
        ENDDO
      ENDDO
    ENDDO
    
    close(nuin)
      
  END SUBROUTINE calculate_soiltype

  SUBROUTINE read_namelists_extpar_HWSD_index(namelist_file,        &
                                                path_HWSD_index_files,&
                                                lookup_table_HWSD,   &
                                                HWSD_data,            &
                                                HWSD_data_deep,       &
                                                HWSD_data_extpar)

    CHARACTER (len=1024), INTENT(IN) :: namelist_file !< filename with namelists for for EXTPAR settings

    ! HWSD idex files
    CHARACTER (len=1024)                    :: path_HWSD_index_files, &
         &                                     lookup_table_HWSD, &   
         &                                     HWSD_data, &   
         &                                     HWSD_data_deep, &    
         &                                     HWSD_data_extpar    
                                            
    INTEGER (KIND=i4)                       :: ierr, nuin

    !>Define the namelist group for soil raw data
    NAMELIST /HWSD_index_files/ path_HWSD_index_files, lookup_table_HWSD, HWSD_data, HWSD_data_deep, HWSD_data_extpar

    nuin = free_un()  ! functioin free_un returns free Fortran unit number
    OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF
    
    READ(nuin, NML=HWSD_index_files, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist HWSD_index_files',__FILE__, __LINE__) 
    ENDIF

    CLOSE(nuin)

 END SUBROUTINE read_namelists_extpar_HWSD_index

END MODULE mo_soil_consistency
