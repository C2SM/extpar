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

!> kind parameters are defined in MODULE data_parameters
USE mo_kind, ONLY: wp, &
                   i4, &
                   i8

!> abort_extpar defined in MODULE utilities_extpar
USE mo_utilities_extpar, ONLY: abort_extpar, free_un

USE mo_io_units,         ONLY: filename_max

USE mo_grid_structures,  ONLY: target_grid_def
USE mo_soil_tg_fields,   ONLY: soiltype_fao, soiltype_hwsd!, soiltype_deep

IMPLICIT NONE

PUBLIC :: calculate_soiltype, &
          read_namelists_extpar_HWSD_index


CONTAINS

  SUBROUTINE calculate_soiltype(tg,            &
     &                          ldeep_soil,    &
     &                          soiltype_fao,  &
     &                          soiltype_hwsd,  &
     &                          fr_sand,       &
     &                          fr_silt,       &
     &                          fr_clay,       &
     &                          fr_oc,         &
     &                          fr_bd,         &
     &                          fr_sand_deep,  &
     &                          fr_silt_deep,  &
     &                          fr_clay_deep,  &
     &                          fr_oc_deep,    &
     &                          fr_bd_deep     )
!     &                          soiltype_deep)

    TYPE(target_grid_def), INTENT(IN)      :: tg
    LOGICAL, INTENT(IN)                    :: ldeep_soil
    INTEGER (KIND=i4), INTENT(INOUT)       :: soiltype_hwsd(:,:,:)!(1:tg%ie,1:tg%je,1:tg%ke)
    INTEGER (KIND=i4), INTENT(INOUT)       :: soiltype_fao(:,:,:)!(1:tg%ie,1:tg%je,1:tg%ke)
   
    REAL(KIND=wp), INTENT(OUT)  :: fr_sand(1:tg%ie,1:tg%je,1:tg%ke) !< fraction sand due to HWSD
    REAL(KIND=wp), INTENT(OUT)  :: fr_silt(1:tg%ie,1:tg%je,1:tg%ke) !< fraction silt due to HWSD
    REAL(KIND=wp), INTENT(OUT)  :: fr_clay(1:tg%ie,1:tg%je,1:tg%ke) !< fraction clay due to HWSD
    REAL(KIND=wp), INTENT(OUT)  :: fr_oc(1:tg%ie,1:tg%je,1:tg%ke) !< fraction oc due to HWSD
    REAL(KIND=wp), INTENT(OUT)  :: fr_bd(1:tg%ie,1:tg%je,1:tg%ke) !< fraction bd due to HWSD

    REAL(KIND=wp)  :: texture(1:tg%ie,1:tg%je,1:tg%ke) ! Texture
    REAL(KIND=wp)  :: coarse(1:tg%ie,1:tg%je,1:tg%ke) ! Texture
    REAL(KIND=wp)  :: medium(1:tg%ie,1:tg%je,1:tg%ke) ! Texture
    REAL(KIND=wp)  :: fine(1:tg%ie,1:tg%je,1:tg%ke) ! Texture


!    INTEGER (KIND=i4), INTENT(INOUT), OPTIONAL:: soiltype_deep(:,:,:)!(1:tg%ie,1:tg%je,1:tg%ke)
    REAL(KIND=wp), INTENT(OUT),  OPTIONAL     :: fr_sand_deep(1:tg%ie,1:tg%je,1:tg%ke) !< fraction sand due to HWSD
    REAL(KIND=wp), INTENT(OUT), OPTIONAL      :: fr_silt_deep(1:tg%ie,1:tg%je,1:tg%ke) !< fraction silt due to HWSD
    REAL(KIND=wp), INTENT(OUT), OPTIONAL      :: fr_clay_deep(1:tg%ie,1:tg%je,1:tg%ke) !< fraction clay due to HWSD
    REAL(KIND=wp), INTENT(OUT), OPTIONAL      :: fr_oc_deep(1:tg%ie,1:tg%je,1:tg%ke) !< fraction oc due to HWSD
    REAL(KIND=wp), INTENT(OUT), OPTIONAL      :: fr_bd_deep(1:tg%ie,1:tg%je,1:tg%ke) !< fraction bd due to HWSD


    CHARACTER (len=filename_max) :: namelist_file !< filename with namelists for for EXTPAR settings

    ! HWSD idex files
    CHARACTER (len=filename_max) :: path_HWSD_index_files
    CHARACTER (len=filename_max) :: lookup_table_HWSD   
    CHARACTER (len=filename_max) :: HWSD_data  
    CHARACTER (len=filename_max) :: HWSD_data_deep     
    CHARACTER (len=filename_max) :: HWSD_data_extpar 

    CHARACTER (len=filename_max) :: path_lookup_table_HWSD
    CHARACTER (len=filename_max) :: path_HWSD_data
    CHARACTER (len=filename_max) :: path_HWSD_data_deep
    CHARACTER (len=filename_max) :: path_HWSD_data_extpar

    INTEGER (KIND=i4) :: ic, ic_deep

    INTEGER, PARAMETER:: n_soil=16102
    INTEGER           :: HWSD_ID(n_soil),stype
    REAL (KIND=wp)    :: HWSD_SU(n_soil), HWSD_TERRA(n_soil)

    INTEGER, PARAMETER:: n_soil_db=48148
    INTEGER           :: i_soil_db
    INTEGER           :: HWSD_ID_DB(n_soil_db),HWSD_SU_DB(n_soil_db)
    INTEGER           :: HWSD_SU_DB_S(n_soil_db)
    INTEGER           :: T_SAND(n_soil_db),T_SILT(n_soil_db),T_CLAY(n_soil_db)
    INTEGER           :: S_SAND(n_soil_db),S_SILT(n_soil_db),S_CLAY(n_soil_db)
    INTEGER (KIND=i4) :: nuin
    INTEGER (KIND=i4) :: i,j,k !<counters
    
    REAL(KIND=wp)     :: T_OC(n_soil_db),T_BD(n_soil_db),nfac
    REAL(KIND=wp)     :: S_OC(n_soil_db),S_BD(n_soil_db),nfac_deep
    CHARACTER*10      :: styp_out
    
    REAL(KIND=wp)     :: S,C,OM,D,zsandf,zsiltf,zclayf,zomf, zrog,zrocg, &
         &               zkw,avG,mvG, lnalpha,lnnm1,lnKs, zporv,         &
         &               zpwp, zfcap
    REAL(KIND=wp)     :: S_deep,C_deep,OM_deep,D_deep,                   &
         &               zsandf_deep,zsiltf_deep,zclayf_deep,            &
         &               zomf_deep, zrog_deep,zrocg_deep,                &
         &               zkw_deep,avG_deep,mvG_deep,                     &
         &               lnalpha_deep,lnnm1_deep,lnKs_deep, zporv_deep,  &
         &               zpwp_deep, zfcap_deep

    REAL(KIND=wp),     PARAMETER:: horizon_mid = 65.                    ! horizon mid-point of the subsoil in cm 
    REAL(KIND=wp),     PARAMETER:: minimum     = 0.1
    REAL(KIND=wp),     PARAMETER:: psi_fcap    = 0.1                    !bar
    REAL(KIND=wp),     PARAMETER:: psi_pwp     = 15.                    !bar
    REAL(KIND=wp),     PARAMETER:: zadp        = 0.01
    REAL(KIND=wp),     PARAMETER:: topsoil     = 1.
    REAL(KIND=wp),     PARAMETER:: subsoil     = 0.

      PRINT *,'Open HWSD index files...'

      namelist_file = 'INPUT_SOIL'

      CALL read_namelists_extpar_HWSD_index(namelist_file,        &
                                            path_HWSD_index_files,&
                                            lookup_table_HWSD,   &
                                            HWSD_data,            &
                                            HWSD_data_deep,       &
                                            HWSD_data_extpar)

! TODO: read namelist!!! mes

      nuin = free_un()
      path_lookup_table_HWSD = TRIM(path_HWSD_index_files)//TRIM(lookup_table_HWSD)
      PRINT*, 'path_lookup_table_HWSD: ', TRIM(path_lookup_table_HWSD)
      OPEN(nuin,file=TRIM(path_lookup_table_HWSD), status='old')
      READ(nuin,*) !header

      DO i=1,n_soil
        READ(nuin,*) HWSD_ID(i),HWSD_SU(i),HWSD_TERRA(i)
      END DO
      CLOSE(nuin)
      
      nuin = free_un()
      path_HWSD_data = TRIM(path_HWSD_index_files)//TRIM(HWSD_data)
      PRINT*, 'path_HWSD_data: ', TRIM(path_HWSD_data)
      OPEN(nuin,file=TRIM(path_HWSD_data), status='old')
      READ(nuin,*) !header
      
      DO i=1,n_soil_db
        READ(nuin,*) HWSD_SU_DB(i),T_SAND(i),T_SILT(i),T_CLAY(i),T_OC(i),T_BD(i)
      END DO
      CLOSE(nuin)

      IF (ldeep_soil) THEN
        nuin = free_un()
        path_HWSD_data_deep = TRIM(path_HWSD_index_files)//TRIM(HWSD_data_deep)
        PRINT*, 'path_HWSD_data_deep: ', TRIM(path_HWSD_data_deep)
        OPEN(nuin,file=TRIM(path_HWSD_data_deep), status='old')
        READ(nuin,*) !header
        DO i=1,n_soil_db
          READ(nuin,*) HWSD_SU_DB_S(i),S_SAND(i),S_SILT(i),S_CLAY(i),S_OC(i),S_BD(i)
        END DO
        CLOSE(nuin)
      ENDIF

      nuin = free_un()
      path_HWSD_data_extpar = TRIM(path_HWSD_index_files) // TRIM(HWSD_data_extpar)
      PRINT*, 'path_HWSD_data_extpar: ', TRIM(path_HWSD_data_extpar)
      OPEN(nuin,file=TRIM(path_HWSD_data_extpar), status='unknown')

       print*, 'MIN/MAX soiltype_HWSD : ', MINVAL(soiltype_hwsd), MAXVAL(soiltype_hwsd)

! Reset soil texture
       texture=0._wp
       coarse=0._wp
       fine=0._wp
       medium=0._wp


      DO k=1,tg%ke
        DO j=1,tg%je
          DO i=1,tg%ie
            
            ic=soiltype_hwsd(i,j,k)
            soiltype_hwsd(i,j,k)=HWSD_TERRA(ic)
!            IF (ldeep_soil) THEN
!              ic_deep = soiltype_deep(i,j,k)
!              soiltype_deep(i,j,k) = HWSD_TERRA(ic_deep)
!            ENDIF
            
            DO i_soil_db=1,n_soil_db
              IF(INT(HWSD_SU(ic))==HWSD_SU_DB(i_soil_db)) THEN
                
                nfac=(100.0-T_OC(i_soil_db))/100.
                
                fr_sand(i,j,k)=T_SAND(i_soil_db)
                fr_silt(i,j,k)=T_SILT(i_soil_db)
                fr_clay(i,j,k)=T_CLAY(i_soil_db)
                fr_oc(i,j,k)=T_OC(i_soil_db)
                fr_bd(i,j,k)=T_BD(i_soil_db)
                
                
              END IF
              IF (ldeep_soil) THEN
                IF(INT(HWSD_SU(ic))==HWSD_SU_DB_S(i_soil_db)) THEN
                
                  nfac_deep=(100.0-S_OC(i_soil_db))/100.
                
                  fr_sand_deep(i,j,k)=S_SAND(i_soil_db)
                  fr_silt_deep(i,j,k)=S_SILT(i_soil_db)
                  fr_clay_deep(i,j,k)=S_CLAY(i_soil_db)
                  fr_oc_deep(i,j,k)=S_OC(i_soil_db)
                  fr_bd_deep(i,j,k)=S_BD(i_soil_db)
                                
                END IF
              ENDIF
              
!              IF(HWSD_TERRA(ic)>0._wp.AND.HWSD_TERRA(ic).ne.12._wp) soiltype_fao(i,j,k)=HWSD_TERRA(ic) ! cities are in landuse data
              
              IF(HWSD_TERRA(ic)>0._wp.AND.HWSD_TERRA(ic).le.9._wp) soiltype_fao(i,j,k)=HWSD_TERRA(ic)
              IF(HWSD_TERRA(ic)>9._wp) soiltype_fao(i,j,k)=5 ! for undef soiltypes (<9 and 255) use loam

              IF(INT(HWSD_SU(ic))==HWSD_SU_DB(i_soil_db)) EXIT ! leave loop - only first entry will be considered
              IF (ldeep_soil) THEN
                IF(INT(HWSD_SU(ic_deep))==HWSD_SU_DB(i_soil_db)) EXIT ! leave loop - only first entry will be considered
              ENDIF


            END DO
                        
            
          ENDDO
        ENDDO
      ENDDO
      
      close(nuin)
      
      
    END SUBROUTINE calculate_soiltype

    !-----------------------------------------------------------------------------------------------------------------------!

    SUBROUTINE read_namelists_extpar_HWSD_index(namelist_file,        &
                                                path_HWSD_index_files,&
                                                lookup_table_HWSD,   &
                                                HWSD_data,            &
                                                HWSD_data_deep,       &
                                                HWSD_data_extpar)

      CHARACTER (len=1024), INTENT(IN) :: namelist_file !< filename with namelists for for EXTPAR settings

      ! HWSD idex files
      CHARACTER (len=1024) :: path_HWSD_index_files
      CHARACTER (len=1024) :: lookup_table_HWSD   
      CHARACTER (len=1024) :: HWSD_data   
      CHARACTER (len=1024) :: HWSD_data_deep    
      CHARACTER (len=1024) :: HWSD_data_extpar    

      !>Define the namelist group for soil raw data
      NAMELIST /HWSD_index_files/ path_HWSD_index_files, lookup_table_HWSD, HWSD_data, HWSD_data_deep, HWSD_data_extpar

      INTEGER           :: nuin !< unit number
      INTEGER (KIND=i4) :: ierr !< error flag


      nuin = free_un()  ! functioin free_un returns free Fortran unit number
      OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
      
      READ(nuin, NML=HWSD_index_files, IOSTAT=ierr)

      CLOSE(nuin)
  

 END SUBROUTINE read_namelists_extpar_HWSD_index


END MODULE mo_soil_consistency
