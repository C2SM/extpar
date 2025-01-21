MODULE mo_python_routines

  USE mo_logging
  USE mo_kind,                  ONLY:  i4, wp

  USE mo_utilities_extpar,      ONLY: free_un ! function to get free unit number

  USE mo_io_utilities,          ONLY: check_netcdf

  USE mo_io_units,              ONLY: filename_max

  USE netcdf,                     ONLY:   &
    &                                   nf90_open,              &
    &                                   nf90_close,             &
    &                                   nf90_nowrite

  USE mo_python_tg_fields,        ONLY: alb_interpol
  USE mo_soil_tg_fields,          ONLY: soiltype_fao
  USE mo_target_grid_data,        ONLY: tg, &
    &                                   lon_geo,lat_geo

  USE mo_bilinterpol,             ONLY: calc_weight_bilinear_interpol, &
    &                                   calc_value_bilinear_interpol

  USE mo_python_data,             ONLY: zalso, &
    &                                   max_tiles_isa

  USE mo_icon_grid_data,          ONLY: icon_grid_region

  IMPLICIT NONE

  PUBLIC :: &
  ! emiss
       &    read_namelists_extpar_emiss, &
  ! cru
       &    read_namelists_extpar_t_clim, &
  ! ndvi
       &    read_namelists_extpar_ndvi, &
  ! edgar
       &    read_namelists_extpar_edgar, &
  ! albedo
       &    read_namelists_extpar_alb, &
       &    open_netcdf_ALB_data, &
       &    const_check_interpol_alb, &
  ! era
       &    read_namelists_extpar_era, &
  ! ahf
       &    read_namelists_extpar_ahf, &
  ! isa
       &    read_namelists_extpar_isa, &
  ! HHS      
       &    read_namelists_extpar_hhs, &
       &    read_namelists_extpar_hhs_ormc, &
       &    read_namelists_extpar_hhs_alfa, &
       &    read_namelists_extpar_hhs_critw, &
       &    read_namelists_extpar_hhs_fieldc, &
       &    read_namelists_extpar_hhs_n, &
       &    read_namelists_extpar_hhs_satf, &
       &    read_namelists_extpar_hhs_stc, &
       &    read_namelists_extpar_hhs_wcav, &
       &    read_namelists_extpar_hhs_wcpf2, &
       &    read_namelists_extpar_hhs_wcpf3, &              
       &    read_namelists_extpar_hhs_wcpf42, &
       &    read_namelists_extpar_hhs_wcres, &
       &    read_namelists_extpar_hhs_wcsat                     
       
       CONTAINS

  !---------------------------------------------------------------------------
  !> subroutine to read namelist for EMISS data settings for EXTPAR 
  SUBROUTINE read_namelists_extpar_emiss(namelist_file, &
       &                                 raw_data_emiss_path, &
       &                                 raw_data_emiss_filename, &
       &                                 emiss_buffer_file, &
       &                                 emiss_output_file)

    CHARACTER (len=*), INTENT(IN)            :: namelist_file !< filename with namelists for for EXTPAR settings

    CHARACTER (len=filename_max),INTENT(OUT) :: raw_data_emiss_path, &         !< path to raw data
         &                                      raw_data_emiss_filename, &  !< filename EMISS raw data
         &                                      emiss_buffer_file, &  !< name for EMISS buffer file
         &                                      emiss_output_file !< name for EMISS output file

    INTEGER(KIND=i4)                         :: nuin, ierr

    !> namelist with filenames for EMISS data input
    NAMELIST /emiss_raw_data/ raw_data_emiss_path, raw_data_emiss_filename
    !> namelist with filenames for EMISS data output
    NAMELIST /emiss_io_extpar/ emiss_buffer_file, emiss_output_file

    nuin = free_un()  ! functioin free_un returns free Fortran unit number
    OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=emiss_raw_data, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot read in namelist emiss_raw_data - reason: ', ierr
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=emiss_io_extpar, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot read in namelist emiss_io_extpar - reason: ', ierr
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF

    CLOSE(nuin)

  END SUBROUTINE read_namelists_extpar_emiss

  !---------------------------------------------------------------------------
  !> subroutine to read namelist for t_clim data settings for EXTPAR 
  SUBROUTINE read_namelists_extpar_t_clim(namelist_file,            &
       &                                  it_cl_type,               &
       &                                  raw_data_t_clim_path,     &
       &                                  raw_data_t_clim_filename, &
       &                                  t_clim_buffer_file,       &
       &                                  t_clim_output_file)

    CHARACTER (len=*), INTENT(IN)             :: namelist_file !< filename with namelists for for EXTPAR settings
    INTEGER (KIND=i4), INTENT(OUT)            :: it_cl_type    !< integer switch to choose a land use raw data set

    CHARACTER (len=filename_max), INTENT(OUT) :: raw_data_t_clim_path, &        !< path to raw data
         &                                       raw_data_t_clim_filename, &    !< filename temperature climatology raw data
         &                                       t_clim_buffer_file, & !< name for temperature climatology buffer
         &                                       t_clim_output_file !< name for temperature climatology output file
    

    INTEGER (KIND=i4)                         :: nuin, ierr

    !> namelist with filename for temperature climatlogy data output
    NAMELIST /t_clim_raw_data/ raw_data_t_clim_path, raw_data_t_clim_filename, it_cl_type

    !> namelist with filename for temperature climatlogy data output
    NAMELIST /t_clim_io_extpar/ t_clim_buffer_file, t_clim_output_file

    it_cl_type = -1
    
    raw_data_t_clim_path = ''
    raw_data_t_clim_filename = ''

    t_clim_buffer_file = ''
    t_clim_output_file = ''

    OPEN(NEWUNIT=nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('CRU namelist open error ', __FILE__, __LINE__)
    ENDIF
    READ(nuin, NML=t_clim_raw_data, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot read in namelist t_clim_raw_data - reason: ', ierr
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF
    READ(nuin, NML=t_clim_io_extpar, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot read in namelist t_clim_io_extpar - reason: ', ierr
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF
    CLOSE(nuin)
    
  END SUBROUTINE read_namelists_extpar_t_clim

  !> subroutine to read namelist for NDVI data settings for EXTPAR 
  SUBROUTINE read_namelists_extpar_ndvi(namelist_file, &
       &                                raw_data_ndvi_path, &
       &                                raw_data_ndvi_filename, &
       &                                ndvi_buffer_file, &
       &                                ndvi_output_file)



    CHARACTER (len=*), INTENT(IN) :: namelist_file !< filename with namelists for for EXTPAR settings

    CHARACTER (len=filename_max)             :: raw_data_ndvi_path, &        !< path to raw data
         &                                      raw_data_ndvi_filename, & !< filename NDVI raw data
         &                                      ndvi_buffer_file, & !< name for NDVI buffer file
         &                                      ndvi_output_file !< name for NDVI output file
       
    INTEGER (KIND=i4)                        :: ierr, nuin

    !> namelist with filenames for NDVI data input
    NAMELIST /ndvi_raw_data/ raw_data_ndvi_path, raw_data_ndvi_filename
    !> namelist with filenames for NDVI data output
    NAMELIST /ndvi_io_extpar/ ndvi_buffer_file, ndvi_output_file

    nuin = free_un()  ! functioin free_un returns free Fortran unit number
    OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=ndvi_raw_data, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot read in namelist ndvi_raw_data - reason: ', ierr
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=ndvi_io_extpar, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot read in namelist ndvi_io_extpar - reason: ', ierr
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF

    CLOSE(nuin)

  END SUBROUTINE read_namelists_extpar_ndvi

  !> subroutine to read namelist for EDGAR data settings for EXTPAR 
  SUBROUTINE read_namelists_extpar_edgar(namelist_file, edgar_buffer_file)

    CHARACTER (len=*), INTENT(IN)             :: namelist_file                      !< filename with namelists
    CHARACTER (len=filename_max), INTENT(OUT) :: edgar_buffer_file                  !< name for EDGAR buffer file

    INTEGER (KIND=i4)                         :: ierr, nuin

    !> namelist with filenames for EDGAR data output
    NAMELIST /edgar_io_extpar/ edgar_buffer_file

    nuin = free_un()  ! function free_un returns free Fortran unit number
    OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=edgar_io_extpar, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot read in namelist edgar_io_extpar - reason: ', ierr
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF

    CLOSE(nuin)

  END SUBROUTINE read_namelists_extpar_edgar

  !---------------------------------------------------------------------------
  !> subroutine to read namelist including albedo data settings for EXTPAR 
  SUBROUTINE read_namelists_extpar_alb(namelist_file, &
      &                                  raw_data_alb_path, &
      &                                  raw_data_alb_filename, &
      &                                  raw_data_alnid_filename, &
      &                                  raw_data_aluvd_filename, &
      &                                  ialb_type,         &
      &                                  alb_buffer_file, &
      &                                  alb_output_file)
    
    CHARACTER (len=*), INTENT(IN)  :: namelist_file !< filename with namelists for for EXTPAR settings

    CHARACTER (len=filename_max)   :: raw_data_alb_path, &        !< path to raw data
      &                               raw_data_alb_filename, & !< filenames albedo raw data
      &                               raw_data_alnid_filename, &
      &                               raw_data_aluvd_filename, &
      &                               alb_buffer_file, & !< name for albedo buffer file
      &                               alb_output_file !< name for albedo output file

    INTEGER (KIND=i4)              :: ialb_type, nuin, ierr

  !> namelist with filenames for albedo data input
    NAMELIST /alb_raw_data/ raw_data_alb_path, raw_data_alb_filename, ialb_type
    NAMELIST /alnid_raw_data/ raw_data_alb_path, raw_data_alnid_filename
    NAMELIST /aluvd_raw_data/ raw_data_alb_path, raw_data_aluvd_filename
  !> namelist with filenames for albedo data output
    NAMELIST /alb_io_extpar/ alb_buffer_file, alb_output_file


    nuin = free_un()  ! functioin free_un returns free Fortran unit number

    OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=alb_raw_data, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot read in namelist alb_raw_data - reason: ', ierr
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF

    IF ((ialb_type < 1).OR.(ialb_type > 3)) THEN
      WRITE(message_text,*) 'ialb_type must be in the range 1-3. It is now:',ialb_type
      CALL logging%error(message_text,__FILE__,__LINE__)
    ENDIF

    IF ((ialb_type /= 2).AND.(ialb_type /= 3)) THEN
      READ(nuin, NML=alnid_raw_data, IOSTAT=ierr)
      IF (ierr /= 0) THEN
        WRITE(message_text,*)'Cannot read in namelist alnid_raw_data - reason: ', ierr
        CALL logging%error(message_text,__FILE__, __LINE__) 
      ENDIF

      READ(nuin, NML=aluvd_raw_data, IOSTAT=ierr)
      IF (ierr /= 0) THEN
        WRITE(message_text,*)'Cannot read in namelist aluvd_raw_data - reason: ', ierr
        CALL logging%error(message_text,__FILE__, __LINE__) 
      ENDIF

    ENDIF
    READ(nuin, NML=alb_io_extpar, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot read in namelist alb_io_extpar - reason: ', ierr
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF

    CLOSE(nuin)
    
  END SUBROUTINE read_namelists_extpar_alb

  !> subroutine to read namelist for ERA data 
  SUBROUTINE read_namelists_extpar_era(namelist_file, &
       &                               iera_type, &
       &                               era_buffer_file)

    CHARACTER (len=*), INTENT(IN) :: namelist_file !< filename with namelists for for EXTPAR settings

    CHARACTER (len=filename_max)             :: era_buffer_file !< name for ERA buffer file
       
    INTEGER (KIND=i4)                        :: iera_type, ierr, nuin

    !> namelist with filenames for ERA data output
    NAMELIST /era_raw_data/ iera_type
    NAMELIST /era_io_extpar/ era_buffer_file

    nuin = free_un()  ! functioin free_un returns free Fortran unit number
    OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=era_raw_data, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist era_raw_data',__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=era_io_extpar, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist era_io_extpar',__FILE__, __LINE__) 
    ENDIF

    CLOSE(nuin)

  END SUBROUTINE read_namelists_extpar_era

  !> subroutine to read namelist for AHF data settings for EXTPAR 
  SUBROUTINE read_namelists_extpar_ahf(namelist_file, &
                                      iahf_type,    & !_br 14.04.16
                                      raw_data_ahf_path, &
                                      raw_data_ahf_filename, &
                                      ahf_buffer_file)
  
    CHARACTER (len=1024), INTENT(IN) :: namelist_file !< filename with namelists for for EXTPAR settings

    CHARACTER (len=1024)            :: raw_data_ahf_path, &        !< path to raw data
         &                             raw_data_ahf_filename, & !< filename AHF raw data
         &                             ahf_buffer_file !< name for AHF buffer file

    INTEGER (KIND=i4)               :: iahf_type, &  !< ID of dataset used !_br 14.04.16
         &                             nuin, & !< unit number
         &                             ierr !< error flag

    !> namelist with filenames for AHF data input
    NAMELIST /ahf_raw_data/ raw_data_ahf_path, raw_data_ahf_filename, iahf_type !_br 14.04.16
    !> namelist with filenames for AHF data output
    NAMELIST /ahf_io_extpar/ ahf_buffer_file
    

    nuin = free_un()  ! functioin free_un returns free Fortran unit number
    OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF
    
    READ(nuin, NML=ahf_raw_data, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist ahf_raw_data',__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=ahf_io_extpar, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist ahf_io_extpar',__FILE__, __LINE__) 
    ENDIF
    
    CLOSE(nuin)
  
  END SUBROUTINE read_namelists_extpar_ahf

  SUBROUTINE read_namelists_extpar_isa(namelist_file,  &
       &                               isa_type,    &
       &                               raw_data_isa_path, &
       &                               raw_data_isa_filename, &
       &                               ntiles_isa, &
       &                               isa_buffer_file      )

  
    CHARACTER (len=*), INTENT(IN) :: namelist_file !< filename with namelists for for EXTPAR settings

    CHARACTER (len=filename_max), INTENT(OUT) :: raw_data_isa_path, &        !< path to raw data
         &                                       raw_data_isa_filename(1:max_tiles_isa), & !< filename isa raw data
         &                                       isa_buffer_file !< name for isa buffer file

    INTEGER, INTENT(OUT)                      :: ntiles_isa

    INTEGER (KIND=i4)                         :: isa_type, &  !< ID of dataset used !_br 14.04.16
         &                                       nuin, & !< unit number
         &                                       ierr !< error flag

    ! NAMELIST /isa_raw_data/ raw_data_isa_path, raw_data_isa_filename, i_isa_data, ilookup_table_isa, ntiles_isa
    NAMELIST /isa_raw_data/ raw_data_isa_path, raw_data_isa_filename, ntiles_isa, isa_type !_br 14.04.16
    !> namelist with filenames for isa data output
    NAMELIST /isa_io_extpar/ isa_buffer_file

    ntiles_isa=1
    nuin = free_un()  ! functioin free_un returns free Fortran unit number
    
    OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=isa_raw_data, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist isa_raw_data',__FILE__, __LINE__) 
    ENDIF
    
    READ(nuin, NML=isa_io_extpar, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist isa_io_extpar',__FILE__, __LINE__) 
    ENDIF
     
    CLOSE(nuin)

  END SUBROUTINE read_namelists_extpar_isa

  !> subroutine to read namelist for hhs ksat data settings for EXTPAR 
  SUBROUTINE read_namelists_extpar_hhs(namelist_file, &
                                      raw_data_hihydrosoil_path, &
                                      raw_data_hihydrosoil_filename_file, &
                                      hihydrosoil_buffer_file)
  
    CHARACTER (len=1024), INTENT(IN) :: namelist_file !< filename with namelists for for EXTPAR settings

    CHARACTER (len=1024)            :: raw_data_hihydrosoil_path, &        !< path to raw data
         &                             raw_data_hihydrosoil_filename_file, & !< filename HHS_KSAT raw data
         &                             hihydrosoil_buffer_file !< name for HHS_KSAT buffer file

    INTEGER (KIND=i4)               :: nuin, & !< unit number
         &                             ierr !< error flag

    !> namelist with filenames for HHS_KSAT data input
    NAMELIST /HiHydroSoil_raw_data/ raw_data_hihydrosoil_path, raw_data_hihydrosoil_filename_file !_br 14.04.16
    !> namelist with filenames for HHS_KSAT data output
    NAMELIST /HiHydroSoil_io_extpar/ hihydrosoil_buffer_file
    
    
    nuin = free_un()  ! functioin free_un returns free Fortran unit number
    OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF
    
    READ(nuin, NML=HiHydroSoil_raw_data, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist hhs_raw_data',__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=HiHydroSoil_io_extpar, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist hhs__io_extpar',__FILE__, __LINE__) 
    ENDIF
    
    CLOSE(nuin)
  
  END SUBROUTINE read_namelists_extpar_hhs

  !> subroutine to read namelist for hhs ormc data settings for EXTPAR 
  SUBROUTINE read_namelists_extpar_hhs_ormc(namelist_file, &
                                      raw_data_hhs_ormc_path, &
                                      raw_data_hhs_ormc_filename, &
                                      hhs_ormc_buffer_file)
  
    CHARACTER (len=1024), INTENT(IN) :: namelist_file !< filename with namelists for for EXTPAR settings

    CHARACTER (len=1024)            :: raw_data_hhs_ormc_path, &        !< path to raw data
         &                             raw_data_hhs_ormc_filename, & !< filename HHS_ORMC raw data
         &                             hhs_ormc_buffer_file !< name for HHS_ORMC buffer file

    INTEGER (KIND=i4)               :: nuin, & !< unit number
         &                             ierr !< error flag

    !> namelist with filenames for HHS_ORMC data input
    NAMELIST /hhs_ormc_raw_data/ raw_data_hhs_ormc_path, raw_data_hhs_ormc_filename !_br 14.04.16
    !> namelist with filenames for HHS_ORMC data output
    NAMELIST /hhs_ormc_io_extpar/ hhs_ormc_buffer_file
    

    nuin = free_un()  ! functioin free_un returns free Fortran unit number
    OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF
    
    READ(nuin, NML=hhs_ormc_raw_data, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist hhs_ormc_raw_data',__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=hhs_ormc_io_extpar, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist hhs_ormc_io_extpar',__FILE__, __LINE__) 
    ENDIF
    
    CLOSE(nuin)
  
  END SUBROUTINE read_namelists_extpar_hhs_ormc

    !> subroutine to read namelist for hhs alfa data settings for EXTPAR 
  SUBROUTINE read_namelists_extpar_hhs_alfa(namelist_file, &
                                      raw_data_hhs_alfa_path, &
                                      raw_data_hhs_alfa_filename, &
                                      hhs_alfa_buffer_file)
  
    CHARACTER (len=1024), INTENT(IN) :: namelist_file !< filename with namelists for for EXTPAR settings

    CHARACTER (len=1024)            :: raw_data_hhs_alfa_path, &        !< path to raw data
         &                             raw_data_hhs_alfa_filename, & !< filename HHS_ALFA raw data
         &                             hhs_alfa_buffer_file !< name for HHS_ALFA buffer file

    INTEGER (KIND=i4)               :: nuin, & !< unit number
         &                             ierr !< error flag

    !> namelist with filenames for HHS_ALFA data input
    NAMELIST /hhs_alfa_raw_data/ raw_data_hhs_alfa_path, raw_data_hhs_alfa_filename !_br 14.04.16
    !> namelist with filenames for HHS_ALFA data output
    NAMELIST /hhs_alfa_io_extpar/ hhs_alfa_buffer_file
    

    nuin = free_un()  ! functioin free_un returns free Fortran unit number
    OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF
    
    READ(nuin, NML=hhs_alfa_raw_data, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist hhs_alfa_raw_data',__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=hhs_alfa_io_extpar, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist hhs_alfa_io_extpar',__FILE__, __LINE__) 
    ENDIF
    
    CLOSE(nuin)
  
  END SUBROUTINE read_namelists_extpar_hhs_alfa

    !> subroutine to read namelist for hhs critw data settings for EXTPAR 
  SUBROUTINE read_namelists_extpar_hhs_critw(namelist_file, &
                                      raw_data_hhs_critw_path, &
                                      raw_data_hhs_critw_filename, &
                                      hhs_critw_buffer_file)
  
    CHARACTER (len=1024), INTENT(IN) :: namelist_file !< filename with namelists for for EXTPAR settings

    CHARACTER (len=1024)            :: raw_data_hhs_critw_path, &        !< path to raw data
         &                             raw_data_hhs_critw_filename, & !< filename HHS_CRITW raw data
         &                             hhs_critw_buffer_file !< name for HHS_CRITW buffer file

    INTEGER (KIND=i4)               :: nuin, & !< unit number
         &                             ierr !< error flag

    !> namelist with filenames for HHS_CRITW data input
    NAMELIST /hhs_critw_raw_data/ raw_data_hhs_critw_path, raw_data_hhs_critw_filename !_br 14.04.16
    !> namelist with filenames for HHS_CRITW data output
    NAMELIST /hhs_critw_io_extpar/ hhs_critw_buffer_file
    

    nuin = free_un()  ! functioin free_un returns free Fortran unit number
    OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF
    
    READ(nuin, NML=hhs_critw_raw_data, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist hhs_critw_raw_data',__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=hhs_critw_io_extpar, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist hhs_critw_io_extpar',__FILE__, __LINE__) 
    ENDIF
    
    CLOSE(nuin)
  
  END SUBROUTINE read_namelists_extpar_hhs_critw

    !> subroutine to read namelist for hhs fieldc data settings for EXTPAR 
  SUBROUTINE read_namelists_extpar_hhs_fieldc(namelist_file, &
                                      raw_data_hhs_fieldc_path, &
                                      raw_data_hhs_fieldc_filename, &
                                      hhs_fieldc_buffer_file)
  
    CHARACTER (len=1024), INTENT(IN) :: namelist_file !< filename with namelists for for EXTPAR settings

    CHARACTER (len=1024)            :: raw_data_hhs_fieldc_path, &        !< path to raw data
         &                             raw_data_hhs_fieldc_filename, & !< filename HHS_FIELDC raw data
         &                             hhs_fieldc_buffer_file !< name for HHS_FIELDC buffer file

    INTEGER (KIND=i4)               :: nuin, & !< unit number
         &                             ierr !< error flag

    !> namelist with filenames for HHS_FIELDC data input
    NAMELIST /hhs_fieldc_raw_data/ raw_data_hhs_fieldc_path, raw_data_hhs_fieldc_filename !_br 14.04.16
    !> namelist with filenames for HHS_FIELDC data output
    NAMELIST /hhs_fieldc_io_extpar/ hhs_fieldc_buffer_file
    

    nuin = free_un()  ! functioin free_un returns free Fortran unit number
    OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF
    
    READ(nuin, NML=hhs_fieldc_raw_data, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist hhs_fieldc_raw_data',__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=hhs_fieldc_io_extpar, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist hhs_fieldc_io_extpar',__FILE__, __LINE__) 
    ENDIF
    
    CLOSE(nuin)
  
  END SUBROUTINE read_namelists_extpar_hhs_fieldc

    !> subroutine to read namelist for hhs n data settings for EXTPAR 
  SUBROUTINE read_namelists_extpar_hhs_n(namelist_file, &
                                      raw_data_hhs_n_path, &
                                      raw_data_hhs_n_filename, &
                                      hhs_n_buffer_file)
  
    CHARACTER (len=1024), INTENT(IN) :: namelist_file !< filename with namelists for for EXTPAR settings

    CHARACTER (len=1024)            :: raw_data_hhs_n_path, &        !< path to raw data
         &                             raw_data_hhs_n_filename, & !< filename HHS_N raw data
         &                             hhs_n_buffer_file !< name for HHS_N buffer file

    INTEGER (KIND=i4)               :: nuin, & !< unit number
         &                             ierr !< error flag

    !> namelist with filenames for HHS_N data input
    NAMELIST /hhs_n_raw_data/ raw_data_hhs_n_path, raw_data_hhs_n_filename !_br 14.04.16
    !> namelist with filenames for HHS_N data output
    NAMELIST /hhs_n_io_extpar/ hhs_n_buffer_file
    

    nuin = free_un()  ! functioin free_un returns free Fortran unit number
    OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF
    
    READ(nuin, NML=hhs_n_raw_data, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist hhs_n_raw_data',__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=hhs_n_io_extpar, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist hhs_n_io_extpar',__FILE__, __LINE__) 
    ENDIF
    
    CLOSE(nuin)
  
  END SUBROUTINE read_namelists_extpar_hhs_n

    !> subroutine to read namelist for hhs satf data settings for EXTPAR 
  SUBROUTINE read_namelists_extpar_hhs_satf(namelist_file, &
                                      raw_data_hhs_satf_path, &
                                      raw_data_hhs_satf_filename, &
                                      hhs_satf_buffer_file)
  
    CHARACTER (len=1024), INTENT(IN) :: namelist_file !< filename with namelists for for EXTPAR settings

    CHARACTER (len=1024)            :: raw_data_hhs_satf_path, &        !< path to raw data
         &                             raw_data_hhs_satf_filename, & !< filename HHS_SATF raw data
         &                             hhs_satf_buffer_file !< name for HHS_SATF buffer file

    INTEGER (KIND=i4)               :: nuin, & !< unit number
         &                             ierr !< error flag

    !> namelist with filenames for HHS_SATF data input
    NAMELIST /hhs_satf_raw_data/ raw_data_hhs_satf_path, raw_data_hhs_satf_filename !_br 14.04.16
    !> namelist with filenames for HHS_SATF data output
    NAMELIST /hhs_satf_io_extpar/ hhs_satf_buffer_file
    

    nuin = free_un()  ! functioin free_un returns free Fortran unit number
    OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF
    
    READ(nuin, NML=hhs_satf_raw_data, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist hhs_satf_raw_data',__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=hhs_satf_io_extpar, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist hhs_satf_io_extpar',__FILE__, __LINE__) 
    ENDIF
    
    CLOSE(nuin)
  
  END SUBROUTINE read_namelists_extpar_hhs_satf

    !> subroutine to read namelist for hhs stc data settings for EXTPAR 
  SUBROUTINE read_namelists_extpar_hhs_stc(namelist_file, &
                                      raw_data_hhs_stc_path, &
                                      raw_data_hhs_stc_filename, &
                                      hhs_stc_buffer_file)
  
    CHARACTER (len=1024), INTENT(IN) :: namelist_file !< filename with namelists for for EXTPAR settings

    CHARACTER (len=1024)            :: raw_data_hhs_stc_path, &        !< path to raw data
         &                             raw_data_hhs_stc_filename, & !< filename HHS_STC raw data
         &                             hhs_stc_buffer_file !< name for HHS_STC buffer file

    INTEGER (KIND=i4)               :: nuin, & !< unit number
         &                             ierr !< error flag

    !> namelist with filenames for HHS_STC data input
    NAMELIST /hhs_stc_raw_data/ raw_data_hhs_stc_path, raw_data_hhs_stc_filename !_br 14.04.16
    !> namelist with filenames for HHS_STC data output
    NAMELIST /hhs_stc_io_extpar/ hhs_stc_buffer_file
    

    nuin = free_un()  ! functioin free_un returns free Fortran unit number
    OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF
    
    READ(nuin, NML=hhs_stc_raw_data, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist hhs_stc_raw_data',__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=hhs_stc_io_extpar, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist hhs_stc_io_extpar',__FILE__, __LINE__) 
    ENDIF
    
    CLOSE(nuin)
  
  END SUBROUTINE read_namelists_extpar_hhs_stc

    !> subroutine to read namelist for hhs wcav data settings for EXTPAR 
  SUBROUTINE read_namelists_extpar_hhs_wcav(namelist_file, &
                                      raw_data_hhs_wcav_path, &
                                      raw_data_hhs_wcav_filename, &
                                      hhs_wcav_buffer_file)
  
    CHARACTER (len=1024), INTENT(IN) :: namelist_file !< filename with namelists for for EXTPAR settings

    CHARACTER (len=1024)            :: raw_data_hhs_wcav_path, &        !< path to raw data
         &                             raw_data_hhs_wcav_filename, & !< filename HHS_WCAV raw data
         &                             hhs_wcav_buffer_file !< name for HHS_WCAV buffer file

    INTEGER (KIND=i4)               :: nuin, & !< unit number
         &                             ierr !< error flag

    !> namelist with filenames for HHS_WCAV data input
    NAMELIST /hhs_wcav_raw_data/ raw_data_hhs_wcav_path, raw_data_hhs_wcav_filename !_br 14.04.16
    !> namelist with filenames for HHS_WCAV data output
    NAMELIST /hhs_wcav_io_extpar/ hhs_wcav_buffer_file
    

    nuin = free_un()  ! functioin free_un returns free Fortran unit number
    OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF
    
    READ(nuin, NML=hhs_wcav_raw_data, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist hhs_wcav_raw_data',__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=hhs_wcav_io_extpar, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist hhs_wcav_io_extpar',__FILE__, __LINE__) 
    ENDIF
    
    CLOSE(nuin)
  
  END SUBROUTINE read_namelists_extpar_hhs_wcav

    !> subroutine to read namelist for hhs wcpf2 data settings for EXTPAR 
  SUBROUTINE read_namelists_extpar_hhs_wcpf2(namelist_file, &
                                      raw_data_hhs_wcpf2_path, &
                                      raw_data_hhs_wcpf2_filename, &
                                      hhs_wcpf2_buffer_file)
  
    CHARACTER (len=1024), INTENT(IN) :: namelist_file !< filename with namelists for for EXTPAR settings

    CHARACTER (len=1024)            :: raw_data_hhs_wcpf2_path, &        !< path to raw data
         &                             raw_data_hhs_wcpf2_filename, & !< filename HHS_WCPF2 raw data
         &                             hhs_wcpf2_buffer_file !< name for HHS_WCPF2 buffer file

    INTEGER (KIND=i4)               :: nuin, & !< unit number
         &                             ierr !< error flag

    !> namelist with filenames for HHS_WCPF2 data input
    NAMELIST /hhs_wcpf2_raw_data/ raw_data_hhs_wcpf2_path, raw_data_hhs_wcpf2_filename !_br 14.04.16
    !> namelist with filenames for HHS_WCPF2 data output
    NAMELIST /hhs_wcpf2_io_extpar/ hhs_wcpf2_buffer_file
    

    nuin = free_un()  ! functioin free_un returns free Fortran unit number
    OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF
    
    READ(nuin, NML=hhs_wcpf2_raw_data, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist hhs_wcpf2_raw_data',__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=hhs_wcpf2_io_extpar, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist hhs_wcpf2_io_extpar',__FILE__, __LINE__) 
    ENDIF
    
    CLOSE(nuin)
  
  END SUBROUTINE read_namelists_extpar_hhs_wcpf2

    !> subroutine to read namelist for hhs wcpf3 data settings for EXTPAR 
  SUBROUTINE read_namelists_extpar_hhs_wcpf3(namelist_file, &
                                      raw_data_hhs_wcpf3_path, &
                                      raw_data_hhs_wcpf3_filename, &
                                      hhs_wcpf3_buffer_file)
  
    CHARACTER (len=1024), INTENT(IN) :: namelist_file !< filename with namelists for for EXTPAR settings

    CHARACTER (len=1024)            :: raw_data_hhs_wcpf3_path, &        !< path to raw data
         &                             raw_data_hhs_wcpf3_filename, & !< filename HHS_WCPF3 raw data
         &                             hhs_wcpf3_buffer_file !< name for HHS_WCPF3 buffer file

    INTEGER (KIND=i4)               :: nuin, & !< unit number
         &                             ierr !< error flag

    !> namelist with filenames for HHS_WCPF3 data input
    NAMELIST /hhs_wcpf3_raw_data/ raw_data_hhs_wcpf3_path, raw_data_hhs_wcpf3_filename !_br 14.04.16
    !> namelist with filenames for HHS_WCPF3 data output
    NAMELIST /hhs_wcpf3_io_extpar/ hhs_wcpf3_buffer_file
    

    nuin = free_un()  ! functioin free_un returns free Fortran unit number
    OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF
    
    READ(nuin, NML=hhs_wcpf3_raw_data, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist hhs_wcpf3_raw_data',__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=hhs_wcpf3_io_extpar, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist hhs_wcpf3_io_extpar',__FILE__, __LINE__) 
    ENDIF
    
    CLOSE(nuin)
  
  END SUBROUTINE read_namelists_extpar_hhs_wcpf3

    !> subroutine to read namelist for hhs wcpf42 data settings for EXTPAR 
  SUBROUTINE read_namelists_extpar_hhs_wcpf42(namelist_file, &
                                      raw_data_hhs_wcpf42_path, &
                                      raw_data_hhs_wcpf42_filename, &
                                      hhs_wcpf42_buffer_file)
  
    CHARACTER (len=1024), INTENT(IN) :: namelist_file !< filename with namelists for for EXTPAR settings

    CHARACTER (len=1024)            :: raw_data_hhs_wcpf42_path, &        !< path to raw data
         &                             raw_data_hhs_wcpf42_filename, & !< filename HHS_WCPF42 raw data
         &                             hhs_wcpf42_buffer_file !< name for HHS_WCPF42 buffer file

    INTEGER (KIND=i4)               :: nuin, & !< unit number
         &                             ierr !< error flag

    !> namelist with filenames for HHS_WCPF42 data input
    NAMELIST /hhs_wcpf42_raw_data/ raw_data_hhs_wcpf42_path, raw_data_hhs_wcpf42_filename !_br 14.04.16
    !> namelist with filenames for HHS_WCPF42 data output
    NAMELIST /hhs_wcpf42_io_extpar/ hhs_wcpf42_buffer_file
    

    nuin = free_un()  ! functioin free_un returns free Fortran unit number
    OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF
    
    READ(nuin, NML=hhs_wcpf42_raw_data, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist hhs_wcpf42_raw_data',__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=hhs_wcpf42_io_extpar, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist hhs_wcpf42_io_extpar',__FILE__, __LINE__) 
    ENDIF
    
    CLOSE(nuin)
  
  END SUBROUTINE read_namelists_extpar_hhs_wcpf42

    !> subroutine to read namelist for hhs wcres data settings for EXTPAR 
  SUBROUTINE read_namelists_extpar_hhs_wcres(namelist_file, &
                                      raw_data_hhs_wcres_path, &
                                      raw_data_hhs_wcres_filename, &
                                      hhs_wcres_buffer_file)
  
    CHARACTER (len=1024), INTENT(IN) :: namelist_file !< filename with namelists for for EXTPAR settings

    CHARACTER (len=1024)            :: raw_data_hhs_wcres_path, &        !< path to raw data
         &                             raw_data_hhs_wcres_filename, & !< filename HHS_WCRES raw data
         &                             hhs_wcres_buffer_file !< name for HHS_WCRES buffer file

    INTEGER (KIND=i4)               :: nuin, & !< unit number
         &                             ierr !< error flag

    !> namelist with filenames for HHS_WCRES data input
    NAMELIST /hhs_wcres_raw_data/ raw_data_hhs_wcres_path, raw_data_hhs_wcres_filename !_br 14.04.16
    !> namelist with filenames for HHS_WCRES data output
    NAMELIST /hhs_wcres_io_extpar/ hhs_wcres_buffer_file
    

    nuin = free_un()  ! functioin free_un returns free Fortran unit number
    OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF
    
    READ(nuin, NML=hhs_wcres_raw_data, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist hhs_wcres_raw_data',__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=hhs_wcres_io_extpar, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist hhs_wcres_io_extpar',__FILE__, __LINE__) 
    ENDIF
    
    CLOSE(nuin)
  
  END SUBROUTINE read_namelists_extpar_hhs_wcres

    !> subroutine to read namelist for hhs wcsat data settings for EXTPAR 
  SUBROUTINE read_namelists_extpar_hhs_wcsat(namelist_file, &
                                      raw_data_hhs_wcsat_path, &
                                      raw_data_hhs_wcsat_filename, &
                                      hhs_wcsat_buffer_file)
  
    CHARACTER (len=1024), INTENT(IN) :: namelist_file !< filename with namelists for for EXTPAR settings

    CHARACTER (len=1024)            :: raw_data_hhs_wcsat_path, &        !< path to raw data
         &                             raw_data_hhs_wcsat_filename, & !< filename HHS_WCSAT raw data
         &                             hhs_wcsat_buffer_file !< name for HHS_WCSAT buffer file

    INTEGER (KIND=i4)               :: nuin, & !< unit number
         &                             ierr !< error flag

    !> namelist with filenames for HHS_WCSAT data input
    NAMELIST /hhs_wcsat_raw_data/ raw_data_hhs_wcsat_path, raw_data_hhs_wcsat_filename !_br 14.04.16
    !> namelist with filenames for HHS_WCSAT data output
    NAMELIST /hhs_wcsat_io_extpar/ hhs_wcsat_buffer_file
    

    nuin = free_un()  ! functioin free_un returns free Fortran unit number
    OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF
    
    READ(nuin, NML=hhs_wcsat_raw_data, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist hhs_wcsat_raw_data',__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=hhs_wcsat_io_extpar, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist hhs_wcsat_io_extpar',__FILE__, __LINE__) 
    ENDIF
    
    CLOSE(nuin)
  
  END SUBROUTINE read_namelists_extpar_hhs_wcsat
  
  
  !> open netcdf-file and get netcdf unit file number
  SUBROUTINE open_netcdf_ALB_data(path_alb_file, &
                                          ncid)
    CHARACTER (len=*), INTENT(IN) :: path_alb_file     !< filename with path for albedo raw data
    INTEGER, INTENT(OUT)          :: ncid                       !< netcdf unit file number
    !! open netcdf file 
    CALL check_netcdf( nf90_open(TRIM(path_alb_file),NF90_NOWRITE, ncid))

  END SUBROUTINE open_netcdf_ALB_data

  SUBROUTINE const_check_interpol_alb(alb_field_mom_d,fr_land_lu,alb_min)


    REAL(KIND=wp), INTENT(INOUT) :: alb_field_mom_d(:,:,:,:)

    REAL(KIND=wp), INTENT(IN)    :: fr_land_lu(:,:,:), &
         &                          alb_min

    INTEGER (KIND=i4), PARAMETER :: mpy=12     !< month per year
    
    INTEGER (KIND=i4)            :: i_miss, &
         &                          t,k,j,i,i2,j2, &
         &                          igrid_type 

    REAL (KIND=wp)               :: lon_geo_w,lon_geo_e,lat_geo_n,lat_geo_s, &
         &                          alb_sw,alb_nw,alb_se,alb_ne, &
         &                          bwlon,bwlat

    igrid_type = tg%igrid_type

    !preparing albedo values for interpolation, original data is copied to
    !alb_interpol to easier remove the procedure if necessary.
    !at the end alb_interpol is written in alb_field_mom

    DO t=1, mpy
      DO k=1,tg%ke
        DO j=1,tg%je
          DO i=1,tg%ie
            alb_interpol(i,j,k,t) = alb_field_mom_d(i,j,k,t)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO t=1, mpy
      DO k=1,tg%ke
        DO j=1,tg%je
          DO i=1,tg%ie
            IF (igrid_type.eq.2) THEN !COSMO interpolation
              !albedo < 0.07 and no water point
              IF ((alb_interpol(i,j,k,t).lt.alb_min).AND.(fr_land_lu(i,j,k).ge.0.5)) THEN 
                  !avoiding 'out of range' errors at the border
                j2 = j
                i2 = i
                IF (j.EQ.(tg%je)) THEN
                  j2 = j-1
                ENDIF
                IF (j.EQ.1) THEN
                  j2 = j+1
                ENDIF
                IF (i.EQ.(tg%ie)) THEN
                  i2 = i-1
                ENDIF
                IF (i.EQ.1) THEN
                  i2 = i+1
                ENDIF

                !coordinates and albedo values for neighbor points

                lon_geo_w = lon_geo(i2-1,j2,k)
                lon_geo_e = lon_geo(i2+1,j2,k)
                lat_geo_n = lat_geo(i2,j2+1,k)
                lat_geo_s = lat_geo(i2,j2-1,k)

                alb_sw = alb_interpol(i2-1,j2-1,k,t)
                alb_se = alb_interpol(i2+1,j2-1,k,t)
                alb_nw = alb_interpol(i2-1,j2+1,k,t)
                alb_ne = alb_interpol(i2+1,j2+1,k,t)

                !calculate weights for interpolation
                !these weights might give wrong results, if there are no 4 valid
                !points surrounding the interpolated point(i,j)
                CALL calc_weight_bilinear_interpol(lon_geo(i,j,k), &
                                                    lat_geo(i,j,k), &
                                                    lon_geo_w,      &
                                                    lon_geo_e,      &
                                                    lat_geo_n,     &
                                                    lat_geo_s,     &
                                                    bwlon,         &
                                                    bwlat)
                 
                i_miss = 0

                IF (fr_land_lu(i2-1,j2-1,k).LT.0.5) THEN
                  alb_sw = 0.
                  i_miss = i_miss + 1
                ENDIF
                IF (fr_land_lu(i2+1,j2-1,k).LT.0.5) THEN
                  alb_se = 0.
                  i_miss = i_miss + 1
                ENDIF
                IF (fr_land_lu(i2-1,j2+1,k).LT.0.5) THEN
                  alb_nw = 0.
                  i_miss = i_miss + 1
                ENDIF
                IF (fr_land_lu(i2+1,j2+1,k).LT.0.5) THEN
                  alb_ne = 0.
                  i_miss = i_miss + 1
                ENDIF


                IF (i_miss.EQ.4) THEN
                !    if there are no valid interpolation values, the next step is skipped  
                  i_miss = 5
                  GOTO 100
                ENDIF
                 
                alb_interpol(i,j,k,t) = calc_value_bilinear_interpol(bwlon, bwlat, &
                                        alb_sw, alb_se, alb_ne, alb_nw)*(4/(4-i_miss))

    100         IF (alb_interpol(i,j,k,t).LT.alb_min.and. soiltype_fao(i,j,k).LE.9 .AND. soiltype_fao(i,j,k).GE.0) THEN
                  !values that are still too small, will receive a soiltype dependent albedo
                  alb_interpol(i,j,k,t) = zalso(soiltype_fao(i,j,k),t)*fr_land_lu(i,j,k) + &
                                            alb_min*(1.-fr_land_lu(i,j,k))
                ELSEIF (alb_interpol(i,j,k,t).LT.alb_min .AND. (soiltype_fao(i,j,k).GT.9 &
                        & .OR. soiltype_fao(i,j,k).EQ.0)) THEN
                  alb_interpol(i,j,k,t) = alb_min
                ENDIF
  !override wrong values due to bad interpolation (occurs only at some borderpoints)
                IF (alb_interpol(i,j,k,t).GT.0.7) THEN
                  alb_interpol(i,j,k,t) = zalso(soiltype_fao(i,j,k),t)*fr_land_lu(i,j,k) + &
                                            alb_min*(1.-fr_land_lu(i,j,k))
                ENDIF

                !water points:
              ELSEIF ((alb_interpol(i,j,k,t).LT.alb_min).AND.(fr_land_lu(i,j,k).LT.0.5)) THEN
                alb_interpol(i,j,k,t) = alb_min
              ELSEIF ((alb_interpol(i,j,k,t).GE.alb_min).AND.(fr_land_lu(i,j,k).LT.0.5)) THEN
                alb_interpol(i,j,k,t) = alb_min
              ENDIF

            ELSEIF (igrid_type.EQ.1) THEN !ICON interpolation 

              IF (fr_land_lu(i,j,k).LT.0.01) THEN
                !water point
                alb_interpol(i,j,k,t) = 0.07
              ELSEIF (alb_interpol(i,j,k,t).LE.0.) THEN
                !land point with albedo = 0
   
                alb_se = alb_interpol(icon_grid_region%cells%neighbor_index(i,1),j,k,t)
                alb_nw = alb_interpol(icon_grid_region%cells%neighbor_index(i,2),j,k,t)
                alb_ne = alb_interpol(icon_grid_region%cells%neighbor_index(i,3),j,k,t)

                IF (alb_se.GT.0.AND.alb_nw.GT.0.AND.alb_ne.GT.0) THEN
                  alb_interpol(i,j,k,t) = (alb_se+alb_nw+alb_ne)/3
                ELSE IF (alb_se.GT.0.AND.alb_nw.GT.0) THEN
                  alb_interpol(i,j,k,t) = (alb_se+alb_nw)/2
                ELSE IF (alb_nw.GT.0.AND.alb_ne.GT.0) THEN
                  alb_interpol(i,j,k,t) = (alb_nw+alb_ne)/2
                ELSE IF (alb_ne.GT.0.AND.alb_se.GT.0) THEN
                  alb_interpol(i,j,k,t) = (alb_ne+alb_se)/2
                ELSE
                  alb_interpol(i,j,k,t) = zalso(soiltype_fao(i,j,k),t)*fr_land_lu(i,j,k) + &
                                            0.07*(1.-fr_land_lu(i,j,k))
                ENDIF

              ENDIF     
            ENDIF  !ICON interpolation

             !Gletscher
            IF ((soiltype_fao(i,j,k).EQ.1).AND.(fr_land_lu(i,j,k).GE.0.01)) THEN
              alb_interpol(i,j,k,t) = zalso(soiltype_fao(i,j,k),t)
            ENDIF
   
            !catching false values (occurs only at the borders)
            IF (alb_interpol(i,j,k,t).gt.0.7) THEN
              alb_interpol(i,j,k,t) = zalso(soiltype_fao(i,j,k),t)*fr_land_lu(i,j,k) + &
                                     0.07*(1.-fr_land_lu(i,j,k))
            ENDIF

          ENDDO !i
        ENDDO !j
      ENDDO !k
    ENDDO !t           


    DO t=1, mpy
      DO k=1,tg%ke
        DO j=1,tg%je
          DO i=1,tg%ie
            alb_field_mom_d(i,j,k,t) = 100*alb_interpol(i,j,k,t)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
     
  END SUBROUTINE const_check_interpol_alb

END MODULE mo_python_routines
